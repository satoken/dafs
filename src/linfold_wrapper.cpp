/*
 * LinFold wrapper implementation for DAFS
 */

#include "linfold_wrapper.h"
#include "linearfold/param/profile.h"
#include "ribosum.h"
#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <unordered_map>
#include <vector>

namespace
{
template <typename Options>
void allow_canonical_pairs(Options& options)
{
  // _Fold::Options starts with an empty allowed-pair table.  LinFold's
  // constraint construction consults this table even without explicit
  // structure constraints, so leaving it empty reduces the ensemble to the
  // all-unpaired structure.  set_allowed_pair() is symmetric.
  options.set_allowed_pair('a', 'u');
  options.set_allowed_pair('c', 'g');
  options.set_allowed_pair('g', 'u');
}

const char* model_name(LinFoldWrapper::ModelType model_type)
{
  return model_type == LinFoldWrapper::ModelType::LPV ? "LPV" : "LPC";
}

[[noreturn]] void rethrow_linfold_failure(
    LinFoldWrapper::ModelType model_type, size_t sequence_length,
    const char* detail)
{
  throw std::runtime_error(
      std::string("LinearPartition ") + model_name(model_type) +
      " failed for sequence length " + std::to_string(sequence_length) +
      ": " + detail);
}

int pair_type(char a, char b)
{
  const int x = Ribosum85_60::nucleotide_code(a);
  const int y = Ribosum85_60::nucleotide_code(b);
  if (x < 0 || y < 0) return -1;
  if ((x == 0 && y == 3) || (x == 3 && y == 0) ||
      (x == 1 && y == 2) || (x == 2 && y == 1) ||
      (x == 2 && y == 3) || (x == 3 && y == 2))
    return 4 * x + y;
  return -1;
}

double calculate_profile_impl(const std::vector<std::string>& aligned,
                            const std::string& constraint, unsigned beam,
                            float threshold, BP& bp, bool relax_unsupported_pairs,
                            LinFoldWrapper::ProfileEnergyMode profile_energy_mode,
                            std::unique_ptr<LinFold<ProfileNearestNeighbor>>& engine)
{
  using Profile = ProfileNearestNeighbor;
  const size_t length = aligned.front().size();
  const auto energy_mode = profile_energy_mode ==
      LinFoldWrapper::ProfileEnergyMode::RNAalifold
      ? Profile::EnergyMode::RNAalifold
      : Profile::EnergyMode::Legacy;
  auto profile = std::make_unique<Profile>(aligned, energy_mode);
  auto* profile_ptr = profile.get();
  if (engine)
    engine->reset_param_model(std::move(profile));
  else
    engine = std::make_unique<LinFold<Profile>>(std::move(profile));
  typename LinFold<Profile>::Options options;
  options.beam_size(beam);
  options.alphabets("abcdefghijklmnop");
  options.occupancy_prefix_ = profile_ptr->occupancy_prefix();
  options.probability_cutoff_ = threshold;
  options.cache_position_pair_queries_ = profile_ptr->row_count() > 4;
  if (profile_ptr->row_count() == 2) {
    options.pair_signature_.resize(length + 1);
    for (size_t column = 1; column <= length; ++column)
      options.pair_signature_[column] = profile_ptr->two_row_signature(column);
  }

  // A column's four-bit nucleotide mask makes covarying pairs available to
  // the beam search even when the majority bases are incompatible.
  std::string symbols(length, 'a');
  std::vector<uint32_t> previous_occupied(length + 1, 0);
  bool has_all_gap_column = false;
  for (size_t col = 0; col < length; ++col) {
    unsigned mask = 0;
    bool occupied = false;
    for (const auto& row : aligned) {
      const int code = Ribosum85_60::nucleotide_code(row[col]);
      if (code >= 0) mask |= 1u << code;
      occupied |= row[col] != '-';
    }
    symbols[col] = static_cast<char>('a' + mask);
    previous_occupied[col + 1] = occupied
        ? static_cast<uint32_t>(col + 1) : previous_occupied[col];
    has_all_gap_column |= !occupied;
  }
  if (has_all_gap_column)
    options.previous_occupied_ = std::move(previous_occupied);
  for (unsigned first = 1; first < 16; ++first)
    for (unsigned second = 1; second < 16; ++second) {
      bool compatible = false;
      for (unsigned x = 0; x < 4; ++x)
        for (unsigned y = 0; y < 4; ++y)
          if ((first & (1u << x)) && (second & (1u << y))) {
            constexpr char nucleotides[] = "acgu";
            compatible |= pair_type(nucleotides[x], nucleotides[y]) >= 0;
          }
      if (compatible)
        options.set_allowed_pair(static_cast<char>('a' + first),
                                 static_cast<char>('a' + second));
    }
  options.position_pair([profile_ptr](unsigned i, unsigned j) {
    return profile_ptr->can_pair(i, j);
  });

  if (!constraint.empty()) {
    if (constraint.size() != length)
      throw std::invalid_argument("profile constraint length differs from alignment");
    std::vector<uint32_t> parsed(length + 1, _Fold::Options::ANY);
    std::vector<size_t> stack;
    for (size_t i = 0; i < length; ++i) {
      switch (constraint[i]) {
      case '?': break;
      case '.': parsed[i + 1] = _Fold::Options::UNPAIRED; break;
      case '(' : stack.push_back(i); break;
      case ')' : {
        if (stack.empty()) throw std::invalid_argument("unmatched ')' in profile constraint");
        const size_t left = stack.back();
        stack.pop_back();
        if (!options.allow_hairpin(left+1, i+1) || !profile_ptr->can_pair(left+1, i+1)) {
          if (!relax_unsupported_pairs)
            throw std::invalid_argument("infeasible fixed pair in profile constraint");
          // The mixed posterior can predict minority-supported pairs outside
          // this profile ensemble.  As in per-row refinement, remove these
          // pairs only from the model that cannot realize them.
          parsed[left + 1] = parsed[i + 1] = _Fold::Options::UNPAIRED;
        } else {
          parsed[left + 1] = i + 1;
          parsed[i + 1] = left + 1;
        }
        break;
      }
      default: throw std::invalid_argument("invalid profile constraint character");
      }
    }
    if (!stack.empty()) throw std::invalid_argument("unmatched '(' in profile constraint");
    options.constraints(parsed);
  }

  // Match ViennaRNA's default RNAalifold covariance term.  DAFS's
  // RIBOSUM85-60 score belongs to alignment/final-decoding objectives and is
  // intentionally not part of the profile partition function.
  std::unordered_map<uint64_t, float> pair_scores;
  std::array<float, 36 * 36> two_row_pair_scores;
  two_row_pair_scores.fill(std::numeric_limits<float>::quiet_NaN());
  const bool two_rows = profile_ptr->row_count() == 2;
  const size_t pair_score_cache_limit = std::max<size_t>(4096, length * 64);
  options.pair_score([&](unsigned i, unsigned j) {
    const size_t signature_pair = two_rows
        ? static_cast<size_t>(options.pair_signature_[i]) * 36 +
              options.pair_signature_[j]
        : 0;
    const uint64_t key = (static_cast<uint64_t>(i) << 32) | j;
    if (two_rows) {
      const float score = two_row_pair_scores[signature_pair];
      if (!std::isnan(score)) return score;
    } else {
      const auto found = pair_scores.find(key);
      if (found != pair_scores.end()) return found->second;
    }
    const auto counts = profile_ptr->pair_type_counts(i, j);
    const unsigned n = static_cast<unsigned>(aligned.size());
    // The default ViennaRNA matrix for canonical pair-type substitutions.
    static constexpr unsigned distance[7][7] = {
      {0, 0, 0, 0, 0, 0, 0},
      {0, 0, 2, 2, 1, 2, 2},
      {0, 2, 0, 1, 2, 2, 2},
      {0, 2, 1, 0, 2, 1, 2},
      {0, 1, 2, 2, 0, 2, 1},
      {0, 2, 2, 1, 2, 0, 2},
      {0, 2, 2, 2, 1, 2, 0},
    };
    double covariance = 0;
    for (unsigned a = 1; a <= 6; ++a)
      for (unsigned b = a; b <= 6; ++b)
        covariance += static_cast<double>(counts[a]) * counts[b] * distance[a][b];
    // ViennaRNA stores covariance in centikcal/mol.  LinFold scores are
    // dimensionless log weights, hence divide by 100*RT.
    const double incompatible = static_cast<double>(counts[0]);
    const double gap_gap = static_cast<double>(counts[7]);
    const double centikcal = 100.0 * covariance / n -
                             100.0 * (incompatible + 0.25 * gap_gap);
    // RNAalifold's comparative PF grammar applies half of the centikcal
    // covariance score to each paired transition.  Keep that scale when
    // converting to LinFold's dimensionless log weight.
    const double score = centikcal / (200.0 * ProfileNearestNeighbor::RT);
    const float score_as_float = static_cast<float>(score);
    if (two_rows)
      two_row_pair_scores[signature_pair] = score_as_float;
    else if (pair_scores.size() < pair_score_cache_limit)
      pair_scores.emplace(key, score_as_float);
    return score_as_float;
  });

  const auto log_partition = engine->compute_inside(symbols, options);
  if (!std::isfinite(log_partition))
    throw std::runtime_error("profile partition function is not finite");
  engine->compute_outside(symbols, options);
  const auto bpp = engine->compute_basepairing_probabilities(symbols, options);
  if (bpp.size() != length + 1)
    throw std::runtime_error("unexpected profile BPP dimensions");
  bp.assign(length, {});
  for (size_t i = 1; i <= length; ++i)
    for (const auto& [j, probability] : bpp[i])
      if (j > i && j <= length && probability > threshold)
        bp[i - 1].emplace_back(j - 1, probability);
  return log_partition;
}
}

LinFoldWrapper::LinFoldWrapper(float th, ModelType model_type, uint32_t beam_size,
                               ProfileEnergyMode profile_energy_mode)
  : Fold::Model(th), beam_size_(beam_size), model_type_(model_type),
    profile_energy_mode_(profile_energy_mode),
    linfold_turner_(nullptr), linfold_contra_(nullptr)
{
}

double LinFoldWrapper::calculate_profile(
    const ALN& alignment, const std::vector<Fasta>& sequences,
    BP& bp, const std::string& constraint, bool relax_unsupported_pairs) const
{
  if (alignment.empty())
    throw std::invalid_argument("empty profile alignment");
  const size_t length = alignment.front().second.size();
  std::vector<std::string> aligned;
  aligned.reserve(alignment.size());
  for (const auto& [index, present] : alignment) {
    if (index >= sequences.size() || present.size() != length)
      throw std::invalid_argument("invalid profile alignment");
    const auto& sequence = sequences[index].seq();
    std::string row(length, '-');
    size_t position = 0;
    for (size_t col = 0; col < length; ++col)
      if (present[col]) {
        if (position >= sequence.size())
          throw std::invalid_argument("profile alignment exceeds sequence");
        row[col] = sequence[position++];
      }
    if (position != sequence.size())
      throw std::invalid_argument("profile alignment does not consume sequence");
    aligned.push_back(std::move(row));
  }
  try {
    std::lock_guard<std::mutex> lock(profile_mutex_);
    return calculate_profile_impl(aligned, constraint, beam_size_, threshold(), bp,
                                  relax_unsupported_pairs, profile_energy_mode_,
                                  linfold_profile_);
  } catch (const std::exception& e) {
    throw std::runtime_error(std::string("LinearAlifold profile folding failed: ") + e.what());
  }
}

void 
LinFoldWrapper::calculate(const std::string& seq, BP& bp)
{
  try {
    std::vector<std::vector<std::pair<u_int32_t, float>>> bpp;
    
    if (model_type_ == ModelType::LPV) {
      // Initialize LinFold with TurnerNearestNeighbor parameters for this sequence
      auto turner_params = std::make_unique<TurnerNearestNeighbor>(seq);
      if (linfold_turner_)
        linfold_turner_->reset_param_model(std::move(turner_params));
      else
        linfold_turner_ = std::make_unique<LinFold<TurnerNearestNeighbor>>(std::move(turner_params));
      
      // Set up options with beam size
      LinFold<TurnerNearestNeighbor>::Options opt;
      opt.beam_size(beam_size_);
      allow_canonical_pairs(opt);
      
      // First compute inside algorithm (forward pass)
      linfold_turner_->compute_inside(seq, opt);
      
      // Then compute outside algorithm (backward pass)
      linfold_turner_->compute_outside(seq, opt);
      
      // Now compute base-pairing probabilities
      bpp = linfold_turner_->compute_basepairing_probabilities(seq, opt);
    } else { // ModelType::LPC
      // Initialize LinFold with CONTRAfoldNearestNeighbor parameters for this sequence
      auto contra_params = std::make_unique<CONTRAfoldNearestNeighbor>(seq);
      if (linfold_contra_)
        linfold_contra_->reset_param_model(std::move(contra_params));
      else
        linfold_contra_ = std::make_unique<LinFold<CONTRAfoldNearestNeighbor>>(std::move(contra_params));
      
      // Set up options with beam size
      LinFold<CONTRAfoldNearestNeighbor>::Options opt;
      opt.beam_size(beam_size_);
      allow_canonical_pairs(opt);
      
      // First compute inside algorithm (forward pass)
      linfold_contra_->compute_inside(seq, opt);
      
      // Then compute outside algorithm (backward pass)
      linfold_contra_->compute_outside(seq, opt);
      
      // Now compute base-pairing probabilities
      bpp = linfold_contra_->compute_basepairing_probabilities(seq, opt);
    }
    
    // Convert LinFold output format to BP format used by DAFS
    uint L = seq.size();
    bp.clear();
    bp.resize(L);
    
    // LinFold returns a 1-based (L+1)-row matrix and 1-based partner
    // coordinates; DAFS BP is 0-based with exactly L rows.
    if (bpp.size() != L + 1)
      throw std::runtime_error("unexpected base-pair probability dimensions");
    for (uint i = 1; i <= L; ++i)
    {
      for (const auto& pair : bpp[i])
      {
        uint j = pair.first;
        float prob = pair.second;
        
        // Only store probabilities above threshold
        if (prob > threshold() && j > i && j <= L)
        {
          bp[i - 1].push_back(std::make_pair(j - 1, prob));
        }
      }
    }
  } catch (const std::exception& e) {
    rethrow_linfold_failure(model_type_, seq.size(), e.what());
  } catch (...) {
    rethrow_linfold_failure(model_type_, seq.size(), "unknown exception");
  }
}

void 
LinFoldWrapper::calculate(const std::string& seq, const std::string& str, BP& bp)
{
  try {
    if (str.size() != seq.size())
      throw std::invalid_argument(
          "structure constraint length differs from sequence length");

    std::vector<std::vector<std::pair<u_int32_t, float>>> bpp;
    
    // Convert structure constraint to LinFold format
    // In DAFS: '.' = unpaired, '(' ')' = paired, '?' = any
    // In LinFold constraint format, we need to set up the constraint structure
    std::vector<uint32_t> constraint(str.size() + 1, _Fold::Options::ANY);
    
    // Parse the structure string to build constraints
    std::vector<int> stack;
    for (size_t i = 0; i < str.size(); ++i)
    {
      if (str[i] == '(')
      {
        stack.push_back(i);
      }
      else if (str[i] == ')')
      {
        if (stack.empty())
          throw std::invalid_argument("unmatched ')' in structure constraint");
        int j = stack.back();
        stack.pop_back();
        // LinFold uses 1-based indexing for constraints
        constraint[j + 1] = i + 1;
        constraint[i + 1] = j + 1;
      }
      else if (str[i] == '.')
      {
        constraint[i + 1] = _Fold::Options::UNPAIRED;
      }
      // '?' remains as ANY
    }
    if (!stack.empty())
      throw std::invalid_argument("unmatched '(' in structure constraint");
    
    if (model_type_ == ModelType::LPV) {
      // Initialize LinFold with TurnerNearestNeighbor parameters for this sequence
      auto turner_params = std::make_unique<TurnerNearestNeighbor>(seq);
      if (linfold_turner_)
        linfold_turner_->reset_param_model(std::move(turner_params));
      else
        linfold_turner_ = std::make_unique<LinFold<TurnerNearestNeighbor>>(std::move(turner_params));
      
      // Set up options with beam size and constraint
      LinFold<TurnerNearestNeighbor>::Options opt;
      opt.beam_size(beam_size_);
      opt.constraints(constraint);
      allow_canonical_pairs(opt);
      
      // First compute inside algorithm (forward pass)
      linfold_turner_->compute_inside(seq, opt);
      
      // Then compute outside algorithm (backward pass)
      linfold_turner_->compute_outside(seq, opt);
      
      // Now compute base-pairing probabilities with constraints
      bpp = linfold_turner_->compute_basepairing_probabilities(seq, opt);
    } else { // ModelType::LPC
      // Initialize LinFold with CONTRAfoldNearestNeighbor parameters for this sequence
      auto contra_params = std::make_unique<CONTRAfoldNearestNeighbor>(seq);
      if (linfold_contra_)
        linfold_contra_->reset_param_model(std::move(contra_params));
      else
        linfold_contra_ = std::make_unique<LinFold<CONTRAfoldNearestNeighbor>>(std::move(contra_params));
      
      // Set up options with beam size and constraint
      LinFold<CONTRAfoldNearestNeighbor>::Options opt;
      opt.beam_size(beam_size_);
      opt.constraints(constraint);
      allow_canonical_pairs(opt);
      
      // First compute inside algorithm (forward pass)
      linfold_contra_->compute_inside(seq, opt);
      
      // Then compute outside algorithm (backward pass)
      linfold_contra_->compute_outside(seq, opt);
      
      // Now compute base-pairing probabilities with constraints
      bpp = linfold_contra_->compute_basepairing_probabilities(seq, opt);
    }
    
    // Convert LinFold output format to BP format used by DAFS
    uint L = seq.size();
    bp.clear();
    bp.resize(L);
    
    if (bpp.size() != L + 1)
      throw std::runtime_error("unexpected base-pair probability dimensions");
    for (uint i = 1; i <= L; ++i)
    {
      for (const auto& pair : bpp[i])
      {
        uint j = pair.first;
        float prob = pair.second;
        
        // Only store probabilities above threshold
        if (prob > threshold() && j > i && j <= L)
        {
          bp[i - 1].push_back(std::make_pair(j - 1, prob));
        }
      }
    }
  } catch (const std::exception& e) {
    rethrow_linfold_failure(model_type_, seq.size(), e.what());
  } catch (...) {
    rethrow_linfold_failure(model_type_, seq.size(), "unknown exception");
  }
}
