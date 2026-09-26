#include "nussinov.h"
#include "relaxed_bounds.h"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <limits>
#include <numeric>
#include <random>
#include <stack>
#include <unordered_map>
#include <unordered_set>

const uint Fold::Decoder::n_support_brackets = 4 + 26;
const char* Fold::Decoder::left_brackets = "([{<ABCDEFGHIJKLMNOPQRSTUVWXYZ";
const char* Fold::Decoder::right_brackets = ")]}>abcdefghijklmnopqrstuvwxyz";

namespace {
bool close(float lhs, float rhs)
{
  return std::fabs(lhs-rhs) <= 1e-5f;
}

float exact_cached_score(const VVF& pair_score)
{
  const uint length = pair_score.size();
  if (length == 0)
    return 0.0f;
  VVF best(length, VF(length, 0.0f));
  for (uint width = 2; width < length; ++width) {
    for (uint left = 0; left + width < length; ++left) {
      const uint right = left + width;
      float value = best[left][right-1];
      for (uint pair_left = left; pair_left + 2 <= right; ++pair_left) {
        const float prefix = pair_left == left
            ? 0.0f : best[left][pair_left-1];
        const float inside = pair_left + 1 == right
            ? 0.0f : best[pair_left+1][right-1];
        value = std::max(value, prefix + inside +
                                pair_score[pair_left][right]);
      }
      best[left][right] = value;
    }
  }
  return best[0][length-1];
}

struct ReferenceLinearNussinovCell
{
  uint start;
  double score;
  uint pair_left;
  uint inside_start;
  bool paired;
};

using ReferenceLinearNussinovBeam =
    std::vector<ReferenceLinearNussinovCell>;

const ReferenceLinearNussinovCell*
reference_find_cell(const ReferenceLinearNussinovBeam& beam, uint start)
{
  const auto it = std::lower_bound(
      beam.begin(), beam.end(), start,
      [](const ReferenceLinearNussinovCell& cell, uint value) {
        return cell.start < value;
      });
  return it != beam.end() && it->start == start ? &*it : nullptr;
}

const ReferenceLinearNussinovCell*
reference_find_suffix(const ReferenceLinearNussinovBeam& beam,
                      uint minimum_start)
{
  const ReferenceLinearNussinovCell* best = nullptr;
  for (const auto& cell : beam) {
    if (cell.start < minimum_start)
      continue;
    if (!best || cell.score > best->score ||
        (cell.score == best->score && cell.start > best->start))
      best = &cell;
  }
  return best;
}

bool reference_pair_is_usable(uint left, uint right, uint length)
{
  return left < length && right < length && left < right &&
         right - left >= LinearNussinov::cached_support_minimum_pair_span;
}

struct ReferenceLinearNussinovCertificate
{
  std::vector<std::vector<double>> prefix_potentials;
  double additive_upper_bound = 0.0;
  double pruned_upper_bound = -std::numeric_limits<double>::infinity();
  size_t pruned_states = 0;

  double completion_bound(const ReferenceLinearNussinovCell& cell,
                          uint right, uint length) const
  {
    double outside = std::numeric_limits<double>::infinity();
    for (const auto& prefix : prefix_potentials) {
      outside = std::min(outside,
          prefix[cell.start] + prefix[length] - prefix[right+1]);
    }
    return cell.score + outside;
  }

  void record(const ReferenceLinearNussinovCell& cell,
              uint right, uint length)
  {
    pruned_upper_bound = std::max(
        pruned_upper_bound, completion_bound(cell, right, length));
    ++pruned_states;
  }
};

std::vector<double>
reference_tightened_vertex_cover(
    const std::vector<std::vector<std::pair<uint, double>>>& adjacency,
    const std::vector<double>& initial, bool reverse_first)
{
  std::vector<double> potential = initial;
  std::vector<double> best = initial;
  double best_total = std::accumulate(best.begin(), best.end(), 0.0);
  constexpr uint sweeps = 4;
  for (uint sweep = 0; sweep < sweeps; ++sweep) {
    const bool reverse = reverse_first != ((sweep & 1u) != 0);
    for (uint step = 0; step < potential.size(); ++step) {
      const uint i = reverse
          ? static_cast<uint>(potential.size()) - 1 - step : step;
      double tightened = 0.0;
      for (const auto& [j, value] : adjacency[i])
        tightened = std::max(tightened, value - potential[j]);
      potential[i] = tightened;
    }
    const double total =
        std::accumulate(potential.begin(), potential.end(), 0.0);
    if (total < best_total) {
      best_total = total;
      best = potential;
    }
  }
  return best;
}

ReferenceLinearNussinovCertificate
reference_make_certificate(const VVF& pair_score,
                           const VVU& pairs_by_right)
{
  const uint length = pair_score.size();
  std::vector<double> left(length, 0.0);
  std::vector<double> right(length, 0.0);
  std::vector<double> incident(length, 0.0);
  std::vector<std::vector<std::pair<uint, double>>> adjacency(length);
  for (uint j = 0; j < length; ++j) {
    for (const uint i : pairs_by_right[j]) {
      if (!reference_pair_is_usable(i, j, length))
        continue;
      const double value = std::max(0.0,
          static_cast<double>(pair_score[i][j]));
      if (!(value > 0.0))
        continue;
      left[i] = std::max(left[i], value);
      right[j] = std::max(right[j], value);
      incident[i] = std::max(incident[i], value);
      incident[j] = std::max(incident[j], value);
      adjacency[i].push_back({j, value});
      adjacency[j].push_back({i, value});
    }
  }

  std::vector<double> half_incident(length, 0.0);
  for (uint i = 0; i < length; ++i)
    half_incident[i] = 0.5 * incident[i];
  std::vector<std::vector<double>> covers;
  covers.push_back(std::move(left));
  covers.push_back(std::move(right));
  covers.push_back(half_incident);
  covers.push_back(reference_tightened_vertex_cover(
      adjacency, half_incident, false));
  covers.push_back(reference_tightened_vertex_cover(
      adjacency, half_incident, true));

  ReferenceLinearNussinovCertificate certificate;
  certificate.additive_upper_bound =
      std::numeric_limits<double>::infinity();
  for (const auto& cover : covers) {
    std::vector<double> prefix(length+1, 0.0);
    for (uint i = 0; i < length; ++i)
      prefix[i+1] = prefix[i] + cover[i];
    certificate.additive_upper_bound = std::min(
        certificate.additive_upper_bound, prefix.back());
    certificate.prefix_potentials.push_back(std::move(prefix));
  }
  if (length == 0)
    certificate.additive_upper_bound = 0.0;
  return certificate;
}

LinearNussinovResult
reference_linear_decode(const VVF& pair_score, const VVU& pairs_by_right,
                        uint beam_size, VU& ss, bool certify)
{
  const uint length = pair_score.size();
  ss.assign(length, -1u);
  if (length == 0)
    return {};

  ReferenceLinearNussinovCertificate certificate;
  if (certify)
    certificate = reference_make_certificate(pair_score, pairs_by_right);

  std::vector<ReferenceLinearNussinovBeam> beams(length);
  const uint beam_limit = std::max(1u, beam_size);
  for (uint right = 0; right < length; ++right) {
    std::unordered_map<uint, ReferenceLinearNussinovCell> candidates;
    candidates.reserve(beam_limit * 2 + pairs_by_right[right].size());
    const auto offer = [&](uint start, double score, bool paired,
                           uint left, uint inside_start) {
      const auto it = candidates.find(start);
      if (it == candidates.end()) {
        candidates.emplace(start, ReferenceLinearNussinovCell{
            start, score, left, inside_start, paired});
      } else if (score > it->second.score) {
        it->second = ReferenceLinearNussinovCell{
            start, score, left, inside_start, paired};
      }
    };

    offer(right, 0.0, false, -1u, -1u);
    if (right > 0)
      for (const auto& previous : beams[right-1])
        offer(previous.start, previous.score, false, -1u, -1u);

    for (const uint left : pairs_by_right[right]) {
      if (!reference_pair_is_usable(left, right, length))
        continue;
      const double local_score = pair_score[left][right];
      if (!(local_score > 0.0))
        continue;
      const ReferenceLinearNussinovCell* inside =
          reference_find_suffix(beams[right-1], left+1);
      if (!inside)
        continue;
      offer(left, inside->score + local_score, true, left, inside->start);
      if (left > 0) {
        for (const auto& prefix_interval : beams[left-1])
          offer(prefix_interval.start,
                prefix_interval.score + inside->score + local_score,
                true, left, inside->start);
      }
    }

    ReferenceLinearNussinovBeam beam;
    beam.reserve(candidates.size());
    for (const auto& entry : candidates)
      beam.push_back(entry.second);
    std::sort(beam.begin(), beam.end(),
              [](const ReferenceLinearNussinovCell& lhs,
                 const ReferenceLinearNussinovCell& rhs) {
                return lhs.start > rhs.start;
              });
    std::unordered_set<uint> dominated;
    double best_suffix_score = -std::numeric_limits<double>::infinity();
    ReferenceLinearNussinovBeam nondominated;
    nondominated.reserve(beam.size());
    for (const auto& cell : beam) {
      if (cell.start != 0 && cell.score <= best_suffix_score) {
        dominated.insert(cell.start);
        continue;
      }
      best_suffix_score = std::max(best_suffix_score, cell.score);
      nondominated.push_back(cell);
    }
    beam.swap(nondominated);

    const auto priority = [&](const ReferenceLinearNussinovCell& cell) {
      if (cell.start == 0)
        return cell.score;
      const ReferenceLinearNussinovCell* prefix =
          reference_find_cell(beams[cell.start-1], 0);
      assert(prefix);
      return prefix->score + cell.score;
    };
    const auto better = [&](const ReferenceLinearNussinovCell& lhs,
                            const ReferenceLinearNussinovCell& rhs) {
      const double lhs_priority = priority(lhs);
      const double rhs_priority = priority(rhs);
      return lhs_priority != rhs_priority ? lhs_priority > rhs_priority
                                          : lhs.start < rhs.start;
    };
    const auto root_it = std::find_if(
        beam.begin(), beam.end(),
        [](const ReferenceLinearNussinovCell& cell) {
          return cell.start == 0;
        });
    assert(root_it != beam.end());
    const ReferenceLinearNussinovCell root = *root_it;
    std::sort(beam.begin(), beam.end(), better);
    if (beam.size() > beam_limit)
      beam.resize(beam_limit);
    if (std::none_of(beam.begin(), beam.end(),
                     [](const ReferenceLinearNussinovCell& cell) {
                       return cell.start == 0;
                     }))
      beam.back() = root;
    if (certify) {
      std::unordered_set<uint> retained;
      retained.reserve(beam.size());
      for (const auto& cell : beam)
        retained.insert(cell.start);
      for (const auto& [start, cell] : candidates)
        if (retained.find(start) == retained.end() &&
            dominated.find(start) == dominated.end())
          certificate.record(cell, right, length);
    }
    std::sort(beam.begin(), beam.end(),
              [](const ReferenceLinearNussinovCell& lhs,
                 const ReferenceLinearNussinovCell& rhs) {
                return lhs.start < rhs.start;
              });
    beams[right] = std::move(beam);
  }

  std::vector<std::pair<uint, uint>> pending;
  pending.push_back({0, length-1});
  while (!pending.empty()) {
    const auto [start, right] = pending.back();
    pending.pop_back();
    if (start > right)
      continue;
    const ReferenceLinearNussinovCell* cell =
        reference_find_cell(beams[right], start);
    assert(cell);
    if (!cell->paired) {
      if (start < right)
        pending.push_back({start, right-1});
      continue;
    }
    const uint left = cell->pair_left;
    ss[left] = right;
    if (start < left)
      pending.push_back({start, left-1});
    if (cell->inside_start != -1u)
      pending.push_back({cell->inside_start, right-1});
  }

  const ReferenceLinearNussinovCell* root =
      reference_find_cell(beams[length-1], 0);
  assert(root);
  LinearNussinovResult result;
  result.score = static_cast<float>(root->score);
  if (!certify) {
    result.upper_bound = result.score;
    result.pruned_upper_bound = result.score;
    result.additive_upper_bound = result.score;
    return result;
  }
  const double pruned_upper_bound = std::max(
      root->score, certificate.pruned_upper_bound);
  const double upper_bound = std::min(
      pruned_upper_bound, certificate.additive_upper_bound);
  result.pruned_upper_bound =
      std::max(result.score,
               RelaxedBounds::round_up_to_float(pruned_upper_bound));
  result.additive_upper_bound =
      std::max(result.score, RelaxedBounds::round_up_to_float(
                                 certificate.additive_upper_bound));
  result.upper_bound = std::max(
      result.score, RelaxedBounds::round_up_to_float(upper_bound));
  result.pruned_states = certificate.pruned_states;
  return result;
}

bool same_float_bits(float lhs, float rhs)
{
  uint32_t lhs_bits = 0;
  uint32_t rhs_bits = 0;
  std::memcpy(&lhs_bits, &lhs, sizeof(lhs_bits));
  std::memcpy(&rhs_bits, &rhs, sizeof(rhs_bits));
  return lhs_bits == rhs_bits;
}

bool same_linear_result(const LinearNussinovResult& lhs,
                        const LinearNussinovResult& rhs)
{
  return same_float_bits(lhs.score, rhs.score) &&
         same_float_bits(lhs.upper_bound, rhs.upper_bound) &&
         same_float_bits(lhs.pruned_upper_bound, rhs.pruned_upper_bound) &&
         same_float_bits(lhs.additive_upper_bound,
                         rhs.additive_upper_bound) &&
         lhs.pruned_states == rhs.pruned_states;
}
}

int main()
{
  constexpr uint length = 10;
  constexpr float threshold = 0.2f;
  VVF probability(length, VF(length, 0.0f));
  VVF multiplier(length, VF(length, 0.0f));

  probability[0][9] = 0.90f;
  probability[1][8] = 0.80f;
  probability[3][6] = 0.65f;
  probability[0][5] = 0.75f;
  probability[4][9] = 0.70f;
  // The probability is absent, but a negative multiplier makes this pair
  // profitable.  Sparse support must therefore be the union of p and q.
  multiplier[2][7] = -0.80f;

  Nussinov exact(threshold);
  LinearNussinov linear(threshold, length);
  VU exact_structure, linear_structure;
  const float exact_score =
      exact.decode(1.0f, probability, multiplier, exact_structure);
  const float dense_linear_score =
      linear.decode(1.0f, probability, multiplier, linear_structure);
  if (!close(exact_score, dense_linear_score)) {
    std::cerr << "dense score mismatch: exact=" << exact_score
              << " linear=" << dense_linear_score << '\n';
    for (uint i = 0; i < length; ++i)
      if (exact_structure[i] != -1u || linear_structure[i] != -1u)
        std::cerr << i << ": exact=" << exact_structure[i]
                  << " linear=" << linear_structure[i] << '\n';
    return 1;
  }

  SparseFloatMatrix sparse_probability, sparse_multiplier;
  sparse_probability.assign(length, length);
  sparse_multiplier.assign(length, length);
  for (uint i = 0; i < length; ++i)
    for (uint j = i+1; j < length; ++j) {
      sparse_probability.set(i, j, probability[i][j]);
      sparse_multiplier.set(i, j, multiplier[i][j]);
    }

  const float sparse_linear_score =
      linear.decode(1.0f, sparse_probability, sparse_multiplier,
                    linear_structure);
  if (!close(exact_score, sparse_linear_score)) {
    std::cerr << "sparse score mismatch: exact=" << exact_score
              << " linear=" << sparse_linear_score << '\n';
    return 1;
  }
  if (linear_structure[2] != 7) {
    std::cerr << "q-only profitable pair was omitted\n";
    return 1;
  }

  // Final consensus decoding keeps posterior probabilities and profile
  // bonuses separate.  A bonus-only pair must enter the sparse support, and
  // exact, sparse, and full-beam linear decoders must agree on the score.
  SparseFloatMatrix final_probability, final_bonus;
  final_probability.assign(length, length);
  final_bonus.assign(length, length);
  final_probability.set(1, 8, 0.60f);
  final_bonus.set(0, 9, 0.50f);
  VU exact_bonus_structure, sparse_bonus_structure, linear_bonus_structure;
  std::string exact_bonus_brackets, sparse_bonus_brackets,
      linear_bonus_brackets;
  const float exact_bonus_score = exact.decode(
      final_probability, final_bonus,
      exact_bonus_structure, exact_bonus_brackets);
  SparseNussinov sparse_exact(threshold);
  const float sparse_bonus_score = sparse_exact.decode(
      final_probability, final_bonus,
      sparse_bonus_structure, sparse_bonus_brackets);
  const float linear_bonus_score = linear.decode(
      final_probability, final_bonus,
      linear_bonus_structure, linear_bonus_brackets);
  if (!close(exact_bonus_score, 0.70f) ||
      !close(sparse_bonus_score, exact_bonus_score) ||
      !close(linear_bonus_score, exact_bonus_score) ||
      sparse_bonus_structure != exact_bonus_structure ||
      linear_bonus_structure != exact_bonus_structure ||
      exact_bonus_structure[0] != 9 || exact_bonus_structure[1] != 8) {
    std::cerr << "final pair bonus decoding mismatch: exact="
              << exact_bonus_score << " sparse=" << sparse_bonus_score
              << " linear=" << linear_bonus_score << '\n';
    return 1;
  }

  const LinearNussinovResult full_certificate = linear.decode_certified(
      1.0f, sparse_probability, sparse_multiplier, linear_structure);
  if (!close(full_certificate.score, exact_score) ||
      !close(full_certificate.upper_bound, exact_score) ||
      full_certificate.pruned_states != 0) {
    std::cerr << "unpruned certificate mismatch: exact=" << exact_score
              << " beam=" << full_certificate.score
              << " upper=" << full_certificate.upper_bound
              << " pruned=" << full_certificate.pruned_states << '\n';
    return 1;
  }

  // DAFS caches every nonzero posterior entry and can therefore supply a
  // span-two pair.  This pair has always been decodable; its certificate must
  // cover it as well so a feasible beam score never exceeds its upper bound.
  SparseFloatMatrix short_probability, short_multiplier;
  short_probability.assign(length, length);
  short_multiplier.assign(length, length);
  short_probability.set(0, 2, 1.0f);
  VVU short_support(length);
  short_support[2].push_back(0);
  VU short_structure;
  LinearNussinov short_linear(0.0f, length);
  const LinearNussinovResult short_result = short_linear.decode_certified(
      1.0f, short_probability, short_multiplier,
      short_support, short_structure);
  if (!close(short_result.score, 1.0f) || short_structure[0] != 2 ||
      short_result.score > short_result.upper_bound + 1e-5f ||
      short_result.upper_bound >
          short_result.additive_upper_bound + 1e-5f) {
    std::cerr << "span-two certificate mismatch: lower="
              << short_result.score << " upper=" << short_result.upper_bound
              << " additive=" << short_result.additive_upper_bound << '\n';
    return 1;
  }

  // Exercise certificate bookkeeping on two orderings of one finite support.
  // This fixture has stable winners, so its score, traceback, and certificate
  // result are expected to match; it is not a general support-order tie rule.
  // Keep duplicates and malformed entries because cached support can contain
  // both, while the decoder must still use the same finite entries.
  {
    constexpr uint order_length = 18;
    SparseFloatMatrix order_probability, order_multiplier;
    order_probability.assign(order_length, order_length);
    order_multiplier.assign(order_length, order_length);
    VVU forward_support(order_length), reverse_support(order_length);
    for (uint right = 2; right < order_length; ++right) {
      for (uint left = 0; left + 2 <= right; ++left) {
        const float value = 0.08f +
            0.01f * static_cast<float>((3 * left + right) % 9);
        order_probability.set(left, right, value);
        forward_support[right].push_back(left);
        if ((left + right) % 3 == 0)
          forward_support[right].push_back(left);
      }
      forward_support[right].push_back(right);
      if (right > 0)
        forward_support[right].push_back(right - 1);
      forward_support[right].push_back(order_length + right);
      reverse_support[right] = forward_support[right];
      std::reverse(reverse_support[right].begin(),
                   reverse_support[right].end());
    }
    LinearNussinov forward_decoder(0.0f, 3);
    LinearNussinov reverse_decoder(0.0f, 3);
    VU forward_structure, reverse_structure;
    const LinearNussinovResult forward_result =
        forward_decoder.decode_certified(
            1.0f, order_probability, order_multiplier,
            forward_support, forward_structure);
    const LinearNussinovResult reverse_result =
        reverse_decoder.decode_certified(
            1.0f, order_probability, order_multiplier,
            reverse_support, reverse_structure);
    if (!same_linear_result(forward_result, reverse_result) ||
        forward_structure != reverse_structure) {
      std::cerr << "support order changed certified decoding\n";
      return 1;
    }
  }

  // Compare narrow-beam certificates with the exact cubic decoder over many
  // small deterministic random instances.  Every exact derivation either
  // reaches the root or has a first state removed from the beam; the latter
  // must be covered by the recorded completion bound.
  std::mt19937 generator(20260803u);
  std::uniform_real_distribution<float> score_distribution(0.05f, 1.0f);
  std::bernoulli_distribution edge_distribution(0.45);
  size_t strict_improvements = 0;
  size_t pruned_instances = 0;
  for (uint trial = 0; trial < 2000; ++trial) {
    constexpr uint random_length = 12;
    VVF dense_probability(random_length, VF(random_length, 0.0f));
    VVF dense_multiplier(random_length, VF(random_length, 0.0f));
    SparseFloatMatrix random_probability, random_multiplier;
    random_probability.assign(random_length, random_length);
    random_multiplier.assign(random_length, random_length);
    for (uint i = 0; i < random_length; ++i) {
      for (uint j = i+3; j < random_length; ++j) {
        const bool has_probability = edge_distribution(generator);
        const bool has_q_only_edge = !has_probability &&
            (i + 3*j + trial) % 13 == 0;
        if (!has_probability && !has_q_only_edge)
          continue;
        if (has_probability) {
          const float probability = score_distribution(generator);
          dense_probability[i][j] = probability;
          random_probability.set(i, j, probability);
        }
        if (has_q_only_edge || (i+j+trial) % 7 == 0) {
          const float multiplier = -0.25f * score_distribution(generator);
          dense_multiplier[i][j] = multiplier;
          random_multiplier.set(i, j, multiplier);
        }
      }
    }

    Nussinov random_exact(0.0f);
    VU random_exact_structure;
    const float random_exact_score = random_exact.decode(
        1.0f, dense_probability, dense_multiplier, random_exact_structure);
    VVU cached_support(random_length);
    for (uint i = 0; i < random_length; ++i)
      for (uint j = i+3; j < random_length; ++j)
        if (dense_probability[i][j] != 0.0f ||
            dense_multiplier[i][j] != 0.0f)
          cached_support[j].push_back(i);
    for (const uint beam : {1u, 2u, 3u, 5u}) {
      LinearNussinov random_linear(0.0f, beam);
      VU random_structure;
      const LinearNussinovResult result = random_linear.decode_certified(
          1.0f, random_probability, random_multiplier, random_structure);
      VU cached_structure;
      const LinearNussinovResult cached_result =
          random_linear.decode_certified(
              1.0f, random_probability, random_multiplier,
              cached_support, cached_structure);
      const float tolerance =
          2e-5f * std::max(1.0f, std::fabs(random_exact_score));
      if (result.score > random_exact_score + tolerance ||
          result.upper_bound + tolerance < random_exact_score ||
          result.upper_bound > result.additive_upper_bound + tolerance ||
          result.score > result.upper_bound + tolerance) {
        std::cerr << "invalid beam certificate: trial=" << trial
                  << " beam=" << beam
                  << " lower=" << result.score
                  << " exact=" << random_exact_score
                  << " upper=" << result.upper_bound
                  << " additive=" << result.additive_upper_bound << '\n';
        return 1;
      }
      if (!close(cached_result.score, result.score) ||
          !close(cached_result.upper_bound, result.upper_bound) ||
          cached_result.pruned_states != result.pruned_states ||
          cached_structure != random_structure) {
        std::cerr << "cached support changed certified decoding: trial="
                  << trial << " beam=" << beam << '\n';
        return 1;
      }
      strict_improvements +=
          result.upper_bound + tolerance < result.additive_upper_bound;
      pruned_instances += result.pruned_states != 0;
    }
  }
  if (strict_improvements == 0 || pruned_instances == 0) {
    std::cerr << "random certificate tests did not exercise pruning\n";
    return 1;
  }

  // Reuse one decoder across changing lengths, then repeat the first input.
  // This catches stale beam, candidate, certificate, and traceback entries in
  // the retained workspace independently of the randomized same-length loop.
  LinearNussinov reused_linear(0.0f, 3);
  const auto decode_reused = [&](uint reused_length) {
    SparseFloatMatrix reused_probability, reused_multiplier;
    reused_probability.assign(reused_length, reused_length);
    reused_multiplier.assign(reused_length, reused_length);
    VVU reused_support(reused_length);
    for (uint i = 0; i + 2 < reused_length; ++i) {
      const uint j = reused_length - 1 - i;
      if (j < i + 2)
        break;
      const float probability = 0.5f + 0.01f * i;
      reused_probability.set(i, j, probability);
      reused_support[j].push_back(i);
    }
    VU reused_structure;
    const LinearNussinovResult result = reused_linear.decode_certified(
        1.0f, reused_probability, reused_multiplier,
        reused_support, reused_structure);
    return std::make_pair(result, reused_structure);
  };
  const auto first_reused = decode_reused(13);
  decode_reused(5);
  decode_reused(19);
  const auto repeated_reused = decode_reused(13);
  if (!close(first_reused.first.score, repeated_reused.first.score) ||
      !close(first_reused.first.upper_bound,
             repeated_reused.first.upper_bound) ||
      first_reused.first.pruned_states !=
          repeated_reused.first.pruned_states ||
      first_reused.second != repeated_reused.second) {
    std::cerr << "workspace reuse changed LinearNussinov decoding\n";
    return 1;
  }

  // Exercise the cached DAFS path independently with span-two edges.  Also
  // include malformed support entries: the shared decoder/certificate
  // predicate must reject them before matrix access.  The exact DP permits
  // the same minimum span and verifies that every reported bound remains a
  // true upper bound, not merely that the internal inequalities agree.
  for (uint trial = 0; trial < 500; ++trial) {
    constexpr uint random_length = 12;
    SparseFloatMatrix random_probability, random_multiplier;
    random_probability.assign(random_length, random_length);
    random_multiplier.assign(random_length, random_length);
    VVU cached_support(random_length);
    VVF pair_score(random_length, VF(random_length, 0.0f));
    for (uint i = 0; i < random_length; ++i) {
      for (uint j = i+2; j < random_length; ++j) {
        const bool has_probability = edge_distribution(generator);
        const bool has_q_only_edge = !has_probability &&
            (i + 5*j + trial) % 11 == 0;
        if (!has_probability && !has_q_only_edge)
          continue;
        float score = 0.0f;
        if (has_probability) {
          const float probability = score_distribution(generator);
          random_probability.set(i, j, probability);
          score += probability;
        }
        if (has_q_only_edge || (i + j + trial) % 7 == 0) {
          const float multiplier = -0.25f * score_distribution(generator);
          random_multiplier.set(i, j, multiplier);
          score -= multiplier;
        }
        pair_score[i][j] = score;
        cached_support[j].push_back(i);
      }
    }
    for (uint right = 0; right < random_length; ++right) {
      cached_support[right].push_back(right);
      if (right > 0)
        cached_support[right].push_back(right-1);
      cached_support[right].push_back(random_length + right);
    }
    const float exact_score = exact_cached_score(pair_score);
    for (const uint beam : {1u, 2u, 3u, 5u}) {
      LinearNussinov random_linear(0.0f, beam);
      VU random_structure;
      const LinearNussinovResult result = random_linear.decode_certified(
          1.0f, random_probability, random_multiplier,
          cached_support, random_structure);
      const float tolerance =
          2e-5f * std::max(1.0f, std::fabs(exact_score));
      if (result.score > exact_score + tolerance ||
          result.upper_bound + tolerance < exact_score ||
          result.score > result.upper_bound + tolerance ||
          result.upper_bound > result.additive_upper_bound + tolerance) {
        std::cerr << "invalid span-two certificate: trial=" << trial
                  << " beam=" << beam
                  << " lower=" << result.score
                  << " exact=" << exact_score
                  << " upper=" << result.upper_bound
                  << " additive=" << result.additive_upper_bound << '\n';
        return 1;
      }
    }
  }

  // Differentially execute the pre-optimization beam/certificate recurrence
  // above against the reusable dense-scratch implementation.  The reference
  // intentionally keeps unordered candidate/set traversal, while the test
  // inputs include ties, reordered/duplicate support, malformed indices,
  // varying lengths, certified calls, and narrow beams that prune states.
  size_t differential_pruned_instances = 0;
  const auto compare_reference = [&](uint case_id,
                                     const SparseFloatMatrix& probability,
                                     const VVU& support,
                                     const VVF& pair_score,
                                     uint beam_size, bool certify) {
    const uint differential_length = pair_score.size();
    VVU decoder_support;
    if (!certify) {
      decoder_support.resize(differential_length);
      for (uint right = 3; right < differential_length; ++right)
        for (uint left = 0; left + 2 < right; ++left)
          decoder_support[right].push_back(left);
    }
    const VVU& reference_support = certify ? support : decoder_support;
    VU optimized_structure, reference_structure;
    LinearNussinov optimized(0.0f, beam_size);
    LinearNussinovResult optimized_result;
    if (certify) {
      SparseFloatMatrix zero_multiplier;
      zero_multiplier.assign(differential_length, differential_length);
      optimized_result = optimized.decode_certified(
          1.0f, probability, zero_multiplier, support,
          optimized_structure);
    } else {
      const VVF zero_multiplier(
          differential_length, VF(differential_length, 0.0f));
      optimized_result.score = optimized.decode(
          1.0f, probability, zero_multiplier, optimized_structure);
      optimized_result.upper_bound = optimized_result.score;
      optimized_result.pruned_upper_bound = optimized_result.score;
      optimized_result.additive_upper_bound = optimized_result.score;
    }
    const LinearNussinovResult reference_result = reference_linear_decode(
        pair_score, reference_support, beam_size, reference_structure,
        certify);
    if (!same_linear_result(optimized_result, reference_result) ||
        optimized_structure != reference_structure) {
      std::cerr << "reference differential mismatch: case=" << case_id
                << " length=" << differential_length
                << " beam=" << beam_size
                << " certified=" << certify
                << " optimized_score=" << optimized_result.score
                << " reference_score=" << reference_result.score
                << " optimized_upper=" << optimized_result.upper_bound
                << " reference_upper=" << reference_result.upper_bound
                << " optimized_pruned=" << optimized_result.pruned_states
                << " reference_pruned=" << reference_result.pruned_states
                << '\n';
      return false;
    }
    if (certify && reference_result.pruned_states != 0)
      ++differential_pruned_instances;
    return true;
  };

  uint differential_case = 0;
  for (const uint differential_length : {0u, 1u, 2u, 3u, 4u, 6u, 9u,
                                         14u}) {
    SparseFloatMatrix probability;
    probability.assign(differential_length, differential_length);
    VVU support(differential_length);
    VVF pair_score(differential_length,
                   VF(differential_length, 0.0f));
    for (uint left = 0; left < differential_length; ++left) {
      for (uint right = left + 2; right < differential_length; ++right) {
        const float value = (left + 2 * right) % 3 == 0 ? 0.5f : 0.25f;
        probability.set(left, right, value);
        pair_score[left][right] = value;
        support[right].push_back(left);
      }
    }
    for (auto& lefts : support)
      std::reverse(lefts.begin(), lefts.end());
    for (uint right = 0; right < differential_length; ++right) {
      support[right].push_back(right);
      if (right > 0)
        support[right].push_back(right-1);
      support[right].push_back(differential_length + right);
    }
    for (const uint beam_size : {0u, 1u, 2u, 4u, 16u}) {
      if (!compare_reference(differential_case++, probability, support,
                             pair_score, beam_size, true) ||
          !compare_reference(differential_case++, probability, support,
                             pair_score, beam_size, false))
        return 1;
    }
  }

  std::mt19937 differential_generator(20260919u);
  const std::array<float, 7> differential_values = {
      -0.5f, -0.125f, 0.0f, 0.125f, 0.25f, 0.5f, 1.0f};
  for (uint trial = 0; trial < 160; ++trial) {
    const uint differential_length = trial % 19;
    SparseFloatMatrix probability;
    probability.assign(differential_length, differential_length);
    VVU support(differential_length);
    VVF pair_score(differential_length,
                   VF(differential_length, 0.0f));
    for (uint left = 0; left < differential_length; ++left) {
      for (uint right = left + 2; right < differential_length; ++right) {
        const uint selector =
            (7 * left + 11 * right + 13 * trial) % 17;
        const float value = differential_values[selector %
                                                differential_values.size()];
        if (value != 0.0f) {
          probability.set(left, right, value);
          pair_score[left][right] = value;
        }
        if (selector < 13 || value != 0.0f)
          support[right].push_back(left);
      }
    }
    for (auto& lefts : support)
      std::shuffle(lefts.begin(), lefts.end(), differential_generator);
    for (uint right = 0; right < differential_length; ++right) {
      support[right].push_back(right);
      if (right > 0)
        support[right].push_back(right-1);
      support[right].push_back(differential_length + right);
    }
    for (const uint beam_size : {0u, 1u, 2u, 3u, 7u}) {
      if (!compare_reference(differential_case++, probability, support,
                             pair_score, beam_size, true))
        return 1;
    }
    if ((trial % 4) == 0) {
      if (!compare_reference(differential_case++, probability, support,
                             pair_score, 2u, false))
        return 1;
    }
  }
  if (differential_pruned_instances == 0) {
    std::cerr << "reference differential tests did not exercise pruning\n";
    return 1;
  }

  // A concentrated sparse-support row produces more than the radix cutoff
  // candidate starts at one right endpoint while keeping M = O(L).  Use
  // reversed, duplicated, and malformed support entries so the bounded
  // candidate selection and the flat score cache are checked with a sparse
  // input support; the reference fixture's dense pair-score matrix is
  // intentional.  The two sizes also exercise the same core at a fourfold
  // scale increase.
  for (const uint concentrated_length : {257u, 1025u}) {
    SparseFloatMatrix probability;
    probability.assign(concentrated_length, concentrated_length);
    VVU support(concentrated_length);
    VVF pair_score(concentrated_length,
                   VF(concentrated_length, 0.0f));
    const uint right = concentrated_length - 1;
    for (uint left = right - 2;; --left) {
      const float value = 0.125f +
          0.03125f * static_cast<float>((left * 17) % 23);
      probability.set(left, right, value);
      pair_score[left][right] = value;
      support[right].push_back(left);
      if (left == 0)
        break;
    }
    for (uint left = 0; left + 2 < concentrated_length; left += 31)
      support[right].push_back(left);
    support[right].push_back(right);
    support[right].push_back(right - 1);
    support[right].push_back(concentrated_length + 3);
    support[0].push_back(0);
    support[0].push_back(concentrated_length + 7);

    size_t support_entries = 0;
    for (const auto& row : support)
      support_entries += row.size();
    if (support_entries < concentrated_length - 2 ||
        support_entries > 2 * concentrated_length) {
      std::cerr << "concentrated support is not O(L): length="
                << concentrated_length << " entries=" << support_entries
                << '\n';
      return 1;
    }
    for (const uint beam_size : {1u, 4u, 9u}) {
      if (!compare_reference(differential_case++, probability, support,
                             pair_score, beam_size, true))
        return 1;
    }
  }

  return 0;
}
