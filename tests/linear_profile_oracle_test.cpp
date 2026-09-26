#include "linfold_wrapper.h"
#include <cmath>
#include <cstring>
#include <cstdlib>
#include <iostream>
#include <random>
#include <vector>

// The production base-class definition shares fold.cpp with non-linear models.
void Fold::Model::calculate(const std::vector<Fasta>& fa, std::vector<BP>& bp)
{
  bp.resize(fa.size());
  for (size_t i=0; i<fa.size(); ++i) calculate(fa[i].seq(), bp[i]);
}

extern "C" {
#include <ViennaRNA/fold_compound.h>
#include <ViennaRNA/model.h>
#include <ViennaRNA/part_func.h>
#include <ViennaRNA/structures/problist.h>
}

// Independent DP oracle: an unpruned one-row profile with no covariance
// must equal ViennaRNA's partition ensemble under the same energy model.
int main()
{
  // Pair eligibility is lazily cached for every length.  The two profile
  // lengths must agree on their shared prefix, including gaps and
  // incompatible pairs.
  std::vector<std::string> short_rows(3, std::string(1024, 'A'));
  constexpr char alphabet[] = "AUGC-N";
  for (size_t row=0; row<short_rows.size(); ++row)
    for (size_t column=0; column<1024; ++column)
      short_rows[row][column] = alphabet[(column + 2*row) % 6];
  auto long_rows = short_rows;
  for (auto& row : long_rows) row.push_back('A');
  const ProfileNearestNeighbor cached(short_rows);
  const ProfileNearestNeighbor uncached(long_rows);
  const ProfileNearestNeighbor one_row({short_rows.front()});
  const ProfileNearestNeighbor one_row_long({long_rows.front()});
  for (size_t i=1; i<=1024; i+=17)
    for (size_t j=i+1; j<=1024; j+=19) {
      if (cached.can_pair(i,j) != uncached.can_pair(i,j) ||
          one_row.can_pair(i,j) != one_row_long.can_pair(i,j)) {
        std::cerr << "profile pair eligibility cache differs at "
                  << i << ',' << j << '\n';
        return 1;
      }
    }

  // Two-row successor classes must give the same partition and marginals as
  // direct candidate scanning, including columns with gaps and unknown bases.
  std::vector<std::string> paired_rows(2, std::string(96, 'A'));
  const char paired_bases[] = "ACGU";
  for (size_t column = 0; column < 96; ++column) {
    paired_rows[0][column] = paired_bases[(7 * column + column / 5) % 4];
    paired_rows[1][column] = paired_rows[0][column];
    if (column % 13 == 0) paired_rows[1][column] = '-';
    if (column % 17 == 0) paired_rows[0][column] = '-';
    if (column % 19 == 0) paired_rows[1][column] = 'N';
    if (column % 23 == 0) paired_rows[1][column] =
        paired_bases[(column + 1) % 4];
    if (column >= 40 && column < 55)
      paired_rows[0][column] = paired_rows[1][column] = '-';
  }
  std::string profile_symbols(96, 'a');
  const auto nucleotide_mask = [](char base) -> unsigned {
    switch (base) {
    case 'A': return 1;
    case 'C': return 2;
    case 'G': return 4;
    case 'U': return 8;
    default: return 0;
    }
  };
  for (size_t column = 0; column < profile_symbols.size(); ++column)
    profile_symbols[column] = static_cast<char>('a' +
        (nucleotide_mask(paired_rows[0][column]) |
         nucleotide_mask(paired_rows[1][column])));
  const ProfileNearestNeighbor two_row_profile(paired_rows);
  const auto canonical_pair = [](char left, char right) {
    return (left == 'A' && right == 'U') ||
           (left == 'U' && right == 'A') ||
           (left == 'C' && right == 'G') ||
           (left == 'G' && right == 'C') ||
           (left == 'G' && right == 'U') ||
           (left == 'U' && right == 'G');
  };
  for (size_t i = 1; i <= profile_symbols.size(); ++i)
    for (size_t j = i + 1; j <= profile_symbols.size(); ++j) {
      unsigned incompatible = 0;
      unsigned double_gap = 0;
      bool canonical = false;
      for (const auto& row : paired_rows) {
        const char left = row[i - 1], right = row[j - 1];
        if (left == '-' && right == '-') ++double_gap;
        else if (canonical_pair(left, right)) canonical = true;
        else ++incompatible;
      }
      const bool expected = canonical &&
          2 * incompatible + double_gap < paired_rows.size();
      if (two_row_profile.can_pair(i, j) != expected) {
        std::cerr << "two-row eligibility differs at " << i << ',' << j << '\n';
        return 1;
      }
    }
  const auto fold_two_rows = [&](bool use_successors, bool skip_all_gap) {
    auto profile = std::make_unique<ProfileNearestNeighbor>(paired_rows);
    auto* profile_ptr = profile.get();
    LinFold<ProfileNearestNeighbor> engine(std::move(profile));
    LinFold<ProfileNearestNeighbor>::Options options;
    options.beam_size(100);
    options.alphabets("abcdefghijklmnop");
    options.occupancy_prefix_ = profile_ptr->occupancy_prefix();
    options.probability_cutoff_ = 0.0f;
    options.cache_position_pair_queries_ = false;
    if (skip_all_gap) {
      options.previous_occupied_.resize(profile_symbols.size() + 1);
      for (size_t column = 1; column <= profile_symbols.size(); ++column)
        options.previous_occupied_[column] =
            paired_rows[0][column - 1] != '-' ||
            paired_rows[1][column - 1] != '-'
            ? column : options.previous_occupied_[column - 1];
    }
    if (use_successors) {
      options.pair_signature_.resize(profile_symbols.size() + 1);
      for (size_t column = 1; column <= profile_symbols.size(); ++column)
        options.pair_signature_[column] =
            profile_ptr->two_row_signature(column);
    }
    for (char left = 'a'; left <= 'p'; ++left)
      for (char right = 'a'; right <= 'p'; ++right)
        options.set_allowed_pair(left, right);
    options.position_pair([profile_ptr](unsigned i, unsigned j) {
      return profile_ptr->can_pair(i, j);
    });
    const auto log_z = engine.compute_inside(profile_symbols, options);
    engine.compute_outside(profile_symbols, options);
    const auto bpp = engine.compute_basepairing_probabilities(
        profile_symbols, options);
    return std::make_pair(log_z, bpp);
  };
  if (fold_two_rows(false, false) != fold_two_rows(true, true)) {
    std::cerr << "two-row profile successor/gap lookup changed partition/BPP\n";
    return 1;
  }

  // Bit-mask pair counts must match a row-wise oracle for profiles on both
  // sides of the 64-row fast-path boundary, including T, N and gaps.
  std::mt19937 mask_random(57117);
  constexpr char profile_alphabet[] = "ACGUT-N";
  const auto oriented_type = [&](char left, char right) -> unsigned {
    if (left == 'T') left = 'U';
    if (right == 'T') right = 'U';
    if (left == '-' && right == '-') return 7;
    if (left == 'C' && right == 'G') return 1;
    if (left == 'G' && right == 'C') return 2;
    if (left == 'G' && right == 'U') return 3;
    if (left == 'U' && right == 'G') return 4;
    if (left == 'A' && right == 'U') return 5;
    if (left == 'U' && right == 'A') return 6;
    return 0;
  };
  for (size_t row_count : {3u, 10u, 64u, 65u}) {
    std::vector<std::string> rows(row_count, std::string(96, 'A'));
    for (auto& row : rows)
      for (char& base : row)
        base = profile_alphabet[mask_random() % 7];
    ProfileNearestNeighbor profile(rows);
    for (size_t i = 1; i <= 96; ++i)
      for (size_t j = i + 1; j <= 96; ++j) {
        std::array<unsigned, 8> expected{};
        for (const auto& row : rows)
          ++expected[oriented_type(row[i - 1], row[j - 1])];
        const unsigned canonical = expected[1] + expected[2] +
            expected[3] + expected[4] + expected[5] + expected[6];
        const bool can_pair = canonical != 0 &&
            2 * expected[0] + expected[7] < row_count;
        if (profile.pair_type_counts(i, j) != expected ||
            profile.can_pair(i, j) != can_pair) {
          std::cerr << "profile bit-mask counts differ at " << row_count
                    << " rows, " << i << ',' << j << '\n';
          return 1;
        }
      }
    std::vector<double> scratch;
    for (size_t m = 2; m <= 20; ++m) {
      const double original = profile.score_helix(45 - (m - 1),
                                                    52 + (m - 1), m);
      const double reused = profile.score_helix_extend(45, 52, m, scratch);
      if (std::memcmp(&original, &reused, sizeof(double)) != 0) {
        std::cerr << "profile helix energy differs at " << row_count
                  << " rows, length " << m << '\n';
        return 1;
      }
    }
  }

  std::vector<std::string> sequences{
      "GAAAC", "GGGAAACCC", "GUGAAACAC", "GGAAACGAAACC",
      "GGGAAACCCGGGAAACCC", "GGGGAAACCCGGGAAACCCC",
      "GGGGGGGAAAAACCCCCCC", "AAAAA", "A", "AC"};
  std::mt19937 random(178979);
  const char bases[] = "ACGU";
  for (unsigned sample=0; sample<120; ++sample) {
    std::string sequence(sample<80 ? 5+sample%20 : 25+sample%16, 'A');
    for (auto& base : sequence) base = bases[random()%4];
    sequences.push_back(std::move(sequence));
  }
  for (unsigned length : {64u, 80u}) {
    std::string sequence(length, 'A');
    for (auto& base : sequence) base = bases[random()%4];
    sequences.push_back(std::move(sequence));
  }
  // Reuse one folding engine across varying lengths and sequences; every
  // result must still agree with the independent ViennaRNA oracle.
  LinFoldWrapper fold(0.0f, LinFoldWrapper::ModelType::LPV, 0);
  LinFoldWrapper fold_rnaalifold(
      0.0f, LinFoldWrapper::ModelType::LPV, 0,
      LinFoldWrapper::ProfileEnergyMode::RNAalifold);
  for (const auto& sequence : sequences) {
    vrna_md_t model;
    vrna_md_set_default(&model);
    model.temperature = 37.0;
    model.dangles = 2;
    model.special_hp = 1;
    model.pf_smooth = 0;
    model.compute_bpp = 1;
    auto* compound = vrna_fold_compound(sequence.c_str(), &model, VRNA_OPTION_PF);
    // Keep the generic triloop rule but disable sequence-specific bonuses.
    compound->exp_params->Triloops[0] = '\0';
    compound->exp_params->Tetraloops[0] = '\0';
    compound->exp_params->Hexaloops[0] = '\0';
    const double free_energy = vrna_pf(compound, nullptr);
    auto* probabilities = vrna_plist_from_probs(compound, 0.0);
    std::vector<std::vector<double>> expected(sequence.size(),
        std::vector<double>(sequence.size(), 0.0));
    for (auto* pair=probabilities; pair->i; ++pair)
      expected[pair->i-1][pair->j-1] = pair->p;
    free(probabilities);
    vrna_fold_compound_free(compound);

    const std::vector<Fasta> rows{{"one", sequence}};
    const ALN alignment{{0u, std::vector<bool>(sequence.size(), true)}};
    BP actual;
    const double log_z = fold.calculate_profile(alignment, rows, actual);
    const double expected_log_z = -free_energy / 0.6163207755;
    if (std::abs(log_z-expected_log_z) > 3e-5*(1+std::abs(expected_log_z))) {
      std::cerr << "logZ differs from ViennaRNA for " << sequence
                << ": " << log_z << " versus " << expected_log_z << '\n';
      return 1;
    }
    std::vector<std::vector<double>> observed(sequence.size(),
        std::vector<double>(sequence.size(), 0.0));
    for (size_t i=0; i<actual.size(); ++i)
      for (const auto& [j,p] : actual[i]) observed[i][j] = p;
    for (size_t i=0; i<sequence.size(); ++i)
      for (size_t j=i+1; j<sequence.size(); ++j)
        if (std::abs(observed[i][j]-expected[i][j]) > 3e-5) {
          std::cerr << "BPP differs from ViennaRNA for " << sequence
                    << " pair " << i << ',' << j << ": "
                    << observed[i][j] << " versus " << expected[i][j] << '\n';
          return 1;
        }

    // The RNAalifold-compatible profile mode must retain ViennaRNA's
    // sequence-dependent special tri-/tetra-/hexaloop terms.  Recompute the
    // same one-row ensemble without clearing those tables and compare both
    // logZ and every base-pair marginal against the independent oracle.
    vrna_md_t rnaalifold_model;
    vrna_md_set_default(&rnaalifold_model);
    rnaalifold_model.temperature = 37.0;
    rnaalifold_model.dangles = 2;
    rnaalifold_model.special_hp = 1;
    rnaalifold_model.pf_smooth = 0;
    rnaalifold_model.compute_bpp = 1;
    auto* rnaalifold_compound = vrna_fold_compound(
        sequence.c_str(), &rnaalifold_model, VRNA_OPTION_PF);
    const double rnaalifold_free_energy =
        vrna_pf(rnaalifold_compound, nullptr);
    auto* rnaalifold_probabilities = vrna_plist_from_probs(
        rnaalifold_compound, 0.0);
    std::vector<std::vector<double>> rnaalifold_expected(sequence.size(),
        std::vector<double>(sequence.size(), 0.0));
    for (auto* pair = rnaalifold_probabilities; pair->i; ++pair)
      rnaalifold_expected[pair->i - 1][pair->j - 1] = pair->p;
    free(rnaalifold_probabilities);
    vrna_fold_compound_free(rnaalifold_compound);

    BP rnaalifold_actual;
    const double rnaalifold_log_z = fold_rnaalifold.calculate_profile(
        alignment, rows, rnaalifold_actual);
    const double rnaalifold_expected_log_z =
        -rnaalifold_free_energy / 0.6163207755;
    if (std::abs(rnaalifold_log_z - rnaalifold_expected_log_z) >
        3e-5 * (1 + std::abs(rnaalifold_expected_log_z))) {
      std::cerr << "RNAalifold-compatible logZ differs from ViennaRNA for "
                << sequence << ": " << rnaalifold_log_z << " versus "
                << rnaalifold_expected_log_z << '\n';
      return 1;
    }
    std::vector<std::vector<double>> rnaalifold_observed(sequence.size(),
        std::vector<double>(sequence.size(), 0.0));
    for (size_t i = 0; i < rnaalifold_actual.size(); ++i)
      for (const auto& [j, probability] : rnaalifold_actual[i])
        rnaalifold_observed[i][j] = probability;
    for (size_t i = 0; i < sequence.size(); ++i)
      for (size_t j = i + 1; j < sequence.size(); ++j)
        if (std::abs(rnaalifold_observed[i][j] -
                     rnaalifold_expected[i][j]) > 3e-5) {
          std::cerr << "RNAalifold-compatible BPP differs from ViennaRNA for "
                    << sequence << " pair " << i << ',' << j << ": "
                    << rnaalifold_observed[i][j] << " versus "
                    << rnaalifold_expected[i][j] << '\n';
          return 1;
        }
  }
  std::cout << "ViennaRNA partition/BPP oracle: " << sequences.size() << " cases passed\n";
}
