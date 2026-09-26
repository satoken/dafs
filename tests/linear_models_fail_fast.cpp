#include "linearalign.h"
#include "linearalign/BeamAlign.h"
#include "linfold_wrapper.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <exception>
#include <functional>
#include <iostream>
#include <limits>
#include <string>
#include <array>
#include <unordered_map>
#include <utility>
#include <vector>

extern "C" {
#include <ViennaRNA/partfunc/global.h>
#include <ViennaRNA/structures/problist.h>
}

// The production definitions live in fold.cpp/align.cpp together with the
// non-linear model implementations.  Keep this focused test independent of
// those optional libraries while satisfying the two base-class vtables.
void Fold::Model::calculate(const std::vector<Fasta>& fa, std::vector<BP>& bp)
{
  bp.resize(fa.size());
  for (size_t i = 0; i < fa.size(); ++i)
    calculate(fa[i].seq(), bp[i]);
}

void Align::Model::calculate(
    const std::vector<Fasta>& fa, std::vector<std::vector<MP>>& mp)
{
  mp.resize(fa.size(), std::vector<MP>(fa.size()));
  for (size_t i = 0; i < fa.size(); ++i)
    for (size_t j = i + 1; j < fa.size(); ++j)
      calculate(fa[i].seq(), fa[j].seq(), mp[i][j]);
}

namespace
{
template <typename Function>
bool throws_with(Function&& function, const std::string& expected)
{
  try {
    function();
  } catch (const std::exception& e) {
    if (std::string(e.what()).find(expected) != std::string::npos)
      return true;
    std::cerr << "unexpected exception: " << e.what() << '\n';
    return false;
  }
  std::cerr << "expected an exception containing: " << expected << '\n';
  return false;
}

struct ReferenceMaxNode
{
  double score = std::numeric_limits<double>::lowest();
  std::uint64_t previous = 0;
  char operation = 0;
};

struct ReferenceMaxResult
{
  double score = std::numeric_limits<double>::lowest();
  std::vector<unsigned> mapping;
};

// Preserve the former per-call unordered_map implementation as an oracle.
// It intentionally remains independent of BeamAlign's private storage so the
// optimized implementation is checked against the old update, pruning, and
// traceback order.
ReferenceMaxResult reference_max_alignment(
    unsigned length1, unsigned length2, int beam,
    const MatchScoreFunction& match_score)
{
  const unsigned stride = length2 + 1;
  const auto key_of = [stride](unsigned i, unsigned k) {
    return static_cast<std::uint64_t>(i) * stride + k;
  };
  std::vector<std::unordered_map<std::uint64_t, ReferenceMaxNode>> layers(
      length1 + length2 + 1);
  layers[0][0].score = 0.0;

  const auto update_max = [](ReferenceMaxNode& node, double score,
                             std::uint64_t previous, char operation) {
    if (node.score < score) {
      node.score = score;
      node.previous = previous;
      node.operation = operation;
    }
  };

  for (unsigned step = 0; step < length1 + length2; ++step) {
    auto& layer = layers[step];
    if (beam > 0 && layer.size() > static_cast<std::size_t>(beam)) {
      std::vector<std::pair<double, std::uint64_t>> ranked;
      ranked.reserve(layer.size());
      for (const auto& [key, node] : layer)
        ranked.emplace_back(node.score, key);
      std::sort(ranked.begin(), ranked.end(),
                [](const auto& lhs, const auto& rhs) {
                  return lhs.first != rhs.first
                       ? lhs.first > rhs.first
                       : lhs.second < rhs.second;
                });
      for (std::size_t rank = beam; rank < ranked.size(); ++rank)
        layer.erase(ranked[rank].second);
    }

    std::vector<std::uint64_t> keys;
    keys.reserve(layer.size());
    for (const auto& [key, node] : layer)
      keys.push_back(key);
    std::sort(keys.begin(), keys.end());
    for (const std::uint64_t key : keys) {
      const ReferenceMaxNode& node = layer.at(key);
      const unsigned i = key / stride;
      const unsigned k = key % stride;
      if (i < length1 && k < length2) {
        const std::uint64_t next = key_of(i + 1, k + 1);
        update_max(layers[step + 2][next],
                   node.score + match_score(i, k), key, 'M');
      }
      if (i < length1) {
        const std::uint64_t next = key_of(i + 1, k);
        update_max(layers[step + 1][next], node.score, key, 'X');
      }
      if (k < length2) {
        const std::uint64_t next = key_of(i, k + 1);
        update_max(layers[step + 1][next], node.score, key, 'Y');
      }
    }
  }

  ReferenceMaxResult result;
  result.mapping.assign(length1, -1u);
  std::uint64_t key = key_of(length1, length2);
  unsigned step = length1 + length2;
  const auto terminal = layers[step].find(key);
  if (terminal == layers[step].end())
    return result;

  result.score = terminal->second.score;
  while (step > 0) {
    const ReferenceMaxNode& node = layers[step].at(key);
    const unsigned i = key / stride;
    const unsigned k = key % stride;
    if (node.operation == 'M') {
      result.mapping[i - 1] = k - 1;
      step -= 2;
    } else {
      --step;
    }
    key = node.previous;
  }
  return result;
}
}

int main()
{
  // max_alignment reuses its layer and ranking buffers.  Alternate dimensions
  // and then repeat the first call to catch stale nodes or traceback state.
  BeamAlign beam_aligner(100);
  std::vector<unsigned> mapping;
  const auto diagonal_score = [](int i, int j) {
    return i == j ? 2.0 : -1.0;
  };
  const double first_alignment =
      beam_aligner.max_alignment(6, 6, mapping, diagonal_score);
  const std::vector<unsigned> expected_mapping = {0, 1, 2, 3, 4, 5};
  if (first_alignment != 12.0 || mapping != expected_mapping)
    return 1;
  beam_aligner.max_alignment(
      3, 8, mapping,
      [](int i, int j) { return i + 2 == j ? 3.0 : -2.0; });
  const double repeated_alignment =
      beam_aligner.max_alignment(6, 6, mapping, diagonal_score);
  if (repeated_alignment != first_alignment || mapping != expected_mapping)
    return 1;

  const auto make_scores = [](unsigned length1, unsigned length2,
                              const auto& score_at) {
    std::vector<double> scores(static_cast<size_t>(length1) * length2);
    for (unsigned i = 0; i < length1; ++i)
      for (unsigned k = 0; k < length2; ++k)
        scores[static_cast<size_t>(i) * length2 + k] = score_at(i, k);
    return scores;
  };
  const auto differential_case = [&](unsigned length1, unsigned length2,
                                      int beam, const std::vector<double>& scores,
                                      const char* label) {
    const MatchScoreFunction score_at = [&](int i, int k) {
      return scores[static_cast<size_t>(i) * length2 + k];
    };
    const ReferenceMaxResult expected =
        reference_max_alignment(length1, length2, beam, score_at);
    beam_aligner.beam = beam;
    std::vector<unsigned> actual_mapping;
    const double actual_score =
        beam_aligner.max_alignment(length1, length2, actual_mapping, score_at);
    if (beam_aligner.max_score()!=actual_score) return false;
    if ((beam<=0 || static_cast<unsigned>(beam)>=std::min(length1,length2)+1) &&
        beam_aligner.max_pruned_states()!=0) return false;
    if (length1>2 && length2>2 && beam==1 && beam_aligner.max_pruned_states()==0) return false;
    if (actual_score != expected.score || actual_mapping != expected.mapping) {
      std::cerr << "max_alignment differential mismatch in " << label
                << " (lengths=" << length1 << ',' << length2
                << ", beam=" << beam << ")\n";
      return false;
    }
    return true;
  };

  // Exercise exact ties, negative and mixed scores, all beam modes, pruning,
  // and changing dimensions on the same reused BeamAlign object.
  const std::vector<double> all_zero =
      make_scores(7, 7, [](unsigned, unsigned) { return 0.0; });
  for (const int beam : {0, 1, 2, 5, 100})
    if (!differential_case(7, 7, beam, all_zero, "all-zero ties"))
      return 1;

  // Invalid match offers must remain below finite gap paths and must not
  // evade a positive-beam cutoff.  The sentinel is a finite lowest() score,
  // while the callback itself supplies NaN.
  beam_aligner.beam = 1;
  std::vector<unsigned> invalid_mapping;
  const double invalid_score = beam_aligner.max_alignment(
      3, 3, invalid_mapping,
      [](int, int) { return std::numeric_limits<double>::quiet_NaN(); });
  if (invalid_score != 0.0 ||
      !std::all_of(invalid_mapping.begin(), invalid_mapping.end(),
                   [](unsigned value) { return value == -1u; })) {
    std::cerr << "NaN match offers corrupted the finite gap path\n";
    return 1;
  }

  const std::vector<double> negative = make_scores(
      8, 5, [](unsigned i, unsigned k) {
        return -static_cast<double>((i * 7 + k * 11) % 13) -
               (i == k ? 0.25 : 0.0);
      });
  for (const int beam : {0, 1, 3, 6})
    if (!differential_case(8, 5, beam, negative, "negative scores"))
      return 1;

  const std::vector<double> mixed = make_scores(
      5, 9, [](unsigned i, unsigned k) {
        return static_cast<double>((i * 19 + k * 7) % 17) - 8.0;
      });
  for (const int beam : {0, 1, 2, 4, 100})
    if (!differential_case(5, 9, beam, mixed, "mixed scores"))
      return 1;

  const std::vector<double> repeated_values = make_scores(
      3, 4, [](unsigned i, unsigned k) {
        return (i + k) % 2 == 0 ? -2.0 : 1.0;
      });
  for (const auto dimensions : {std::pair<unsigned, unsigned>{0, 0},
                                {0, 6}, {6, 0}}) {
    const std::vector<double> zero_dimension_scores;
    for (const int beam : {0, 1, 4})
      if (!differential_case(dimensions.first, dimensions.second, beam,
                             zero_dimension_scores, "zero dimension"))
        return 1;
  }
  for (const int beam : {0, 1, 2, 4})
    if (!differential_case(3, 4, beam, repeated_values, "repeated values"))
      return 1;

  // The optimized path gathers match offers before merging gaps, but it must
  // not change the observable match-score callback order.  Compare it with
  // the independent reference while recording every callback.
  std::vector<std::pair<int, int>> reference_calls;
  const MatchScoreFunction reference_callback = [&](int i, int k) {
    reference_calls.emplace_back(i, k);
    return static_cast<double>((i * 5 + k * 3) % 11) - 5.0;
  };
  const ReferenceMaxResult callback_expected =
      reference_max_alignment(6, 5, 0, reference_callback);
  std::vector<std::pair<int, int>> actual_calls;
  const MatchScoreFunction actual_callback = [&](int i, int k) {
    actual_calls.emplace_back(i, k);
    return static_cast<double>((i * 5 + k * 3) % 11) - 5.0;
  };
  beam_aligner.beam = 0;
  std::vector<unsigned> callback_mapping;
  const double callback_score =
      beam_aligner.max_alignment(6, 5, callback_mapping, actual_callback);
  if (callback_score != callback_expected.score ||
      callback_mapping != callback_expected.mapping ||
      actual_calls != reference_calls) {
    std::cerr << "max_alignment callback order mismatch\n";
    return 1;
  }

  // Deterministic small fuzz cases vary sparse beam shapes and deliberately
  // reuse many tied score values against the independent map oracle.
  std::uint32_t random_state = 0x9e3779b9u;
  for (unsigned case_index = 0; case_index < 32; ++case_index) {
    const unsigned length1 = 1 + (case_index * 7) % 14;
    const unsigned length2 = 1 + (case_index * 11 + 3) % 14;
    const std::vector<double> random_scores = make_scores(
        length1, length2, [&](unsigned, unsigned) {
          random_state = random_state * 1664525u + 1013904223u;
          return static_cast<double>(random_state % 17) - 8.0;
        });
    for (const int beam : {0, 1, 2, 5, 100})
      if (!differential_case(length1, length2, beam, random_scores,
                             "deterministic sparse fuzz"))
        return 1;
  }

  const std::vector<double> wide_repeated = make_scores(
      101, 103, [](unsigned i, unsigned k) {
        return static_cast<double>((i * 3 + k * 5) % 9) - 4.0;
      });
  for (const int beam : {0, 1, 32, 100, 1000})
    if (!differential_case(101, 103, beam, wide_repeated,
                           "wide pruning and cutoff ties"))
      return 1;

  LinFoldWrapper fold(0.01f, LinFoldWrapper::ModelType::LPC, 10);
  BP bp;
  if (!throws_with(
          [&] { fold.calculate("ACGU", "???", bp); },
          "LinearPartition LPC failed"))
    return 1;

  LinearAlign align(0.01f, 10);
  align.setHMMParameters(nullptr, nullptr);
  MP mp;
  if (!throws_with(
          [&] { align.calculate("ACGU", "ACGU", mp); },
          "LinearAlign failed"))
    return 1;

  const std::string sequence = "GGGAAACCC";
  for (const auto model : {LinFoldWrapper::ModelType::LPV,
                           LinFoldWrapper::ModelType::LPC}) {
    LinFoldWrapper working_fold(0.01f, model, 100);
    working_fold.calculate(sequence, bp);
    size_t pair_count = 0;
    std::vector<float> marginal(sequence.size(), 0.0f);
    if (bp.size() != sequence.size())
      return 1;
    for (size_t i = 0; i < bp.size(); ++i) {
      for (const auto& [j, probability] : bp[i]) {
        if (j <= i || j >= sequence.size() || !std::isfinite(probability) ||
            probability <= 0.0f || probability > 1.0f)
          return 1;
        marginal[i] += probability;
        marginal[j] += probability;
        ++pair_count;
      }
    }
    if (pair_count == 0)
      return 1;
    for (const float probability : marginal)
      if (probability > 1.0001f)
        return 1;

    // An independent two-structure ensemble: GAAAC is either unpaired or
    // a single GC hairpin with energy 5.4 kcal/mol (Turner 2004, 37 C).
    LinFoldWrapper profile_fold(0.0f, model, 100);
    const std::vector<Fasta> hairpin_sequences{{"hairpin", "GAAAC"}};
    const ALN hairpin_row{{0u, std::vector<bool>(5, true)}};
    BP profile_bp;
    const double hairpin_weight = std::exp(-5.4 / 0.6163207755);
    const double log_z = profile_fold.calculate_profile(hairpin_row, hairpin_sequences, profile_bp);
    if (std::abs(log_z - std::log1p(hairpin_weight)) > 1e-7 ||
        profile_bp[0].size() != 1 || profile_bp[0][0].first != 4 ||
        std::abs(profile_bp[0][0].second - hairpin_weight / (1+hairpin_weight)) > 1e-7) {
      std::cerr << "profile partition differs from the two-structure oracle: logZ=" << log_z
                << " expected=" << std::log1p(hairpin_weight) << " pairs=" << profile_bp[0].size() << '\n';
      return 1;
    }
    const double forced_log_z = profile_fold.calculate_profile(
        hairpin_row, hairpin_sequences, profile_bp, "(...)");
    if (std::abs(forced_log_z + 5.4 / 0.6163207755) > 1e-6 ||
        profile_bp[0].size() != 1 || std::abs(profile_bp[0][0].second - 1) > 1e-6)
      return 1;
    if (profile_fold.calculate_profile(hairpin_row, hairpin_sequences,
                                       profile_bp, ".....") != 0.0)
      return 1;
    for (const auto& row : profile_bp) if (!row.empty()) return 1;

    // Two-row analytic ensembles also exercise ViennaRNA RNAalifold's
    // covariance normalization and type-7 thermodynamics for a double-gap
    // projected pair.  Conserved pair types have zero covariance bonus.
    const auto check_two_rows = [&](const std::string& second,
                                    const std::vector<bool>& present,
                                    double effective_energy) {
      const std::vector<Fasta> rows{{"one", "GAAAC"}, {"two", second}};
      const ALN alignment{{0u, std::vector<bool>(5, true)}, {1u, present}};
      BP pairs;
      const double z = profile_fold.calculate_profile(alignment, rows, pairs);
      const double weight = std::exp(-effective_energy / 0.6163207755);
      return std::abs(z - std::log1p(weight)) < 1e-7 &&
             pairs[0].size() == 1 && pairs[0][0].first == 4 &&
             std::abs(pairs[0][0].second - weight/(1+weight)) < 1e-7;
    };
    if (!check_two_rows("GAAAC", std::vector<bool>(5, true), 5.4) ||
        !check_two_rows("AAA", {false, true, true, true, false}, 5.9 + 0.125)) {
      std::cerr << "two-row profile differs from analytic ensemble\n";
      return 1;
    }

    const std::vector<Fasta> incompatible_rows{{"one", "GAAAC"},
                                                {"two", "AAAAA"}};
    const ALN incompatible_alignment{{0u, std::vector<bool>(5, true)},
                                      {1u, std::vector<bool>(5, true)}};
    if (profile_fold.calculate_profile(incompatible_alignment,
                                       incompatible_rows, profile_bp) != 0.0 ||
        std::any_of(profile_bp.begin(), profile_bp.end(),
                    [](const auto& row) { return !row.empty(); })) {
      std::cerr << "RNAalifold counterexample pair was not suppressed\n";
      return 1;
    }

    // Short comparative profiles must follow ViennaRNA's default
    // RNAalifold covariance model rather than DAFS's RIBOSUM objective.
    const auto compare_vienna_profile = [&](const std::vector<std::string>& rows) {
      const size_t length = rows.front().size();
      std::vector<Fasta> profile_rows;
      ALN profile_alignment;
      std::vector<const char*> vienna_rows;
      for (size_t row = 0; row < rows.size(); ++row) {
        profile_rows.emplace_back("row" + std::to_string(row), rows[row]);
        profile_alignment.emplace_back(static_cast<unsigned>(row),
                                       std::vector<bool>(length, true));
        vienna_rows.push_back(rows[row].c_str());
      }
      vienna_rows.push_back(nullptr);
      LinFoldWrapper comparative_fold(0.0f, LinFoldWrapper::ModelType::LPV, 100);
      BP comparative_bp;
      const double linear_log_z = comparative_fold.calculate_profile(
          profile_alignment, profile_rows, comparative_bp);
      vrna_ep_t* plist = nullptr;
      const float vienna_g = vrna_pf_alifold(vienna_rows.data(), nullptr, &plist);
      constexpr double rt = 0.6163207755;
      const double vienna_log_z = -vienna_g / rt;
      if (std::abs(linear_log_z - vienna_log_z) > 0.01) {
        std::cerr << "LinearAlifold logZ differs from RNAalifold: "
                  << linear_log_z << " versus " << vienna_log_z << '\n';
        free(plist);
        return false;
      }
      for (size_t i = 0; i < length; ++i) {
        for (size_t j = i + 1; j < length; ++j) {
          float linear_probability = 0.0f;
          for (const auto& [partner, probability] : comparative_bp[i])
            if (partner == j) linear_probability = probability;
          float vienna_probability = 0.0f;
          for (auto* pair = plist; pair && pair->i; ++pair)
            if (pair->i == i + 1 && pair->j == j + 1)
              vienna_probability = pair->p;
          if (std::abs(linear_probability - vienna_probability) > 0.01f) {
            std::cerr << "LinearAlifold BPP differs from RNAalifold at "
                      << i << ',' << j << '\n';
            free(plist);
            return false;
          }
        }
      }
      free(plist);
      return true;
    };
    if (!compare_vienna_profile({"GGGAAACCC", "GGGAAACCC"}) ||
        !compare_vienna_profile({"GGGAAACCC", "GGGAAUCCC"}) ||
        !compare_vienna_profile({"GGGAAACCC", "CCCAAAGGG"}) ||
        !compare_vienna_profile({"GGGAAAAACCC", "GGGAAAAACCC"}))
      return 1;

    const std::vector<Fasta> minority_rows{{"one", "GAAAC"},
                                          {"two", "AAAAA"}, {"three", "AAAAA"}};
    const ALN minority_alignment{{0u, std::vector<bool>(5, true)},
                                  {1u, std::vector<bool>(5, true)},
                                  {2u, std::vector<bool>(5, true)}};
    if (!throws_with([&] {
          profile_fold.calculate_profile(minority_alignment, minority_rows,
                                          profile_bp, "(...)");
        }, "infeasible fixed pair")) return 1;
    if (profile_fold.calculate_profile(minority_alignment, minority_rows,
                                       profile_bp, "(...)", true) != 0.0)
      return 1;
    for (const auto& row : profile_bp) if (!row.empty()) return 1;

    // An all-gap column changes alignment coordinates, not this ensemble.
    ALN gap_column{{0u, {true, true, false, true, true, true}}};
    const double gap_log_z = profile_fold.calculate_profile(
        gap_column, hairpin_sequences, profile_bp);
    if (std::abs(gap_log_z-log_z) > 1e-7 ||
        profile_bp[0].size() != 1 || profile_bp[0][0].first != 5 ||
        std::abs(profile_bp[0][0].second-hairpin_weight/(1+hairpin_weight)) > 1e-7)
      return 1;

    if (!throws_with([&] {
          profile_fold.calculate_profile(hairpin_row, hairpin_sequences,
                                          profile_bp, "(??)?");
        }, "infeasible fixed pair"))
      return 1;

    // Exercise the extended-helix branch beyond the default max_helix=30.
    const std::string stem_sequence = std::string(34, 'G') + "AAA" + std::string(34, 'C');
    const std::string stem_constraint = std::string(34, '(') + "..." + std::string(34, ')');
    const std::vector<Fasta> stem_sequences{{"long-stem", stem_sequence}};
    const ALN stem_alignment{{0u, std::vector<bool>(stem_sequence.size(), true)}};
    const double stem_log_z = profile_fold.calculate_profile(
        stem_alignment, stem_sequences, profile_bp, stem_constraint);
    if (!std::isfinite(stem_log_z)) return 1;
    for (size_t i=0; i<34; ++i)
      if (profile_bp[i].size()!=1 || profile_bp[i][0].first!=stem_sequence.size()-1-i ||
          std::abs(profile_bp[i][0].second-1)>1e-5)
        return 1;

    // 45 alignment columns contain only 22.5 nucleotides on average.  The
    // profile loop limit must admit this forced internal loop; the short row
    // sees a stack and the long row a bulge beyond the table's length 30.
    const std::vector<Fasta> long_gap_sequences{
        {"short", "GGAAACC"}, {"long", "G" + std::string(45, 'A') + "GAAACC"}};
    ALN long_gap_alignment{{0u, std::vector<bool>(52, true)},
                            {1u, std::vector<bool>(52, true)}};
    for (size_t i=1; i<=45; ++i) long_gap_alignment[0].second[i] = false;
    const double long_gap_log_z = profile_fold.calculate_profile(
        long_gap_alignment, long_gap_sequences, profile_bp,
        "(" + std::string(45, '.') + "(...))");
    if (!std::isfinite(long_gap_log_z)) return 1;
    for (const auto i : {0u, 46u}) {
      const auto j = i==0 ? 51u : 50u;
      if (profile_bp[i].size()!=1 || profile_bp[i][0].first!=j ||
          std::abs(profile_bp[i][0].second-1)>1e-5)
        return 1;
    }

    const std::vector<Fasta> covarying{{"one", sequence},
                                      {"two", "GGGAAAUCC"}};
    const ALN two_rows{{0u, std::vector<bool>(sequence.size(), true)},
                       {1u, std::vector<bool>(sequence.size(), true)}};
    working_fold.calculate_profile(two_rows, covarying, profile_bp);
    if (profile_bp.size() != sequence.size() || profile_bp[0].empty())
      return 1;
    working_fold.calculate_profile(two_rows, covarying, profile_bp,
                                   std::string(sequence.size(), '.'));
    for (const auto& row : profile_bp)
      if (!row.empty()) return 1;
    working_fold.calculate_profile(two_rows, covarying, profile_bp,
                                   std::string(1, '(') + std::string(7, '?') + ")");
    bool forced_pair_found = false;
    for (const auto& [j, probability] : profile_bp[0])
      if (j == sequence.size() - 1 && probability > 0.5f)
        forced_pair_found = true;
    if (!forced_pair_found) {
      std::cerr << "profile constraint did not retain its forced pair\n";
      return 1;
    }

    const std::vector<Fasta> gapped_sequences{{"one", sequence},
                                             {"short", "GGAAUCCC"}};
    ALN gapped{{0u, std::vector<bool>(sequence.size(), true)},
               {1u, std::vector<bool>(sequence.size(), true)}};
    gapped[1].second[0] = false;
    working_fold.calculate_profile(gapped, gapped_sequences, profile_bp);
    if (profile_bp.size() != sequence.size()) return 1;
    for (size_t i = 0; i < profile_bp.size(); ++i)
      for (const auto& [j, probability] : profile_bp[i])
        if (j <= i || j >= sequence.size() ||
            !std::isfinite(probability) || probability <= 0 || probability > 1)
          return 1;
  }

  LinearAlign working_align(0.01f, 100);
  working_align.calculate(sequence, "GGGAAUCCC", mp);
  size_t match_count = 0;
  for (size_t i = 0; i < mp.size(); ++i) {
    for (const auto& [j, probability] : mp[i]) {
      if (j >= sequence.size() || !std::isfinite(probability) ||
          probability <= 0.0f || probability > 1.0f)
        return 1;
      ++match_count;
    }
  }
  if (match_count == 0)
    return 1;

  if (LinearAlign::parseScoreModel("TurboFold") !=
          LinearAlign::ScoreModel::LinearTurboFold ||
      LinearAlign::parseScoreModel("contralign") !=
          LinearAlign::ScoreModel::CONTRAlign ||
      LinearAlign::parseScoreModel("ProbCons-RNA") !=
          LinearAlign::ScoreModel::ProbConsRNA)
    return 1;
  if (!throws_with(
          [] { LinearAlign::parseScoreModel("unknown"); },
          "unknown LinearAlign score model"))
    return 1;

  std::array<double, 3> posterior_sums{};
  const std::array<LinearAlign::ScoreModel, 3> score_models = {
      LinearAlign::ScoreModel::LinearTurboFold,
      LinearAlign::ScoreModel::CONTRAlign,
      LinearAlign::ScoreModel::ProbConsRNA};
  for (size_t model_index = 0; model_index < score_models.size(); ++model_index) {
    LinearAlign selected_align(0.0001f, 100, score_models[model_index]);
    selected_align.calculate(sequence, "GGGAAUCCC", mp);
    size_t selected_match_count = 0;
    for (const auto& row : mp) {
      for (const auto& [j, probability] : row) {
        if (j >= sequence.size() || !std::isfinite(probability) ||
            probability <= 0.0f || probability > 1.0001f)
          return 1;
        posterior_sums[model_index] += probability;
        ++selected_match_count;
      }
    }
    if (selected_match_count == 0)
      return 1;
  }
  if (std::abs(posterior_sums[0] - posterior_sums[1]) < 1e-6 ||
      std::abs(posterior_sums[0] - posterior_sums[2]) < 1e-6 ||
      std::abs(posterior_sums[1] - posterior_sums[2]) < 1e-6)
    return 1;

  return 0;
}
