#pragma once

// Alignment energy adapter for LinFold's partition grammar.  ViennaRNA
// supplies only constant-time nearest-neighbor energy primitives.  Search,
// beam pruning, partition sums and outside marginals stay in LinFold.

#include <algorithm>
#include <array>
#include <cstdlib>
#include <memory>
#include <cctype>
#include <cstdint>
#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

extern "C" {
#include <ViennaRNA/model.h>
#include <ViennaRNA/params/basic.h>
#include <ViennaRNA/eval/hairpin.h>
#include <ViennaRNA/eval/internal.h>
#include <ViennaRNA/eval/multibranch.h>
#include <ViennaRNA/eval/exterior.h>
}

class ProfileNearestNeighbor
{
public:
  using ScoreType = double;

  enum class EnergyMode {
    Legacy,
    RNAalifold
  };

  // RT at the temperature used by the existing profile folding path.
  static constexpr double RT = 0.6163207755;

  explicit ProfileNearestNeighbor(
      const std::vector<std::string>& aligned,
      EnergyMode energy_mode = EnergyMode::Legacy)
    : length_(aligned.empty() ? 0 : aligned.front().size()),
      occupancy_prefix_(length_ + 1, 0.0),
      energy_mode_(energy_mode)
  {
    if (aligned.empty())
      throw std::invalid_argument("empty folding profile");

    vrna_md_t model;
    vrna_md_set_default(&model);
    model.temperature = 37.0;
    model.dangles = 2;
    // RNAalifold's normal Vienna model enables sequence-dependent special
    // tri-/tetra-/hexaloop terms.  The legacy adapter deliberately disabled
    // them because it passed no loop sequence to vrna_E_hairpin().
    model.special_hp = energy_mode_ == EnergyMode::RNAalifold ? 1 : 0;
    parameters_.reset(vrna_params(&model));
    if (!parameters_) throw std::bad_alloc();

    rows_.reserve(aligned.size());
    if (aligned.size() > 2 && aligned.size() <= 64)
      column_masks_.resize(length_ + 1);
    for (const auto& input : aligned) {
      if (input.size() != length_)
        throw std::invalid_argument("unequal folding profile lengths");

      Row row;
      row.codes.resize(length_ + 1, kUnknown);
      row.prefix.resize(length_ + 1, 0);
      row.previous.resize(length_ + 1, kNoNeighbor);
      row.next.resize(length_ + 1, kNoNeighbor);
      row.nongap_sequence.reserve(length_);
      for (size_t column = 1; column <= length_; ++column) {
        row.codes[column] = nucleotide_code(input[column - 1]);
        if (!column_masks_.empty() && row.codes[column] != kUnknown)
          column_masks_[column][row.codes[column]] |=
              uint64_t{1} << rows_.size();
        row.prefix[column] = row.prefix[column - 1] +
                            (row.codes[column] == kGap ? 0u : 1u);
        if (row.codes[column] != kGap)
          row.nongap_sequence.push_back(sequence_code(input[column - 1]));
        occupancy_prefix_[column] +=
            row.codes[column] == kGap ? 0.0 : 1.0;
      }
      int previous = kNoNeighbor;
      for (size_t column = 1; column <= length_; ++column) {
        row.previous[column] = previous;
        if (row.codes[column] != kGap)
          previous = row.codes[column];
      }
      int next = kNoNeighbor;
      for (size_t column = length_; column > 0; --column) {
        row.next[column] = next;
        if (row.codes[column] != kGap)
          next = row.codes[column];
      }
      rows_.push_back(std::move(row));
    }

    average_denominator_ = 100.0 * RT * static_cast<double>(rows_.size());
    const double row_count = static_cast<double>(rows_.size());
    for (size_t column = 1; column <= length_; ++column)
      occupancy_prefix_[column] += occupancy_prefix_[column - 1];
    for (double& value : occupancy_prefix_)
      value /= row_count;

    // One and two rows have direct predicates; 3..64 use column bit masks.
    // Larger profiles retain a bounded, length-proportional lazy cache.
    if (rows_.size() > 64)
      pair_eligibility_cache_limit_ =
          std::max<size_t>(4096, length_ * 64);

    if (rows_.size() == 2) {
      for (size_t left = 0; left < 36; ++left)
        for (size_t right = 0; right < 36; ++right) {
          const int first = pair_type(left / 6, right / 6);
          const int second = pair_type(left % 6, right % 6);
          const unsigned incompatible =
              static_cast<unsigned>(first == kIncompatiblePair) +
              static_cast<unsigned>(second == kIncompatiblePair);
          const unsigned double_gap =
              static_cast<unsigned>(first == kDoubleGapPair) +
              static_cast<unsigned>(second == kDoubleGapPair);
          const bool canonical =
              (first >= kFirstCanonicalPair && first <= kLastCanonicalPair) ||
              (second >= kFirstCanonicalPair && second <= kLastCanonicalPair);
          two_row_eligible_[left * 36 + right] =
              canonical && 2 * incompatible + double_gap < 2;
        }
    }

    // Pair eligibility is evaluated lazily.  The LinearAlifold-style
    // candidate scan usually asks for only beam-bounded pairs; precomputing
    // every column pair would make the profile path quadratic in memory.
  }

  // Cumulative mean nongap occupancy through each 1-based alignment column.
  // occupancy_prefix()[0] is zero.
  const std::vector<double>& occupancy_prefix() const noexcept
  {
    return occupancy_prefix_;
  }

  // Convenience overload for callers that need one prefix entry.
  double occupancy_prefix(size_t column) const
  {
    return occupancy_prefix_.at(column);
  }

  size_t length() const noexcept { return length_; }
  size_t row_count() const noexcept { return rows_.size(); }

  // For a two-row profile, eligibility depends only on these two column
  // codes.  Equal signatures therefore share the same eligible right columns.
  // The caller supplies a valid 1-based column and a two-row profile.
  uint8_t two_row_signature(size_t column) const noexcept
  {
    return static_cast<uint8_t>(6 * rows_[0].codes[column] +
                                rows_[1].codes[column]);
  }

  // A consensus pair is searchable when the RNAalifold majority rule holds:
  // 2 * incompatible + double-gap < number of rows, and at least one row
  // contains a canonical pair.
  bool can_pair(size_t i, size_t j) const
  {
    if (i == 0 || i >= j || j > length_)
      return false;

    if (rows_.size() == 1) {
      const int type = pair_type(rows_.front().codes[i], rows_.front().codes[j]);
      return type >= kFirstCanonicalPair && type <= kLastCanonicalPair;
    }
    if (rows_.size() == 2)
      return two_row_eligible_[
          static_cast<size_t>(two_row_signature(i)) * 36 +
          two_row_signature(j)] != 0;
    if (!column_masks_.empty()) {
      const auto counts = pair_type_counts(i, j);
      const unsigned canonical = counts[1] + counts[2] + counts[3] +
                                 counts[4] + counts[5] + counts[6];
      return canonical != 0 &&
             2 * counts[0] + counts[7] < rows_.size();
    }
    if (pair_eligibility_cache_limit_ == 0)
      return calculate_pair_eligibility(i, j);
    const auto key = std::make_pair(i, j);
    const auto found = pair_eligibility_cache_.find(key);
    if (found != pair_eligibility_cache_.end())
      return found->second != 0;
    const bool eligible = calculate_pair_eligibility(i, j);
    if (pair_eligibility_cache_.size() < pair_eligibility_cache_limit_)
      pair_eligibility_cache_.emplace(key, static_cast<uint8_t>(eligible));
    return eligible;
  }

  // RNAalifold's six oriented canonical pair types, incompatible pairs and
  // double gaps.  The bit masks partition rows exactly, including unknowns.
  std::array<unsigned, 8> pair_type_counts(size_t i, size_t j) const
  {
    std::array<unsigned, 8> counts{};
    if (!column_masks_.empty()) {
      const auto& left = column_masks_[i];
      const auto& right = column_masks_[j];
      const auto count = [](uint64_t mask) {
        return static_cast<unsigned>(__builtin_popcountll(mask));
      };
      counts[1] = count(left[2] & right[3]); // CG
      counts[2] = count(left[3] & right[2]); // GC
      counts[3] = count(left[3] & right[4]); // GU
      counts[4] = count(left[4] & right[3]); // UG
      counts[5] = count(left[1] & right[4]); // AU
      counts[6] = count(left[4] & right[1]); // UA
      counts[7] = count(left[5] & right[5]);
      counts[0] = static_cast<unsigned>(rows_.size()) -
          (counts[1] + counts[2] + counts[3] + counts[4] +
           counts[5] + counts[6] + counts[7]);
    } else {
      for (const Row& row : rows_)
        ++counts[pair_type(row.codes[i], row.codes[j])];
    }
    return counts;
  }

private:
  bool calculate_pair_eligibility(size_t i, size_t j) const
  {

    size_t incompatible = 0;
    size_t double_gap = 0;
    bool canonical_pair = false;
    for (const Row& row : rows_) {
      const int type = pair_type(row.codes[i], row.codes[j]);
      if (type == kDoubleGapPair)
        ++double_gap;
      else if (type == kIncompatiblePair)
        ++incompatible;
      else if (type >= kFirstCanonicalPair && type <= kLastCanonicalPair)
        canonical_pair = true;
    }
    // RNAalifold suppresses a pair when counterexamples plus double-gaps
    // reach half of the profile (pscore == NONE).  Keep the same strict
    // eligibility condition here so the linear and nonlinear profile paths
    // expose the same pair candidates.
    return canonical_pair &&
           2 * incompatible + double_gap < rows_.size();
  }

public:
  ScoreType score_hairpin(size_t i, size_t j) const
  {
    return average([&](const Row& row) -> double {
      const size_t size = count_between(row, i, j);
      // LinearAlifold assigns 6 kcal/mol to short projected hairpins.
      if (size < 3) return 600.0;
      const int type = energy_pair_type(row.codes[i], row.codes[j]);
      // Generic triloops have no terminal mismatch term, even when special
      // sequence-dependent hairpin bonuses are disabled.
      if (size == 3 && energy_mode_ == EnergyMode::Legacy)
        return parameters_->hairpin[3] + (type > 2 ? parameters_->TerminalAU : 0);
      const std::string* loop_sequence = nullptr;
      std::string sequence;
      if (energy_mode_ == EnergyMode::RNAalifold && size < 7 &&
          type >= kFirstCanonicalPair && type <= kLastCanonicalPair) {
        sequence = hairpin_sequence(row, i, size);
        loop_sequence = &sequence;
      }
      return vrna_E_hairpin(size, energy_pair_type(row.codes[i], row.codes[j]),
          table_code(next_code(row, i)), table_code(previous_code(row, j)),
          loop_sequence ? loop_sequence->c_str() : nullptr, parameters_.get());
    });
  }

  ScoreType score_single_loop(size_t i, size_t j, size_t k, size_t l) const
  {
    return average([&](const Row& row) { return single_loop_energy(row, i, j, k, l); });
  }

  ScoreType score_helix(size_t i, size_t j, size_t m) const
  {
    if (m <= 1) return -0.0;
    return average([&](const Row& row) {
      double energy = 0;
      for (size_t step=0; step+1<m; ++step)
        energy += single_loop_energy(row, i+step, j-step, i+step+1, j-step-1);
      return energy;
    });
  }

  // Called with increasing m for one inner pair.  Store each row's newly
  // exposed outer stack once, but sum outer-to-inner for every m exactly as
  // score_helix() does; this preserves floating-point accumulation order.
  ScoreType score_helix_extend(size_t inner_i, size_t inner_j, size_t m,
                                std::vector<double>& scratch) const
  {
    const size_t n = rows_.size();
    const size_t step = m - 2;
    if (scratch.size() < (step + 1) * n)
      scratch.resize((step + 1) * n);
    for (size_t row = 0; row < n; ++row)
      scratch[step * n + row] = single_loop_energy(rows_[row],
          inner_i - (m - 1), inner_j + (m - 1),
          inner_i - (m - 2), inner_j + (m - 2));
    double total = 0;
    for (size_t row = 0; row < n; ++row) {
      double energy = 0;
      for (size_t edge = step + 1; edge > 0; --edge)
        energy += scratch[(edge - 1) * n + row];
      total += energy;
    }
    return -total / average_denominator_;
  }

  ScoreType score_multi_loop(size_t i, size_t j) const
  {
    return average([&](const Row& row) {
      return parameters_->MLclosing + vrna_E_multibranch_stem(
          energy_pair_type(row.codes[j], row.codes[i]),
          previous_code(row, j), next_code(row, i), parameters_.get());
    });
  }

  ScoreType score_multi_paired(size_t i, size_t j) const
  {
    return average([&](const Row& row) {
      return vrna_E_multibranch_stem(energy_pair_type(row.codes[i], row.codes[j]),
          previous_code(row, i), next_code(row, j), parameters_.get());
    });
  }

  ScoreType score_multi_unpaired(size_t i, size_t j) const
  {
    if (i > j) return 0;
    if (parameters_->MLbase == 0) return -0.0;
    return average([&](const Row& row) {
      return static_cast<double>(parameters_->MLbase) * count_inclusive(row, i, j);
    });
  }

  ScoreType score_external_zero() const { return 0; }
  ScoreType score_external_paired(size_t i, size_t j) const
  {
    return average([&](const Row& row) {
      return vrna_E_exterior_stem(energy_pair_type(row.codes[i], row.codes[j]),
          previous_code(row, i), next_code(row, j), parameters_.get());
    });
  }
  ScoreType score_external_unpaired(size_t, size_t) const { return 0; }

  // Profile folding has no trainable thermodynamic parameters.  Keep the
  // count API that LinFold parameter models expose so generic callers can use
  // the adapter without special cases.
  void count_hairpin(size_t, size_t, ScoreType) {}
  void count_single_loop(size_t, size_t, size_t, size_t, ScoreType) {}
  void count_helix(size_t, size_t, size_t, ScoreType) {}
  void count_multi_loop(size_t, size_t, ScoreType) {}
  void count_multi_paired(size_t, size_t, ScoreType) {}
  void count_multi_unpaired(size_t, size_t, ScoreType) {}
  void count_external_zero(ScoreType) {}
  void count_external_paired(size_t, size_t, ScoreType) {}
  void count_external_unpaired(size_t, size_t, ScoreType) {}

private:
  static constexpr int kUnknown = 0;
  static constexpr int kGap = 5;
  static constexpr int kNoNeighbor = -1;
  static constexpr int kIncompatiblePair = 0;
  static constexpr int kDoubleGapPair = 7;
  static constexpr int kFirstCanonicalPair = 1;
  static constexpr int kLastCanonicalPair = 6;

  struct Row {
    // Codes are indexed by 1-based alignment column.  A gap is kept as 5 so
    // that an invalid pair can be represented by energy table type 7.
    std::vector<int> codes;
    // Number of nongap symbols through each alignment column.
    std::vector<size_t> prefix;
    // Nearest nongap code strictly before/after each alignment column.
    std::vector<int> previous;
    std::vector<int> next;
    // The ungapped row is used only for ViennaRNA's short, sequence-specific
    // hairpin lookup.  Keeping it alongside the prefix table avoids scanning
    // alignment gaps when a beam candidate is evaluated.
    std::string nongap_sequence;
  };

  struct PairHash {
    size_t operator()(const std::pair<size_t, size_t>& pair) const noexcept
    {
      size_t value = std::hash<size_t>{}(pair.first);
      value ^= std::hash<size_t>{}(pair.second) +
          static_cast<size_t>(0x9e3779b97f4a7c15ULL) +
          (value << 6) + (value >> 2);
      return value;
    }
  };

  static int nucleotide_code(char nucleotide)
  {
    switch (std::tolower(static_cast<unsigned char>(nucleotide))) {
    case 'a': return 1;
    case 'c': return 2;
    case 'g': return 3;
    case 'u':
    case 't': return 4;
    case '-': return kGap;
    default: return kUnknown;
    }
  }

  static char sequence_code(char nucleotide)
  {
    switch (std::tolower(static_cast<unsigned char>(nucleotide))) {
    case 'a': return 'A';
    case 'c': return 'C';
    case 'g': return 'G';
    case 'u':
    case 't': return 'U';
    default: return 'N';
    }
  }

  // This is the same orientation and numbering as TurnerNearestNeighbor's
  // private complement_pair table: CG=1, GC=2, GU=3, UG=4, AU=5, UA=6.
  static int complement_pair(int left, int right)
  {
    static constexpr int table[5][5] = {
      {0, 0, 0, 0, 0},
      {0, 0, 0, 0, 5},
      {0, 0, 0, 1, 0},
      {0, 0, 2, 0, 3},
      {0, 6, 0, 4, 0},
    };
    return (left >= 0 && left < 5 && right >= 0 && right < 5)
        ? table[left][right]
        : kIncompatiblePair;
  }

  static int pair_type(int left, int right)
  {
    if (left == kGap && right == kGap)
      return kDoubleGapPair;
    if (left == kGap || right == kGap)
      return kIncompatiblePair;
    return complement_pair(left, right);
  }

  // RNAalifold's thermodynamic tables use type 7 for every pair that is not
  // one of the six canonical Turner pairs.  Eligibility keeps the separate
  // type-0 incompatible and type-7 double-gap counts above; energy scoring
  // deliberately follows the upstream type-7 fallback.
  static int energy_pair_type(int left, int right)
  {
    // Profile rows only contain codes 0..5 here.  Keep a fallback for
    // defensive callers while avoiding the second table lookup used by
    // complement_pair() for every internal-loop stack.
    static constexpr int table[6][6] = {
      {7, 7, 7, 7, 7, 7},
      {7, 7, 7, 7, 5, 7},
      {7, 7, 7, 1, 7, 7},
      {7, 7, 2, 7, 3, 7},
      {7, 6, 7, 4, 7, 7},
      {7, 7, 7, 7, 7, 7},
    };
    if (static_cast<unsigned>(left) < 6 &&
        static_cast<unsigned>(right) < 6)
      return table[left][right];
    return kDoubleGapPair;
  }

  static int table_code(int code)
  {
    return code >= 0 && code < 5 ? code : kUnknown;
  }

  template <class F> ScoreType average(F&& energy) const
  {
    double total = 0;
    for (const Row& row : rows_) total += energy(row);
    return -total / average_denominator_;
  }

  static size_t count_between(const Row& row, size_t left, size_t right)
  {
    if (left >= right || left >= row.prefix.size() || right == 0)
      return 0;
    if (right - 1 >= row.prefix.size())
      return 0;
    return row.prefix[right - 1] - row.prefix[left];
  }

  static size_t count_inclusive(const Row& row, size_t left, size_t right)
  {
    if (left == 0 || left > right || right >= row.prefix.size())
      return 0;
    return row.prefix[right] - row.prefix[left - 1];
  }

  static std::string hairpin_sequence(const Row& row, size_t left,
                                      size_t loop_size)
  {
    // prefix[left - 1] is the zero-based start of the closing pair in the
    // ungapped row.  The sequence passed to vrna_E_hairpin() contains both
    // closing nucleotides and all loop nucleotides.
    const size_t start = row.prefix[left - 1];
    return row.nongap_sequence.substr(start, loop_size + 2);
  }

  static int code_at(const Row& row, size_t column)
  {
    return column < row.codes.size() ? row.codes[column] : kUnknown;
  }

  static int next_code(const Row& row, size_t column)
  {
    return column < row.next.size() ? row.next[column] : kNoNeighbor;
  }

  static int previous_code(const Row& row, size_t column)
  {
    return column < row.previous.size() ? row.previous[column] : kNoNeighbor;
  }

  double single_loop_energy(const Row& row, size_t i, size_t j,
                            size_t k, size_t l) const
  {
    const size_t left_unpaired = count_between(row, i, k);
    const size_t right_unpaired = count_between(row, l, j);
    const int outer_type = energy_pair_type(row.codes[i], row.codes[j]);
    const int inner_type = energy_pair_type(row.codes[l], row.codes[k]);
    // With the default model used by this adapter, a zero-by-zero internal
    // loop is exactly the nearest-neighbor stack entry.  The guards retain
    // ViennaRNA's salt and no-GU-closure behavior if model details change.
    if (left_unpaired == 0 && right_unpaired == 0 &&
        parameters_->SaltStack == 0 && !parameters_->model_details.noGUclosure)
      return parameters_->stack[outer_type][inner_type];
    const int si1 = table_code(next_code(row, i));
    const int sj1 = table_code(previous_code(row, j));
    const int sp1 = table_code(previous_code(row, k));
    const int sq1 = table_code(next_code(row, l));
    return vrna_E_internal(left_unpaired, right_unpaired, outer_type,
        inner_type, si1, sj1, sp1, sq1, parameters_.get());
  }

  struct ParameterDeleter {
    void operator()(vrna_param_t* parameters) const { std::free(parameters); }
  };
  std::unique_ptr<vrna_param_t, ParameterDeleter> parameters_;
  EnergyMode energy_mode_;
  std::vector<Row> rows_;
  std::vector<std::array<uint64_t, 6>> column_masks_;
  size_t length_;
  std::vector<double> occupancy_prefix_;
  mutable std::unordered_map<std::pair<size_t, size_t>, uint8_t, PairHash>
      pair_eligibility_cache_;
  std::array<uint8_t, 36 * 36> two_row_eligible_{};
  size_t pair_eligibility_cache_limit_ = 0;
  double average_denominator_ = 0.0;
};
