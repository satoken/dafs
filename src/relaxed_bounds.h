#ifndef DAFS_RELAXED_BOUNDS_H
#define DAFS_RELAXED_BOUNDS_H

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <utility>
#include <vector>

namespace RelaxedBounds
{
struct StructureBound
{
  float all = 0.0f;
  float left = 0.0f;
  float right = 0.0f;
  float incident = 0.0f;
  float cardinality = 0.0f;
  float top_k = 0.0f;
  float vertex_cover = 0.0f;
  float best = 0.0f;
};

struct AlignmentBound
{
  float row = 0.0f;
  float column = 0.0f;
  float best = 0.0f;
};

struct MonotoneAlignmentSolution
{
  float value = 0.0f;
  std::vector<unsigned> selected;
};

struct WeightedEdge
{
  unsigned first;
  unsigned second;
  float value;
};

inline float round_up_to_float(double value)
{
  if (!(value > 0.0))
    return 0.0f;
  float rounded = static_cast<float>(value);
  if (std::isfinite(rounded) && static_cast<double>(rounded) < value)
    rounded = std::nextafter(rounded, std::numeric_limits<float>::infinity());
  return rounded;
}

inline std::uint32_t nonnegative_float_key(float value)
{
  static_assert(sizeof(float) == sizeof(std::uint32_t),
                "32-bit IEEE-style float required");
  std::uint32_t key = 0;
  std::memcpy(&key, &value, sizeof(key));
  return key;
}

// Sum the K largest non-negative float values using a fixed four-pass radix
// sort.  This is worst-case O(m), unlike comparison selection/sorting.
inline double sum_largest_k(std::vector<float> values, std::size_t k)
{
  if (values.empty() || k == 0)
    return 0.0;
  k = std::min(k, values.size());
  std::vector<float> scratch(values.size());
  for (unsigned shift = 0; shift < 32; shift += 8) {
    std::array<std::size_t, 256> counts{};
    for (const float value : values)
      ++counts[(nonnegative_float_key(value) >> shift) & 0xffu];
    std::array<std::size_t, 256> offsets{};
    for (std::size_t bucket = 1; bucket < offsets.size(); ++bucket)
      offsets[bucket] = offsets[bucket-1] + counts[bucket-1];
    for (const float value : values) {
      const std::size_t bucket =
          (nonnegative_float_key(value) >> shift) & 0xffu;
      scratch[offsets[bucket]++] = value;
    }
    values.swap(scratch);
  }

  double result = 0.0;
  for (std::size_t index = values.size() - k; index < values.size(); ++index)
    result += static_cast<double>(values[index]);
  return result;
}

// Produce a feasible dual solution of the fractional matching relaxation:
// min sum_i u_i, subject to u_i+u_j >= score(i,j), u_i >= 0.  Starting from
// the half-maximum cover is feasible.  Each coordinate tightening preserves
// every constraint, so every recorded objective is a rigorous upper bound.
// A fixed number of sweeps keeps the work linear in the sparse edge support.
inline double fractional_vertex_cover(
    const std::vector<WeightedEdge>& edges, unsigned length,
    const std::vector<float>& best_incident)
{
  std::vector<std::vector<std::pair<unsigned, float>>> adjacency(length);
  for (const auto& edge : edges) {
    adjacency[edge.first].push_back({edge.second, edge.value});
    adjacency[edge.second].push_back({edge.first, edge.value});
  }

  std::vector<double> potential(length, 0.0);
  double objective = 0.0;
  for (unsigned i = 0; i < length; ++i) {
    potential[i] = 0.5 * static_cast<double>(best_incident[i]);
    objective += potential[i];
  }
  double best_objective = objective;
  constexpr unsigned sweeps = 4;
  for (unsigned sweep = 0; sweep < sweeps; ++sweep) {
    const bool reverse = (sweep & 1u) != 0;
    for (unsigned step = 0; step < length; ++step) {
      const unsigned i = reverse ? length - 1 - step : step;
      double tightened = 0.0;
      for (const auto& [j, value] : adjacency[i])
        tightened = std::max(
            tightened, static_cast<double>(value) - potential[j]);
      objective += tightened - potential[i];
      potential[i] = tightened;
    }
    best_objective = std::min(best_objective, objective);
  }
  return best_objective;
}

// Relax non-crossing, but retain necessary matching constraints: every base
// can be the left endpoint, right endpoint, or either endpoint of at most one
// selected pair.  Each expression is independently an upper bound, so their
// minimum is an upper bound as well.
template <typename Score>
StructureBound structure_bound(
    const std::vector<std::pair<unsigned, unsigned>>& support,
    unsigned length, Score score)
{
  std::vector<float> best_left(length, 0.0f);
  std::vector<float> best_right(length, 0.0f);
  std::vector<float> best_incident(length, 0.0f);
  std::vector<WeightedEdge> edges;
  std::vector<float> positive_values;
  edges.reserve(support.size());
  positive_values.reserve(support.size());
  double all = 0.0;
  float best_pair = 0.0f;

  for (const auto& [i, j] : support) {
    // Both exact and beam Nussinov decoders use a minimum hairpin loop
    // length of two, hence j-i must exceed two.
    if (i >= length || j >= length || j <= i + 2)
      continue;
    const float value = std::max(0.0f, score(i, j));
    if (value == 0.0f)
      continue;
    edges.push_back({i, j, value});
    positive_values.push_back(value);
    all += static_cast<double>(value);
    best_pair = std::max(best_pair, value);
    best_left[i] = std::max(best_left[i], value);
    best_right[j] = std::max(best_right[j], value);
    best_incident[i] = std::max(best_incident[i], value);
    best_incident[j] = std::max(best_incident[j], value);
  }

  double left = 0.0;
  double right = 0.0;
  double degree_twice = 0.0;
  for (unsigned i = 0; i < length; ++i) {
    left += static_cast<double>(best_left[i]);
    right += static_cast<double>(best_right[i]);
    degree_twice += static_cast<double>(best_incident[i]);
  }
  // With a minimum hairpin loop length of two, a non-empty non-crossing
  // structure with K pairs needs at least 2K+2 positions.
  const unsigned max_pairs = length > 2 ? (length - 2) / 2 : 0;
  const double cardinality =
      static_cast<double>(max_pairs) * static_cast<double>(best_pair);
  const double top_k = sum_largest_k(std::move(positive_values), max_pairs);
  const double vertex_cover =
      fractional_vertex_cover(edges, length, best_incident);

  StructureBound result;
  result.all = round_up_to_float(all);
  result.left = round_up_to_float(left);
  result.right = round_up_to_float(right);
  result.incident = round_up_to_float(0.5 * degree_twice);
  result.cardinality = round_up_to_float(cardinality);
  result.top_k = round_up_to_float(top_k);
  result.vertex_cover = round_up_to_float(vertex_cover);
  result.best = std::min({result.all, result.left, result.right,
                          result.incident, result.cardinality, result.top_k,
                          result.vertex_cover});
  return result;
}

template <typename Score>
float structure(const std::vector<std::pair<unsigned, unsigned>>& support,
                unsigned length, Score score)
{
  return structure_bound(support, length, score).best;
}

// Exact optimizer for the convex left-endpoint relaxation.  The selected
// edges are a subgradient of the returned support function.
template <typename Score>
float structure_left_solution(
    const std::vector<std::pair<unsigned, unsigned>>& support,
    unsigned length, Score score, std::vector<unsigned>& selected,
    unsigned minimum_pair_span = 3)
{
  selected.assign(length, std::numeric_limits<unsigned>::max());
  std::vector<float> best(length, 0.0f);
  for (const auto& [i, j] : support) {
    if (i >= length || j >= length || j <= i ||
        j - i < minimum_pair_span)
      continue;
    const float value = score(i, j);
    if (value > best[i]) {
      best[i] = value;
      selected[i] = j;
    }
  }
  double total = 0.0;
  for (const float value : best)
    total += static_cast<double>(value);
  return round_up_to_float(total);
}

// A monotone alignment is in particular a bipartite matching.  Dropping
// monotonicity while retaining capacity one on either side gives two cheap
// upper bounds; take the tighter one.
template <typename Score>
AlignmentBound alignment_bound(
    const std::vector<std::pair<unsigned, unsigned>>& support,
    unsigned rows, unsigned columns, Score score)
{
  std::vector<float> best_row(rows, 0.0f);
  std::vector<float> best_column(columns, 0.0f);
  for (const auto& [i, k] : support) {
    if (i >= rows || k >= columns)
      continue;
    const float value = std::max(0.0f, score(i, k));
    best_row[i] = std::max(best_row[i], value);
    best_column[k] = std::max(best_column[k], value);
  }

  double row = 0.0;
  double column = 0.0;
  for (const float value : best_row)
    row += static_cast<double>(value);
  for (const float value : best_column)
    column += static_cast<double>(value);
  AlignmentBound result;
  result.row = round_up_to_float(row);
  result.column = round_up_to_float(column);
  result.best = std::min(result.row, result.column);
  return result;
}

template <typename Score>
float alignment(const std::vector<std::pair<unsigned, unsigned>>& support,
                unsigned rows, unsigned columns, Score score)
{
  return alignment_bound(support, rows, columns, score).best;
}

// Exact optimizer for the positive-weight monotone matching relaxation.
// Candidate edges are processed row by row.  A Fenwick tree stores the best
// predecessor whose last column is strictly smaller than the current column;
// delaying all updates from one row until that row has been evaluated prevents
// two edges from the same row being chained together.  Invalid and duplicate
// support coordinates are ignored/merged, so callers may pass a support union
// assembled from more than one sparse source.
//
// The returned mapping is a subgradient of this convex maximum-of-affines
// function.  Its value is accumulated in double precision and rounded
// outward to float for use as a certified upper bound.  The sort is
// O(m log m), which is O(m log L) for an L-by-L sparse matrix, and the DP
// storage is O(m + L).
template <typename Score>
MonotoneAlignmentSolution alignment_monotone_solution(
    const std::vector<std::pair<unsigned, unsigned>>& support,
    unsigned rows, unsigned columns, Score score)
{
  constexpr std::size_t no_edge = std::numeric_limits<std::size_t>::max();
  const unsigned missing = std::numeric_limits<unsigned>::max();

  MonotoneAlignmentSolution result;
  result.selected.assign(rows, missing);

  struct Candidate
  {
    unsigned row;
    unsigned column;
    float score;
    double total = 0.0;
    std::size_t previous = std::numeric_limits<std::size_t>::max();
  };

  std::vector<Candidate> candidates;
  candidates.reserve(support.size());
  for (const auto& [row, column] : support) {
    if (row >= rows || column >= columns)
      continue;
    const float value = score(row, column);
    // A non-positive edge can never improve a maximum matching that may
    // leave rows unmatched.  The comparison also safely ignores NaN input.
    if (!(value > 0.0f))
      continue;
    candidates.push_back({row, column, value});
  }

  std::sort(candidates.begin(), candidates.end(),
            [](const Candidate& lhs, const Candidate& rhs) {
              if (lhs.row != rhs.row)
                return lhs.row < rhs.row;
              if (lhs.column != rhs.column)
                return lhs.column < rhs.column;
              return lhs.score > rhs.score;
            });
  // The support is normally already unique, but canonicalizing here keeps
  // malformed test/support unions deterministic and does not change the
  // sparse asymptotic storage bound.
  std::size_t unique_size = 0;
  for (const Candidate& candidate : candidates) {
    if (unique_size != 0 &&
        candidates[unique_size - 1].row == candidate.row &&
        candidates[unique_size - 1].column == candidate.column)
      continue;
    candidates[unique_size++] = candidate;
  }
  candidates.resize(unique_size);

  struct State
  {
    double value;
    std::size_t edge;
  };
  const State empty{0.0, no_edge};
  const auto better = [&](const State& lhs, const State& rhs) {
    if (lhs.value > rhs.value)
      return lhs;
    if (rhs.value > lhs.value)
      return rhs;
    if (lhs.edge == no_edge)
      return rhs;
    if (rhs.edge == no_edge)
      return lhs;
    // Equal-valued states with a smaller last column dominate the other state
    // for every future strict-prefix query.  The edge index is a final stable
    // tie breaker after the canonical row/column ordering above.
    if (candidates[lhs.edge].column != candidates[rhs.edge].column)
      return candidates[lhs.edge].column < candidates[rhs.edge].column
           ? lhs : rhs;
    return lhs.edge < rhs.edge ? lhs : rhs;
  };

  std::vector<State> fenwick(columns + 1, empty);
  const auto prefix_best = [&](unsigned count) {
    State best = empty;
    std::size_t index = count;
    while (index != 0) {
      best = better(best, fenwick[index]);
      index &= index - 1;
    }
    return best;
  };

  std::size_t begin = 0;
  while (begin < candidates.size()) {
    std::size_t end = begin + 1;
    while (end < candidates.size() &&
           candidates[end].row == candidates[begin].row)
      ++end;

    // Compute the whole row before updating the tree.  Otherwise two edges
    // from this row could incorrectly form one matching path.
    for (std::size_t edge = begin; edge < end; ++edge) {
      const State predecessor = prefix_best(candidates[edge].column);
      candidates[edge].previous = predecessor.edge;
      candidates[edge].total = predecessor.value +
          static_cast<double>(candidates[edge].score);
    }
    for (std::size_t edge = begin; edge < end; ++edge) {
      State candidate_state{candidates[edge].total, edge};
      std::size_t index = static_cast<std::size_t>(candidates[edge].column) + 1;
      while (index <= columns) {
        fenwick[index] = better(fenwick[index], candidate_state);
        index += index & (~index + 1);
      }
    }
    begin = end;
  }

  const State best = prefix_best(columns);
  for (std::size_t edge = best.edge; edge != no_edge;
       edge = candidates[edge].previous)
    result.selected[candidates[edge].row] = candidates[edge].column;
  result.value = round_up_to_float(best.value);
  return result;
}

// Convenience form matching the existing relaxed-bound helpers.
template <typename Score>
float alignment_monotone_solution(
    const std::vector<std::pair<unsigned, unsigned>>& support,
    unsigned rows, unsigned columns, Score score,
    std::vector<unsigned>& selected)
{
  MonotoneAlignmentSolution result = alignment_monotone_solution(
      support, rows, columns, score);
  selected = std::move(result.selected);
  return result.value;
}


// Exact optimizer for the convex row-capacity alignment relaxation.
template <typename Score>
float alignment_row_solution(
    const std::vector<std::pair<unsigned, unsigned>>& support,
    unsigned rows, unsigned columns, Score score,
    std::vector<unsigned>& selected)
{
  selected.assign(rows, std::numeric_limits<unsigned>::max());
  std::vector<float> best(rows, 0.0f);
  for (const auto& [i, k] : support) {
    if (i >= rows || k >= columns)
      continue;
    const float value = score(i, k);
    if (value > best[i]) {
      best[i] = value;
      selected[i] = k;
    }
  }
  double total = 0.0;
  for (const float value : best)
    total += static_cast<double>(value);
  return round_up_to_float(total);
}
}

#endif
