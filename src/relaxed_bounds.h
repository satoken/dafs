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
    unsigned length, Score score, std::vector<unsigned>& selected)
{
  selected.assign(length, std::numeric_limits<unsigned>::max());
  std::vector<float> best(length, 0.0f);
  for (const auto& [i, j] : support) {
    if (i >= length || j >= length || j <= i + 2)
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
