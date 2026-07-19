#ifndef DAFS_RELAXED_BOUNDS_H
#define DAFS_RELAXED_BOUNDS_H

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <utility>
#include <vector>

namespace RelaxedBounds
{
inline float round_up_to_float(double value)
{
  if (!(value > 0.0))
    return 0.0f;
  float rounded = static_cast<float>(value);
  if (std::isfinite(rounded) && static_cast<double>(rounded) < value)
    rounded = std::nextafter(rounded, std::numeric_limits<float>::infinity());
  return rounded;
}

// Relax non-crossing, but retain necessary matching constraints: every base
// can be the left endpoint, right endpoint, or either endpoint of at most one
// selected pair.  Each expression is independently an upper bound, so their
// minimum is an upper bound as well.
template <typename Score>
float structure(const std::vector<std::pair<unsigned, unsigned>>& support,
                unsigned length, Score score)
{
  std::vector<float> best_left(length, 0.0f);
  std::vector<float> best_right(length, 0.0f);
  std::vector<float> best_incident(length, 0.0f);
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
  // structure with K pairs needs at least 2K+2 positions.  Multiplying that
  // cardinality limit by the best edge is a guaranteed-linear substitute for
  // selecting the K largest weights.
  const unsigned max_pairs = length > 2 ? (length - 2) / 2 : 0;
  const double cardinality =
      static_cast<double>(max_pairs) * static_cast<double>(best_pair);
  return round_up_to_float(std::min(
      {all, left, right, 0.5 * degree_twice, cardinality}));
}

// A monotone alignment is in particular a bipartite matching.  Dropping
// monotonicity while retaining capacity one on either side gives two cheap
// upper bounds; take the tighter one.
template <typename Score>
float alignment(const std::vector<std::pair<unsigned, unsigned>>& support,
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
  return round_up_to_float(std::min(row, column));
}
}

#endif
