#include "relaxed_bounds.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <random>
#include <utility>
#include <vector>

namespace
{
using Edge = std::pair<unsigned, unsigned>;

bool close(float lhs, float rhs)
{
  return std::fabs(lhs-rhs) <= 1e-5f * std::max(1.0f, std::fabs(rhs));
}

bool valid_structure_subset(const std::vector<Edge>& edges, std::uint64_t mask)
{
  for (size_t a = 0; a < edges.size(); ++a) {
    if (!(mask & (std::uint64_t{1} << a)))
      continue;
    const auto [i, j] = edges[a];
    if (j <= i + 2)
      return false;
    for (size_t b = a + 1; b < edges.size(); ++b) {
      if (!(mask & (std::uint64_t{1} << b)))
        continue;
      const auto [k, l] = edges[b];
      if (i == k || i == l || j == k || j == l)
        return false;
      if ((i < k && k < j && j < l) || (k < i && i < l && l < j))
        return false;
    }
  }
  return true;
}

bool valid_alignment_subset(const std::vector<Edge>& edges, std::uint64_t mask)
{
  for (size_t a = 0; a < edges.size(); ++a) {
    if (!(mask & (std::uint64_t{1} << a)))
      continue;
    for (size_t b = a + 1; b < edges.size(); ++b) {
      if (!(mask & (std::uint64_t{1} << b)))
        continue;
      const auto [i, k] = edges[a];
      const auto [j, l] = edges[b];
      if (i == j || k == l || (i < j) != (k < l))
        return false;
    }
  }
  return true;
}

bool valid_monotone_subset(const std::vector<Edge>& edges, std::uint64_t mask)
{
  for (size_t a = 0; a < edges.size(); ++a) {
    if (!(mask & (std::uint64_t{1} << a)))
      continue;
    for (size_t b = a + 1; b < edges.size(); ++b) {
      if (!(mask & (std::uint64_t{1} << b)))
        continue;
      const auto [i, k] = edges[a];
      const auto [j, l] = edges[b];
      if (i == j || k == l || (i < j && k >= l) ||
          (j < i && l >= k))
        return false;
    }
  }
  return true;
}

template <typename Valid>
float exact_value(const std::vector<Edge>& edges,
                  const std::vector<float>& scores, Valid valid)
{
  float best = 0.0f;
  const std::uint64_t end = std::uint64_t{1} << edges.size();
  for (std::uint64_t mask = 0; mask < end; ++mask) {
    if (!valid(edges, mask))
      continue;
    float value = 0.0f;
    for (size_t e = 0; e < edges.size(); ++e)
      if (mask & (std::uint64_t{1} << e))
        value += std::max(0.0f, scores[e]);
    best = std::max(best, value);
  }
  return best;
}

template <typename Score>
bool check_monotone_mapping(const std::vector<unsigned>& selected,
                            unsigned columns, Score score, float value)
{
  unsigned previous_column = 0;
  bool have_previous = false;
  float mapped_value = 0.0f;
  for (unsigned row = 0; row < selected.size(); ++row) {
    const unsigned column = selected[row];
    if (column == std::numeric_limits<unsigned>::max())
      continue;
    if (column >= columns || (have_previous && column <= previous_column))
      return false;
    const float edge_value = score(row, column);
    if (!(edge_value > 0.0f))
      return false;
    mapped_value += edge_value;
    previous_column = column;
    have_previous = true;
  }
  return mapped_value <= value + 1e-5f * std::max(1.0f, std::fabs(value));
}

bool check_bound(float exact, float bound, float old_bound, const char* name)
{
  if (bound + 1e-5f < exact || bound > old_bound + 1e-5f) {
    std::cerr << name << ": exact=" << exact << " bound=" << bound
              << " old=" << old_bound << '\n';
    return false;
  }
  return true;
}

bool check_structure_components(float exact,
                                const RelaxedBounds::StructureBound& bound)
{
  const std::vector<std::pair<const char*, float>> components = {
      {"all", bound.all},
      {"left", bound.left},
      {"right", bound.right},
      {"incident", bound.incident},
      {"cardinality", bound.cardinality},
      {"top_k", bound.top_k},
      {"vertex_cover", bound.vertex_cover},
      {"best", bound.best}};
  for (const auto& [name, value] : components) {
    if (value + 1e-5f < exact) {
      std::cerr << name << " is not an upper bound: exact=" << exact
                << " bound=" << value << '\n';
      return false;
    }
  }
  return true;
}

bool check_alignment_components(
    float exact, const RelaxedBounds::AlignmentBound& bound)
{
  const std::vector<std::pair<const char*, float>> components = {
      {"row", bound.row}, {"column", bound.column}, {"best", bound.best}};
  for (const auto& [name, value] : components) {
    if (value + 1e-5f < exact) {
      std::cerr << name << " is not an upper bound: exact=" << exact
                << " bound=" << value << '\n';
      return false;
    }
  }
  return true;
}
}

int main()
{
  const std::vector<Edge> structure_edges = {
      {0, 6}, {0, 5}, {1, 5}, {1, 6}, {2, 6}, {2, 4}};
  const std::vector<float> structure_scores = {4.0f, 3.0f, 3.5f, 2.0f, 2.5f, 100.0f};
  const auto structure_score = [&](unsigned i, unsigned j) {
    const auto it = std::find(structure_edges.begin(), structure_edges.end(), Edge{i, j});
    return structure_scores[static_cast<size_t>(it - structure_edges.begin())];
  };
  const float exact_structure = exact_value(
      structure_edges, structure_scores, valid_structure_subset);
  const auto structure_details = RelaxedBounds::structure_bound(
      structure_edges, 7, structure_score);
  const float structure_bound = structure_details.best;
  std::vector<unsigned> relaxed_structure;
  const float left_solution = RelaxedBounds::structure_left_solution(
      structure_edges, 7, structure_score, relaxed_structure);
  if (!close(left_solution, structure_details.left))
    return 1;

  // The cached LinearNussinov support used by DAFS admits span-two pairs.
  // Its independent certified track must include the same edge set, while
  // the default exact-structure relaxation above retains span three.
  std::vector<unsigned> cached_relaxed_structure;
  const float cached_left_solution = RelaxedBounds::structure_left_solution(
      structure_edges, 7, structure_score, cached_relaxed_structure, 2);
  if (!close(cached_left_solution, 107.5f) ||
      cached_relaxed_structure[2] != 4)
    return 1;
  float old_structure_bound = 0.0f;
  for (size_t e = 0; e < structure_edges.size(); ++e)
    if (structure_edges[e].second > structure_edges[e].first + 2)
      old_structure_bound += std::max(0.0f, structure_scores[e]);
  if (!check_bound(exact_structure, structure_bound, old_structure_bound,
                   "structure"))
    return 1;
  if (!check_structure_components(exact_structure, structure_details))
    return 1;

  const std::vector<Edge> alignment_edges = {
      {0, 0}, {0, 1}, {1, 0}, {1, 1}, {2, 1}, {2, 2}};
  const std::vector<float> alignment_scores = {5.0f, 4.0f, 4.5f, 1.0f, 3.0f, 2.0f};
  const auto alignment_score = [&](unsigned i, unsigned k) {
    const auto it = std::find(alignment_edges.begin(), alignment_edges.end(), Edge{i, k});
    return alignment_scores[static_cast<size_t>(it - alignment_edges.begin())];
  };
  const float exact_alignment = exact_value(
      alignment_edges, alignment_scores, valid_alignment_subset);
  const float alignment_bound = RelaxedBounds::alignment(
      alignment_edges, 3, 3, alignment_score);
  const auto alignment_details = RelaxedBounds::alignment_bound(
      alignment_edges, 3, 3, alignment_score);
  std::vector<unsigned> relaxed_alignment;
  const float row_solution = RelaxedBounds::alignment_row_solution(
      alignment_edges, 3, 3, alignment_score, relaxed_alignment);
  if (!close(row_solution, alignment_details.row))
    return 1;
  float old_alignment_bound = 0.0f;
  for (const float value : alignment_scores)
    old_alignment_bound += std::max(0.0f, value);
  if (!check_bound(exact_alignment, alignment_bound, old_alignment_bound,
                   "alignment"))
    return 1;
  const float exact_monotone_alignment = exact_value(
      alignment_edges, alignment_scores, valid_monotone_subset);
  if (!check_alignment_components(exact_monotone_alignment,
                                  alignment_details))
    return 1;

  // The tight oracle must reject crossings, equal-row reuse, and equal-column
  // reuse while still returning a deterministic maximizing mapping.
  const std::vector<Edge> crossing_edges = {{0, 1}, {1, 0}};
  const std::vector<float> crossing_scores = {5.0f, 4.0f};
  const auto crossing_score = [&](unsigned i, unsigned k) {
    const auto it = std::find(crossing_edges.begin(), crossing_edges.end(),
                              Edge{i, k});
    return crossing_scores[static_cast<size_t>(it - crossing_edges.begin())];
  };
  const auto crossing_solution = RelaxedBounds::alignment_monotone_solution(
      crossing_edges, 2, 2, crossing_score);
  const auto crossing_capacity = RelaxedBounds::alignment_bound(
      crossing_edges, 2, 2, crossing_score);
  if (!close(crossing_solution.value, 5.0f) ||
      crossing_solution.selected[0] != 1 ||
      crossing_solution.selected[1] != std::numeric_limits<unsigned>::max() ||
      !check_monotone_mapping(crossing_solution.selected, 2,
                               crossing_score, crossing_solution.value) ||
      !check_alignment_components(5.0f, crossing_capacity))
    return 1;

  const std::vector<Edge> equal_row_edges = {{0, 0}, {0, 1}, {1, 2}};
  const std::vector<float> equal_row_scores = {3.0f, 4.0f, 2.0f};
  const auto equal_row_score = [&](unsigned i, unsigned k) {
    const auto it = std::find(equal_row_edges.begin(), equal_row_edges.end(),
                              Edge{i, k});
    return equal_row_scores[static_cast<size_t>(it - equal_row_edges.begin())];
  };
  const auto equal_row_solution = RelaxedBounds::alignment_monotone_solution(
      equal_row_edges, 2, 3, equal_row_score);
  const auto equal_row_capacity = RelaxedBounds::alignment_bound(
      equal_row_edges, 2, 3, equal_row_score);
  if (!close(equal_row_solution.value, 6.0f) ||
      equal_row_solution.selected[0] != 1 ||
      equal_row_solution.selected[1] != 2 ||
      !check_alignment_components(6.0f, equal_row_capacity))
    return 1;

  const std::vector<Edge> equal_column_edges = {
      {0, 0}, {1, 0}, {2, 1}, {2, 2}};
  const std::vector<float> equal_column_scores = {5.0f, 6.0f, 1.0f, 2.0f};
  const auto equal_column_score = [&](unsigned i, unsigned k) {
    const auto it = std::find(equal_column_edges.begin(),
                              equal_column_edges.end(), Edge{i, k});
    return equal_column_scores[static_cast<size_t>(it - equal_column_edges.begin())];
  };
  const auto equal_column_solution =
      RelaxedBounds::alignment_monotone_solution(
          equal_column_edges, 3, 3, equal_column_score);
  const auto equal_column_capacity = RelaxedBounds::alignment_bound(
      equal_column_edges, 3, 3, equal_column_score);
  if (!close(equal_column_solution.value, 8.0f) ||
      equal_column_solution.selected[0] != std::numeric_limits<unsigned>::max() ||
      equal_column_solution.selected[1] != 0 ||
      equal_column_solution.selected[2] != 2 ||
      !check_alignment_components(8.0f, equal_column_capacity))
    return 1;

  // A positive q_z edge with p_z=0 is represented here by its already-shifted
  // positive score.  The oracle must not silently drop it from the support.
  const std::vector<Edge> q_only_edges = {{0, 0}, {1, 1}, {0, 1}};
  const std::vector<float> q_only_scores = {2.0f, 2.5f, 3.0f};
  const auto q_only_score = [&](unsigned i, unsigned k) {
    const auto it = std::find(q_only_edges.begin(), q_only_edges.end(),
                              Edge{i, k});
    return q_only_scores[static_cast<size_t>(it - q_only_edges.begin())];
  };
  const auto q_only_solution = RelaxedBounds::alignment_monotone_solution(
      q_only_edges, 2, 2, q_only_score);
  const auto q_only_capacity = RelaxedBounds::alignment_bound(
      q_only_edges, 2, 2, q_only_score);
  if (!close(q_only_solution.value, 4.5f) ||
      q_only_solution.selected[0] != 0 ||
      q_only_solution.selected[1] != 1 ||
      !check_alignment_components(4.5f, q_only_capacity))
    return 1;

  const std::vector<Edge> zero_q_edges = {{0, 0}, {0, 1}, {1, 1}};
  const auto zero_q_score = [](unsigned, unsigned) { return 0.0f; };
  const auto zero_q_capacity = RelaxedBounds::alignment_bound(
      zero_q_edges, 2, 2, zero_q_score);
  if (!close(zero_q_capacity.row, 0.0f) ||
      !close(zero_q_capacity.column, 0.0f) ||
      !close(zero_q_capacity.best, 0.0f) ||
      !check_alignment_components(0.0f, zero_q_capacity))
    return 1;

  // Equal-valued optima must be independent of sparse-support insertion
  // order.  This also exercises the final-column tie breaker.
  const std::vector<Edge> tie_edges = {{0, 0}, {0, 1}, {1, 1}, {1, 2}};
  const auto tie_score = [](unsigned, unsigned) { return 2.0f; };
  const auto tie_solution = RelaxedBounds::alignment_monotone_solution(
      tie_edges, 2, 3, tie_score);
  std::vector<Edge> shuffled_tie_edges = tie_edges;
  std::reverse(shuffled_tie_edges.begin(), shuffled_tie_edges.end());
  const auto shuffled_tie_solution =
      RelaxedBounds::alignment_monotone_solution(
          shuffled_tie_edges, 2, 3, tie_score);
  if (!close(tie_solution.value, 4.0f) ||
      tie_solution.selected != shuffled_tie_solution.selected ||
      tie_solution.selected[0] != 0 || tie_solution.selected[1] != 1)
    return 1;

  // Out-of-range coordinates and duplicate support entries are harmless.
  const std::vector<Edge> malformed_edges = {
      {99, 0}, {0, 99}, {0, 0}, {0, 0}, {1, 1}};
  const auto malformed_score = [](unsigned i, unsigned k) {
    return (i == 0 && k == 0) ? 3.0f :
           (i == 1 && k == 1) ? 2.0f : -1.0f;
  };
  const auto malformed_solution = RelaxedBounds::alignment_monotone_solution(
      malformed_edges, 2, 2, malformed_score);
  if (!close(malformed_solution.value, 5.0f) ||
      !check_monotone_mapping(malformed_solution.selected, 2,
                               malformed_score, malformed_solution.value))
    return 1;

  // Differentially compare the exact positive monotone oracle with exhaustive
  // enumeration on tiny supports, including crossings and negative scores.
  std::mt19937 alignment_generator(20260919u);
  std::uniform_real_distribution<float> alignment_score_distribution(-3.0f,
                                                                       7.0f);
  for (unsigned trial = 0; trial < 200; ++trial) {
    constexpr unsigned rows = 4;
    constexpr unsigned columns = 4;
    std::vector<Edge> candidates;
    for (unsigned i = 0; i < rows; ++i)
      for (unsigned k = 0; k < columns; ++k)
        candidates.push_back({i, k});
    std::shuffle(candidates.begin(), candidates.end(), alignment_generator);
    candidates.resize(8);
    std::vector<float> scores(candidates.size());
    for (float& value : scores)
      value = alignment_score_distribution(alignment_generator);
    const auto random_score = [&](unsigned i, unsigned k) {
      const auto it = std::find(candidates.begin(), candidates.end(),
                                Edge{i, k});
      return scores[static_cast<size_t>(it - candidates.begin())];
    };
    const float exact = exact_value(candidates, scores,
                                    valid_monotone_subset);
    const auto solution = RelaxedBounds::alignment_monotone_solution(
        candidates, rows, columns, random_score);
    const auto loose_details = RelaxedBounds::alignment_bound(
        candidates, rows, columns, random_score);
    if (!close(solution.value, exact) ||
        !check_monotone_mapping(solution.selected, columns, random_score,
                                solution.value))
      return 1;
    if (loose_details.row + 1e-5f < exact ||
        loose_details.column + 1e-5f < exact ||
        loose_details.best + 1e-5f < exact)
      return 1;
  }

  // These supports deliberately contain competing edges, so the new bounds
  // must be strictly tighter than selecting every positive edge.
  if (!(structure_bound < old_structure_bound) ||
      !(alignment_bound < old_alignment_bound))
    return 1;

  // Exhaustively enumerate structures for many small deterministic random
  // supports.  This checks every new candidate bound independently; taking
  // their minimum is safe only if each candidate remains an upper bound.
  std::mt19937 generator(20260803u);
  std::uniform_real_distribution<float> score_distribution(-2.0f, 8.0f);
  for (unsigned trial = 0; trial < 100; ++trial) {
    constexpr unsigned length = 9;
    std::vector<Edge> candidates;
    for (unsigned i = 0; i < length; ++i)
      for (unsigned j = i+3; j < length; ++j)
        candidates.push_back({i, j});
    std::shuffle(candidates.begin(), candidates.end(), generator);
    candidates.resize(10);
    std::vector<float> scores(candidates.size());
    for (float& value : scores)
      value = score_distribution(generator);
    const auto random_score = [&](unsigned i, unsigned j) {
      const auto it = std::find(candidates.begin(), candidates.end(),
                                Edge{i, j});
      return scores[static_cast<size_t>(it - candidates.begin())];
    };
    const float exact = exact_value(candidates, scores,
                                    valid_structure_subset);
    const auto details = RelaxedBounds::structure_bound(
        candidates, length, random_score);
    if (!check_structure_components(exact, details))
      return 1;
  }

  return 0;
}
