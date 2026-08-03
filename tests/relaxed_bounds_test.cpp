#include "relaxed_bounds.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
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
