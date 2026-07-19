#include "relaxed_bounds.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <utility>
#include <vector>

namespace
{
using Edge = std::pair<unsigned, unsigned>;

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
  const float structure_bound = RelaxedBounds::structure(
      structure_edges, 7, structure_score);
  float old_structure_bound = 0.0f;
  for (size_t e = 0; e < structure_edges.size(); ++e)
    if (structure_edges[e].second > structure_edges[e].first + 2)
      old_structure_bound += std::max(0.0f, structure_scores[e]);
  if (!check_bound(exact_structure, structure_bound, old_structure_bound,
                   "structure"))
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

  return 0;
}
