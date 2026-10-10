// MAFFT-compatible guide-tree algorithms reimplemented for DAFS.
// Reference: MAFFT 7.526 core/constants.c, Galign11.c, pairlocalalign.c,
// and mltaln9.c (setnearest / fixed_musclesupg_*_treeout_memsave).
// See THIRD_PARTY_NOTICES for the MAFFT copyright and BSD license.
#include "guide_tree.h"

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>

namespace guide_tree {
namespace {
constexpr unsigned int no_child = std::numeric_limits<unsigned int>::max();
constexpr double infinity = std::numeric_limits<double>::infinity();
constexpr double gap_open = -3599.0;  // MAFFT -f -2.00, nucleotide scaling/rounding
constexpr double gap_extend = -179.0; // MAFFT -g -0.10
const std::string alphabet = "agctuAGCTUnNbdhkmnrsvwyx-O";
using ScoreMatrix = std::array<std::array<int, 26>, 26>;

int rounded(double value)
{
  return static_cast<int>(value < 0 ? value - 0.5 : value + 0.5);
}

ScoreMatrix nucleotide_scores()
{
  // Uniform-frequency Kimura 2-parameter model (kappa=2), 200 PAM steps.
  std::array<std::array<double, 4>, 4> transition{}, pam{};
  for (size_t i = 0; i < 4; ++i) {
    pam[i][i] = 1;
    for (size_t j = 0; j < 4; ++j)
      transition[i][j] = i == j ? 0.99 : (i / 2 == j / 2 ? 0.005 : 0.0025);
  }
  for (size_t step = 0; step < 200; ++step) {
    std::array<std::array<double, 4>, 4> next{};
    for (size_t i = 0; i < 4; ++i)
      for (size_t j = 0; j < 4; ++j)
        for (size_t k = 0; k < 4; ++k)
          next[i][j] += pam[i][k] * transition[k][j];
    pam = next;
  }
  double average = 0;
  for (auto& row : pam)
    for (auto& value : row) {
      value = std::log10(value / 0.25) * 1000;
      average += value / 16;
    }
  for (auto& row : pam)
    for (auto& value : row) value -= average;
  double diagonal = 0;
  for (size_t i = 0; i < 4; ++i) diagonal += pam[i][i] / 4;

  ScoreMatrix scores{};
  for (size_t i = 0; i < 5; ++i)
    for (size_t j = 0; j < 5; ++j)
      scores[i][j] = scores[i + 5][j + 5] =
          rounded(pam[std::min(i, size_t{3})][std::min(j, size_t{3})] *
                  600 / diagonal - 59); // MAFFT -h 0.10 -> offset 59

  const std::array<std::string, 10> ambiguous_bases = {
      "ag", "ct", "gt", "ac", "gc", "at", "cgt", "agt", "act", "acg"};
  const std::string ambiguous = "rykmswbdhv";
  // Keep the reference rounding and update order, including ambiguous pairs.
  for (size_t i = 0; i < alphabet.size(); ++i) {
    for (size_t j = 0; j < ambiguous.size(); ++j) {
      double sum = 0;
      for (const char base : ambiguous_bases[j]) sum += scores[alphabet.find(base)][i];
      scores[i][alphabet.find(ambiguous[j])] = rounded(sum / ambiguous_bases[j].size());
    }
    for (const char symbol : ambiguous)
      scores[alphabet.find(symbol)][i] = scores[i][alphabet.find(symbol)];
  }
  for (size_t j = 0; j < ambiguous.size(); ++j) {
    double sum = 0;
    for (const char base : ambiguous_bases[j]) {
      const size_t index = alphabet.find(base);
      sum += scores[index][index];
    }
    const size_t index = alphabet.find(ambiguous[j]);
    scores[index][index] = rounded(sum / ambiguous_bases[j].size());
  }
  // n and x have zero scores in MAFFT's default (no --nwildcard).
  return scores;
}

std::vector<size_t> encode(const std::string& sequence)
{
  std::vector<size_t> result;
  for (const unsigned char symbol : sequence) {
    const auto index = alphabet.find(static_cast<char>(std::tolower(symbol)));
    if (index == std::string::npos)
      throw std::invalid_argument("Unsupported nucleotide in MAFFT guide tree");
    result.push_back(index);
  }
  return result;
}

double global_score(const std::vector<size_t>& first,
                    const std::vector<size_t>& second, const ScoreMatrix& scores)
{
  if (first.empty() || second.empty()) return 0;
  // Score-only affine DP, including MAFFT's terminal-gap convention: one
  // opening cost, no extension cost at either end. A sentinel match closes
  // the last gap. Only two rows and the vertical gap states are retained.
  const size_t n = first.size(), m = second.size();
  std::vector<double> current(m + 1), previous(m + 1), vertical(m + 1);
  const auto match = [&](size_t i, size_t j) -> double {
    return i < n && j < m ? scores[first[i]][second[j]] : 0;
  };
  const auto first_column = [&](size_t i) {
    return match(i, 0) + (i ? gap_open : 0);
  };
  for (size_t j = 0; j <= m; ++j)
    current[j] = match(0, j) + (j ? gap_open : 0);
  for (size_t j = 1; j <= m; ++j) vertical[j] = current[j - 1];

  double best = 0;
  for (size_t i = 1; i <= n; ++i) {
    current.swap(previous);
    previous[0] = first_column(i - 1);
    current[0] = first_column(i);
    double horizontal = previous[0];
    for (size_t j = 1; j <= m; ++j) {
      const double diagonal = previous[j - 1];
      best = std::max({diagonal, horizontal + gap_open, vertical[j] + gap_open});
      horizontal = std::max(horizontal, diagonal) + (i < n ? gap_extend : 0);
      vertical[j] = std::max(vertical[j], diagonal) + (j < m ? gap_extend : 0);
      current[j] = match(i, j) + best;
    }
  }
  return best;
}

std::vector<Node> empty_tree(size_t count)
{
  if (count == 0 || count > std::numeric_limits<unsigned int>::max() / 2)
    throw std::invalid_argument("Invalid sequence count for guide tree");
  return std::vector<Node>(2 * count - 1, {0.0f, {no_child, no_child}});
}
} // namespace

double mafft_global_distance(const std::string& first, const std::string& second)
{
  static const ScoreMatrix scores = nucleotide_scores();
  const auto a = encode(first), b = encode(second);
  double self_a = 0, self_b = 0;
  for (const auto symbol : a) self_a += scores[symbol][symbol];
  for (const auto symbol : b) self_b += scores[symbol][symbol];
  const double denominator = std::min(self_a, self_b);
  if (denominator == 0) return 2;
  const double score = global_score(a, b, scores);
  return denominator < score ? 0 : (1 - score / denominator) * 2;
}

std::vector<Node> mafft_linkage(std::vector<std::vector<double>> distances)
{
  const size_t n = distances.size();
  auto tree = empty_tree(n);
  for (size_t i = 0; i < n; ++i) {
    if (distances[i].size() != n) throw std::invalid_argument("Invalid MAFFT distance matrix");
    for (size_t j = i + 1; j < n; ++j)
      if (!std::isfinite(distances[i][j]) || distances[i][j] < 0)
        throw std::invalid_argument("Invalid MAFFT guide-tree distance");
  }
  std::vector<bool> active(n, true);
  std::vector<unsigned int> nodes(n), nearest(n, no_child);
  std::iota(nodes.begin(), nodes.end(), 0);
  std::vector<double> minimum(n, infinity);
  const auto distance = [&](size_t i, size_t j) -> double& {
    return distances[std::min(i, j)][std::max(i, j)];
  };
  const auto reset_nearest = [&](size_t i) {
    minimum[i] = infinity;
    nearest[i] = no_child;
    // MAFFT scans higher slots before lower slots and retains the first tie.
    for (size_t j = i + 1; j < n; ++j)
      if (active[j] && distance(i, j) < minimum[i]) {
        minimum[i] = distance(i, j);
        nearest[i] = j;
      }
    for (size_t j = 0; j < i; ++j)
      if (active[j] && distance(i, j) < minimum[i]) {
        minimum[i] = distance(i, j);
        nearest[i] = j;
      }
  };
  for (size_t i = 0; i < n; ++i) reset_nearest(i);
  for (size_t merge = 0; merge + 1 < n; ++merge) {
    size_t last = n - 1;
    while (!active[last]) --last;
    size_t left = 0;
    double best = infinity;
    for (size_t i = 0; i < last; ++i)
      if (active[i] && minimum[i] < best) {
        best = minimum[i];
        left = i;
      }
    size_t right = nearest[left];
    if (right == no_child) throw std::runtime_error("MAFFT guide tree has no active pair");
    if (right < left) std::swap(left, right);
    const unsigned int parent = n + merge;
    tree[parent].second = {nodes[left], nodes[right]};
    nodes[left] = parent;
    minimum[left] = infinity;
    for (size_t i = 0; i < n; ++i) if (active[i] && i != left && i != right) {
      const double a = distance(i, left), b = distance(i, right);
      const double updated = distance(i, left) = std::min(a, b) * 0.9 + (a + b) * 0.05;
      if (updated < minimum[i]) {
        minimum[i] = updated;
        nearest[i] = left;
      }
      if (updated < minimum[left]) {
        minimum[left] = updated;
        nearest[left] = i;
      }
      if (nearest[i] == right) nearest[i] = left;
    }
    active[right] = false;
    for (size_t i = 0; i < n; ++i)
      if (active[i] && nearest[i] == left && distance(i, left) > minimum[i])
        reset_nearest(i);
  }
  return tree;
}

std::vector<Node> mafft(const std::vector<Fasta>& sequences)
{
  const size_t n = sequences.size();
  if (n < 2) return empty_tree(n);
  std::vector<std::vector<double>> distances(n, std::vector<double>(n));
  for (size_t i = 0; i < n; ++i)
    for (size_t j = i + 1; j < n; ++j)
      distances[i][j] = mafft_global_distance(sequences[i].seq(), sequences[j].seq());
  return mafft_linkage(std::move(distances));
}

std::vector<Node> raf(std::vector<std::vector<float>> distances)
{
  const size_t n = distances.size();
  auto tree = empty_tree(n);
  std::vector<unsigned int> nodes(n);
  std::iota(nodes.begin(), nodes.end(), 0);
  for (size_t i = 0; i < n; ++i) {
    if (distances[i].size() != n) throw std::invalid_argument("Invalid RAF distance matrix");
    distances[i][i] = 1e10f;
  }
  for (size_t merge = 0; merge + 1 < n; ++merge) {
    size_t left = 0, right = 0;
    for (size_t i = 0; i < n; ++i)
      for (size_t j = i + 1; j < n; ++j)
        if (distances[i][j] < distances[left][right]) { left = i; right = j; }
    if (left == right || nodes[left] == no_child || nodes[right] == no_child)
      throw std::runtime_error("RAF guide tree has no active pair");
    const unsigned int parent = n + merge;
    tree[parent] = {1.0f - distances[left][right], {nodes[left], nodes[right]}};
    nodes[left] = parent;
    nodes[right] = no_child;
    distances[left][right] = 1e10f;
    for (size_t i = 0; i < n; ++i) if (i != left && i != right) {
      distances[i][left] = distances[left][i] =
          0.9f * std::min(distances[i][left], distances[i][right]) +
          0.1f * (distances[i][left] + distances[i][right]) / 2.0f;
      distances[i][right] = distances[right][i] = 1e10f;
    }
  }
  return tree;
}
} // namespace guide_tree
