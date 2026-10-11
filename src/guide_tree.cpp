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

double band_global_score(const std::vector<size_t>& first,
                         const std::vector<size_t>& second, const ScoreMatrix& scores,
                         size_t width, BeamStats* stats)
{
  // The exact MAFFT recurrence restricted to |i-j| <= width. Its innermost
  // loop is the same row-wise kernel as global_score; it has no beam merge,
  // candidate entries, or top-B selection. All unvisited states stay invalid.
  const size_t n = first.size(), m = second.size();
  if (stats) stats->band_width = std::max(stats->band_width, width);
  std::vector<double> current(m + 1, -infinity), previous(m + 1, -infinity),
                      vertical(m + 1, -infinity);
  const auto match = [&](size_t i, size_t j) -> double {
    return i < n && j < m ? scores[first[i]][second[j]] : 0;
  };
  for (size_t j = 0; j <= std::min(m, width); ++j)
    current[j] = match(0, j) + (j ? gap_open : 0);
  for (size_t j = 1; j <= std::min(m, width + 1); ++j) vertical[j] = current[j - 1];
  double best = -infinity;
  for (size_t i = 1; i <= n; ++i) {
    current.swap(previous);
    previous[0] = i - 1 <= width ? match(i - 1, 0) + (i > 1 ? gap_open : 0) : -infinity;
    current[0] = i <= width ? match(i, 0) + gap_open : -infinity;
    const size_t start = i > width ? i - width : 1;
    const size_t end = std::min(m, i + width);
    if (end < m) current[end + 1] = -infinity;
    double horizontal = start == 1 ? previous[0] : -infinity;
    for (size_t j = start; j <= end; ++j) {
      const double diagonal = previous[j - 1];
      best = std::max({diagonal, horizontal + gap_open, vertical[j] + gap_open});
      horizontal = std::max(horizontal, diagonal) + (i < n ? gap_extend : 0);
      vertical[j] = std::max(vertical[j], diagonal) + (j < m ? gap_extend : 0);
      current[j] = match(i, j) + best;
    }
    if (stats && start <= end) stats->band_cells += end - start + 1;
  }
  return best;
}

std::pair<bool, double> certified_band_score(const std::vector<size_t>& first,
                                            const std::vector<size_t>& second,
                                            const ScoreMatrix& scores, size_t beam,
                                            BeamStats* stats)
{
  const size_t n = first.size(), m = second.size(), length = std::max(n, m);
  // A band cell uses the inexpensive row kernel, whereas a beam entry needs
  // merging/selection and three separate states. Bound the exact accelerator
  // by a fixed multiple of the user's beam; its budget never grows with L.
  constexpr size_t band_factor = 16;
  const size_t limit = beam > length / band_factor ? length : band_factor * beam;
  const size_t difference = n > m ? n - m : m - n;
  if (difference > limit) return {false, 0};
  size_t width = std::min(limit, std::max(beam, difference));
  double best = band_global_score(first, second, scores, width, stats);
  static const int maximum_score = [&] {
    int result = 0;
    for (const auto& row : scores)
      result = std::max(result, *std::max_element(row.begin(), row.end()));
    return result;
  }();
  // If a path leaves this band, at least width+1 bases in one sequence are
  // gaps. It can have at most max(n,m)-width-1 matches, each worth at most
  // maximum_score; all gap scores are nonpositive. This is an upper bound
  // on EVERY excluded path, including terminal gaps and ambiguous bases.
  if (width < length) {
    const double required = double(length - 1) - std::floor(best / maximum_score);
    if (required > double(limit)) return {false, 0};
    if (required > double(width)) {
      width = static_cast<size_t>(required);
      best = band_global_score(first, second, scores, width, stats);
    }
    if (width < length && best < maximum_score * double(length - width - 1))
      return {false, 0};
  }
  if (stats) stats->certified_exact = true;
  return {true, best};
}

struct Entry {
  size_t i;
  double score;
};
struct Frontier {
  // Coordinate-indexed scores make both recurrence evaluation and duplicate
  // detection independent of the ordering of the retained entries.
  std::array<std::vector<Entry>, 3> entries;
  std::array<std::vector<double>, 3> scores;
  size_t lower = std::numeric_limits<size_t>::max(), upper = 0;
};

double beam_global_score(const std::vector<size_t>& first,
                         const std::vector<size_t>& second, const ScoreMatrix& scores,
                         size_t beam, BeamStats* stats)
{
  // Consumed-prefix affine Viterbi DP. Keep three antidiagonals and at most
  // beam entries per state. The lookup arrays are linear in the shorter
  // sequence; only retained coordinates are cleared when a ring slot is reused.
  const size_t n = first.size(), m = second.size();
  std::array<Frontier, 3> ring;
  for (auto& frontier : ring)
    for (auto& values : frontier.scores) values.assign(n + 1, -infinity);
  ring[0].entries[0].push_back({0, 0});
  ring[0].scores[0][0] = 0;
  ring[0].lower = ring[0].upper = 0;
  const auto better = [](const Entry& a, const Entry& b) {
    return a.score != b.score ? a.score > b.score : a.i < b.i;
  };
  for (size_t step = 1; step <= n + m; ++step) {
    auto& current = ring[step % 3];
    const auto& previous = ring[(step + 2) % 3];
    const auto& older = ring[(step + 1) % 3];
    for (size_t state = 0; state < 3; ++state) {
      for (const auto& entry : current.entries[state])
        current.scores[state][entry.i] = -infinity;
      current.entries[state].clear();
    }
    current.lower = std::numeric_limits<size_t>::max();
    current.upper = 0;
    const size_t minimum_i = step > m ? step - m : 0;
    const size_t maximum_i = std::min(step, n);
    // These transitions and boundary rules are identical to the sparse beam
    // recurrence. In particular, a gap cannot switch into the opposite gap.
    const auto evaluate = [&](size_t state, size_t i) {
      if (i < minimum_i || i > maximum_i) return;
      if (current.scores[state][i] != -infinity) return;
      const size_t j = step - i;
      double best;
      if (state == 0) {
        if (step < 2 || i == 0 || j == 0) return;
        best = std::max({older.scores[0][i - 1], older.scores[1][i - 1],
                         older.scores[2][i - 1]}) + scores[first[i - 1]][second[j - 1]];
      } else if (state == 1) {
        if (i == 0 || (i == n && j < m)) return;
        best = std::max(previous.scores[0][i - 1] + gap_open,
                        previous.scores[1][i - 1]);
        if (j > 0 && j < m) best += gap_extend;
      } else {
        if (j == 0 || (j == m && i < n)) return;
        best = std::max(previous.scores[0][i] + gap_open, previous.scores[2][i]);
        if (i > 0 && i < n) best += gap_extend;
      }
      if (best == -infinity) return;
      current.scores[state][i] = best;
      current.entries[state].push_back({i, best});
    };
    size_t lower = previous.lower, upper = previous.upper + 1;
    if (older.lower != std::numeric_limits<size_t>::max()) {
      lower = std::min(lower, older.lower + 1);
      upper = std::max(upper, older.upper + 1);
    }
    lower = std::max(lower, minimum_i);
    upper = std::min(upper, maximum_i);
    // Typical beams occupy a short interval. Evaluate its three states in a
    // single coordinate pass. A hard O(beam) span limit prevents quadratic
    // scanning when retained coordinates are widely separated.
    if (lower <= upper && (upper - lower) / beam < 3) {
      for (size_t i = lower; i <= upper; ++i) {
        evaluate(0, i);
        evaluate(1, i);
        evaluate(2, i);
      }
    } else {
      for (size_t state = 0; state < 3; ++state) {
        const auto& source = state == 0 ? older : previous;
        const size_t shift = state == 2 ? 0 : 1;
        for (size_t s = 0; s < 3; ++s) {
          if ((state == 1 && s == 2) || (state == 2 && s == 1)) continue;
          for (const auto& entry : source.entries[s]) evaluate(state, entry.i + shift);
        }
      }
    }
    for (size_t state = 0; state < 3; ++state) {
      auto& output = current.entries[state];
      if (stats) stats->peak_candidates = std::max(stats->peak_candidates, output.size());
      if (output.size() > beam) {
        // Total ordering by (score, coordinate) preserves exactly the old
        // retained set without copying, sorting, or coordinate-order erasure.
        std::nth_element(output.begin(), output.begin() + beam - 1, output.end(), better);
        for (size_t index = beam; index < output.size(); ++index)
          current.scores[state][output[index].i] = -infinity;
        output.resize(beam);
      }
      for (const auto& entry : output) {
        current.lower = std::min(current.lower, entry.i);
        current.upper = std::max(current.upper, entry.i);
      }
      if (stats) {
        stats->retained_states += output.size();
        stats->peak_retained = std::max(stats->peak_retained, output.size());
      }
    }
  }
  double best = -infinity;
  for (const auto& entries : ring[(n + m) % 3].entries)
    for (const auto& entry : entries) best = std::max(best, entry.score);
  if (!std::isfinite(best)) throw std::runtime_error("MAFFT beam search has no complete path");
  return best;
}

std::vector<Node> empty_tree(size_t count)
{
  if (count == 0 || count > std::numeric_limits<unsigned int>::max() / 2)
    throw std::invalid_argument("Invalid sequence count for guide tree");
  return std::vector<Node>(2 * count - 1, {0.0f, {no_child, no_child}});
}
} // namespace

double mafft_global_distance(const std::string& first, const std::string& second,
                             size_t beam, BeamStats* stats)
{
  if (stats) *stats = BeamStats{};
  static const ScoreMatrix scores = nucleotide_scores();
  auto a = encode(first), b = encode(second);
  double self_a = 0, self_b = 0;
  for (const auto symbol : a) self_a += scores[symbol][symbol];
  for (const auto symbol : b) self_b += scores[symbol][symbol];
  const double denominator = std::min(self_a, self_b);
  if (denominator == 0) return 2;
  // Canonical orientation makes pruning ties independent of which sequence
  // happened to be passed first. Exact DP needs no such tie convention.
  if (beam && (a.size() > b.size() || (a.size() == b.size() && b < a))) a.swap(b);
  const auto certified = beam ? certified_band_score(a, b, scores, beam, stats)
                              : std::pair<bool, double>{false, 0};
  const double score = certified.first ? certified.second
      : beam ? beam_global_score(a, b, scores, beam, stats) : global_score(a, b, scores);
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

std::vector<Node> mafft(const std::vector<Fasta>& sequences, size_t beam)
{
  const size_t n = sequences.size();
  if (n < 2) return empty_tree(n);
  std::vector<std::vector<double>> distances(n, std::vector<double>(n));
  for (size_t i = 0; i < n; ++i)
    for (size_t j = i + 1; j < n; ++j)
      distances[i][j] = mafft_global_distance(sequences[i].seq(), sequences[j].seq(), beam);
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
