/*
 * Copyright (C) 2012 Kengo Sato
 *
 * This file is part of DAFS.
 *
 * DAFS is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * DAFS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with DAFS.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "ipknot.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

namespace {
ipknot::DecoderOptions decoder_options(const VF& thresholds,
                                       ipknot::Backend backend,
                                       int threads, bool allow_isolated, uint beam = 100)
{
  ipknot::DecoderOptions options;
  options.thresholds = thresholds;
  // DAFS sums the level objectives; IPknot otherwise averages them.
  options.weights.assign(thresholds.size(), 1.0f);
  options.backend = backend;
  options.threads = threads;
  options.no_lonely_pairs = !allow_isolated;
  if (beam > static_cast<uint>(std::numeric_limits<int>::max()))
    throw std::invalid_argument("IPknot beam is too large");
  options.dd.beam = std::max(1u, beam);
  return options;
}

void check_length(uint length)
{
  if (length > static_cast<uint>(std::numeric_limits<int>::max()))
    throw std::invalid_argument("IPknot sequence is too long");
}

void check_dense(const VVF& matrix, uint length)
{
  if (matrix.size() != length)
    throw std::invalid_argument("IPknot matrix dimensions differ");
  for (const auto& row : matrix)
    if (row.size() != length)
      throw std::invalid_argument("IPknot requires square matrices");
}

void check_sparse(const SparseFloatMatrix& matrix, uint length)
{
  if (matrix.rows() != length || matrix.columns() != length)
    throw std::invalid_argument("IPknot matrix dimensions differ");
}
} // namespace

IPknot::IPknot(const VF& th, int n_th, bool allow_isolated)
  : IPknot(decoder_options(th, ipknot::Backend::ILP, n_th, allow_isolated))
{
}

IPknot::IPknot(ipknot::DecoderOptions options)
  : options_(std::move(options))
{
  // Validate the original thresholds before any score normalization.
  const ipknot::Decoder validator(options_);
  if (options_.thresholds.size() > n_support_brackets)
    throw std::invalid_argument("Too many IPknot levels for bracket output");
}

LinearIPknot::LinearIPknot(const VF& th, uint beam_size, bool allow_isolated)
  : IPknot(decoder_options(th, ipknot::Backend::DD, 1, allow_isolated, beam_size))
{
}

float IPknot::decode(float w, const VVF& p, const VVF& q, VU& ss)
{
  const uint length = p.size();
  check_length(length);
  check_dense(p, length);
  check_dense(q, length);
  std::vector<ScoredPair> pairs;
  for (uint i = 0; i < length; ++i)
    for (uint j = i + 1; j < length; ++j)
      pairs.push_back({i, j, double(w) * p[i][j] - q[i][j]});
  return decode_scores(w, length, pairs, ss);
}

float IPknot::decode(float w, const VVF& p,
                     const SparseFloatMatrix& q, VU& ss)
{
  const uint length = p.size();
  check_length(length);
  check_dense(p, length);
  check_sparse(q, length);
  std::vector<ScoredPair> pairs;
  for (uint i = 0; i < length; ++i)
    for (uint j = i + 1; j < length; ++j)
      pairs.push_back({i, j, double(w) * p[i][j] - q.get(i, j)});
  return decode_scores(w, length, pairs, ss);
}

float IPknot::decode(float w, const SparseFloatMatrix& p,
                     const VVF& q, VU& ss)
{
  const uint length = p.rows();
  check_length(length);
  check_sparse(p, length);
  check_dense(q, length);
  std::vector<ScoredPair> pairs;
  for (uint i = 0; i < length; ++i)
    for (uint j = i + 1; j < length; ++j)
      pairs.push_back({i, j, double(w) * p.get(i, j) - q[i][j]});
  return decode_scores(w, length, pairs, ss);
}

float IPknot::decode(float w, const SparseFloatMatrix& p,
                     const SparseFloatMatrix& q, VU& ss)
{
  const uint length = p.rows();
  check_length(length);
  check_sparse(p, length);
  check_sparse(q, length);
  std::vector<ScoredPair> pairs;
  pairs.reserve(p.nonzeros() + q.nonzeros());
  // Scan the support union, including q-only pairs made profitable by DD.
  p.for_each_nonzero([&](uint i, uint j, float probability) {
    if (i < j)
      pairs.push_back({i, j, double(w) * probability - q.get(i, j)});
  });
  q.for_each_nonzero([&](uint i, uint j, float multiplier) {
    if (i < j && p.get(i, j) == 0.0f)
      pairs.push_back({i, j, -double(multiplier)});
  });
  return decode_scores(w, length, pairs, ss);
}

float IPknot::decode(const VVF& p, VU& ss, std::string& str)
{
  SparseFloatMatrix empty_q;
  empty_q.assign(p.size(), p.size());
  const float score = decode(1.0f, p, empty_q, ss);
  make_brackets(ss, str);
  return score;
}

float IPknot::decode(const SparseFloatMatrix& p, VU& ss, std::string& str)
{
  SparseFloatMatrix empty_q;
  empty_q.assign(p.rows(), p.columns());
  const float score = decode(1.0f, p, empty_q, ss);
  make_brackets(ss, str);
  return score;
}

float IPknot::decode_scores(float w, uint length,
                            const std::vector<ScoredPair>& pairs, VU& ss)
{
  if (!std::isfinite(w))
    throw std::invalid_argument("IPknot weight must be finite");
  auto options = options_;
  double minimum_threshold = std::numeric_limits<double>::infinity();
  double maximum_threshold = -std::numeric_limits<double>::infinity();
  for (float threshold : options.thresholds) {
    minimum_threshold = std::min(minimum_threshold, double(w) * threshold);
    maximum_threshold = std::max(maximum_threshold, double(w) * threshold);
  }
  const double offset = std::min(0.0, minimum_threshold);
  double maximum_score = std::max(1.0, maximum_threshold);
  for (const auto& pair : pairs) {
    if (!std::isfinite(pair.score))
      throw std::invalid_argument("IPknot pair scores must be finite");
    maximum_score = std::max(maximum_score, pair.score);
  }
  const double scale = maximum_score - offset;
  if (scale > std::numeric_limits<float>::max())
    throw std::invalid_argument("IPknot pair scores are too large");

  // The public API accepts probabilities in [0,1]. A shared affine transform
  // preserves every level's DAFS objective, including scores above one and
  // negative multipliers: (p' - th'[level]) * scale = w*(p-th[level]) - q.
  for (float& threshold : options.thresholds)
    threshold = (double(w) * threshold - offset) / scale;
  options.weights.assign(options.thresholds.size(), static_cast<float>(scale));
  std::vector<ipknot::PairProbability> probabilities;
  probabilities.reserve(pairs.size());
  for (const auto& pair : pairs)
    if (pair.score > minimum_threshold)
      probabilities.push_back({static_cast<int>(pair.left),
                               static_cast<int>(pair.right),
                               static_cast<float>((pair.score - offset) / scale)});

  const auto result = ipknot::Decoder(options).decode(length, probabilities);
  ss.assign(length, -1u);
  levels_.assign(length, -1u);
  for (uint i = 0; i < length; ++i) {
    if (result.bpseq[i] > static_cast<int>(i))
      ss[i] = result.bpseq[i]; // DAFS stores only the left endpoint's partner.
    if (result.levels[i] >= 0)
      levels_[i] = result.levels[i];
  }
  last_structure_ = ss;
  // Evaluate in the original units, without affine-transform rounding.
  double score = 0;
  for (const auto& pair : pairs)
    if (ss[pair.left] == pair.right)
      score += pair.score - double(w) * options_.thresholds[levels_[pair.left]];
  return static_cast<float>(score);
}

VU IPknot::decompose_levels(const VU& ss)
{
  // Projected structures have no cached levels. Greedy coloring of active
  // nested arcs takes O(length * bracket_count), with a fixed alphabet.
  VU levels(ss.size(), -1u);
  std::vector<std::vector<std::pair<uint, uint>>> active(n_support_brackets);
  for (uint left = 0; left < ss.size(); ++left) {
    const uint right = ss[left];
    if (right == -1u || right <= left || right >= ss.size())
      continue;
    for (uint level = 0; level < n_support_brackets; ++level) {
      auto& stack = active[level];
      while (!stack.empty() && stack.back().second <= left)
        stack.pop_back();
      if (stack.empty() || stack.back().second >= right) {
        levels[left] = levels[right] = level;
        stack.emplace_back(left, right);
        break;
      }
    }
  }
  return levels;
}

void IPknot::make_brackets(const VU& ss, std::string& str) const
{
  const VU levels = last_structure_ == ss && levels_.size() == ss.size()
      ? levels_ : decompose_levels(ss);
  str.assign(ss.size(), '.');
  for (uint i = 0; i < ss.size(); ++i) {
    const uint j = ss[i];
    if (j != -1u && i < j && j < ss.size() && levels[i] < n_support_brackets) {
      str[i] = left_brackets[levels[i]];
      str[j] = right_brackets[levels[i]];
    }
  }
}
