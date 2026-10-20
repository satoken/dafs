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

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "nussinov.h"
#include "relaxed_bounds.h"
#include <ipknot/dual_decomposition.h>
#include <array>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <limits>
#include <numeric>
#include <stack>
#include <tuple>
#include <exception>
#include <stdexcept>
#ifdef _OPENMP
#include <omp.h>
#endif

static
std::string
make_brackets(const VU& ss);

namespace {
void validate_dp_threads(int threads)
{
  if (threads < 0)
    throw std::invalid_argument("DP thread count must be nonnegative");
#ifndef _OPENMP
  if (threads > 1)
    throw std::invalid_argument("Multiple DP threads require OpenMP support");
#endif
}

int dp_worker_limit(int threads, uint length)
{
#ifdef _OPENMP
  return std::max(1, static_cast<int>(std::min<size_t>(
      length, threads ? threads : omp_get_max_threads())));
#else
  return 1;
#endif
}

// Missing posterior/multiplier entries are not pair candidates. Retained
// candidates with non-positive scores still matter: they can support a
// profitable neighboring pair in a helix.
float nussinov_pair_score(float w, float probability, float multiplier,
                         float threshold, bool allow_isolated)
{
  if (!allow_isolated && probability == 0.0f && multiplier == 0.0f)
    return -std::numeric_limits<float>::infinity();
  return w * (probability - threshold) - multiplier;
}

// F is a complete structure with no isolated pairs. P has paired endpoints
// that may need an outer stacked pair. Thus P=s+max(F_inside,P_inside),
// while a complete paired block is Q=s+P_inside. F concatenates Q blocks.
template <typename PairScore>
float decode_stacked_nussinov(uint length, PairScore pair_score, VU& ss,
                              int threads, int& used_workers, size_t& parallel_steps)
{
  used_workers = 1;
  parallel_steps = 0;
  ss.assign(length, -1u);
  if (length == 0)
    return 0.0f;
  const double missing = -std::numeric_limits<double>::infinity();
  std::vector<std::vector<double>> complete(length, std::vector<double>(length));
  std::vector<std::vector<double>> paired(length, std::vector<double>(length, missing));
  VVU trace(length, VU(length, -1u));
  std::vector<std::vector<bool>> extend(length, std::vector<bool>(length));
  // Each interval length depends only on shorter intervals. Closed blocks
  // belong to their right endpoint, so workers never mutate the same list.
  std::vector<std::vector<std::pair<uint, double>>> closed(length);
  for (uint right = 0; right < length; ++right) closed[right].reserve(right);
  const int workers = dp_worker_limit(threads, length);
#ifdef _OPENMP
#pragma omp parallel num_threads(workers) if(workers > 1)
#endif
  {
#ifdef _OPENMP
#pragma omp single
    used_workers = omp_get_num_threads();
#endif
    for (uint span = 1; span < length; ++span) {
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
      for (uint left = 0; left < length-span; ++left) {
        const uint right = left + span;
        auto& blocks = closed[right];
        complete[left][right] = complete[left][right-1];
        if (span >= 3) {
          const double local = pair_score(left, right);
          if (std::isfinite(local)) {
            const double inner_pair = paired[left+1][right-1];
            extend[left][right] = inner_pair > complete[left+1][right-1];
            paired[left][right] = local + std::max(
                inner_pair, complete[left+1][right-1]);
            const double block = local + inner_pair;
            if (block > 0.0) blocks.emplace_back(left, block);
          }
        }
        for (const auto [start, block] : blocks) {
          const double prefix = start > left ? complete[left][start-1] : 0.0;
          if (prefix + block > complete[left][right]) {
            complete[left][right] = prefix + block;
            trace[left][right] = start;
          }
        }
      }
    }
  }
  if (used_workers > 1) parallel_steps = length - 1;
  std::vector<std::tuple<uint, uint, bool>> pending{{0, length-1, false}};
  while (!pending.empty()) {
    const auto [left, right, is_paired] = pending.back();
    pending.pop_back();
    if (is_paired) {
      ss[left] = right;
      pending.emplace_back(left+1, right-1, extend[left][right]);
    } else if (trace[left][right] != -1u) {
      const uint start = trace[left][right];
      ss[start] = right;
      if (left < start)
        pending.emplace_back(left, start-1, false);
      pending.emplace_back(start+1, right-1, true);
    } else if (left < right) {
      pending.emplace_back(left, right-1, false);
    }
  }
  return static_cast<float>(complete[0][length-1]);
}
} // namespace

Nussinov::Nussinov(float th, bool allow_isolated, int dp_threads)
    : th_(th), allow_isolated_(allow_isolated), dp_threads_(dp_threads)
{
  validate_dp_threads(dp_threads);
}

float
Nussinov::
decode(float w, const VVF& p, const VVF& q, VU& ss)
{
  if (!allow_isolated_)
    return decode_stacked_nussinov(p.size(), [&](uint i, uint j) {
      return nussinov_pair_score(w, p[i][j], q[i][j], th_, false);
    }, ss, dp_threads_, used_workers_, parallel_steps_);
  uint L=p.size();
  used_workers_ = 1;
  parallel_steps_ = 0;
  if (L == 0) { ss.clear(); return 0.0f; }
  assert(p[0].size()==L);

  // calculate scoring matrices for the current step
  VVF sm(L, VF(L, 0.0));
  for (uint i=0; i!=L-1; ++i)
    for (uint j=i+1; j!=L; ++j)
      sm[i][j] = w*(p[i][j]-th_)-q[i][j];

  VVF dp(L, VF(L, 0.0));
  VVU tr(L, VU(L, 0));
  const int workers = dp_worker_limit(dp_threads_, L);
#ifdef _OPENMP
#pragma omp parallel num_threads(workers) if(workers > 1)
#endif
  {
#ifdef _OPENMP
#pragma omp single
    used_workers_ = omp_get_num_threads();
#endif
    for (uint l=1; l<L; ++l)
    {
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
      for (uint i=0; i<L-l; ++i)
      {
        uint j=i+l;
        float v=0.0;
        int t=0;
        if (i+1<j)
        {
          v=dp[i+1][j];
          t=1;
        }
        if (i<j-1 && v<dp[i][j-1])
        {
          v=dp[i][j-1];
          t=2;
        }
        if (i+1<j-1 && v<dp[i+1][j-1]+sm[i][j])
        {
          v=dp[i+1][j-1]+sm[i][j];
          t=3;
        }
        for (uint k=i+1; k<j; ++k)
        {
          if (v<dp[i][k]+dp[k+1][j])
          {
            v=dp[i][k]+dp[k+1][j];
            t=k-i+3;
          }
        }
        dp[i][j]=v;
        tr[i][j]=t;
      }
    }
  }
  if (used_workers_ > 1) parallel_steps_ = L - 1;

  // trace back
  ss.resize(L);
  std::fill(ss.begin(), ss.end(), -1u);
  std::stack<std::pair<uint,uint> > st;
  st.push(std::make_pair(0, L-1));
  while (!st.empty())
  {
    const std::pair<uint,uint> p=st.top(); st.pop();
    const int i=p.first, j=p.second;
    switch (tr[i][j])
    {
      case 0:
        break;
      case 1:
        st.push(std::make_pair(i+1, j));
        break;
      case 2:
        st.push(std::make_pair(i, j-1));
        break;
      case 3:
        ss[i]=j;
        st.push(std::make_pair(i+1, j-1));
        break;
      default:
        const int k=i+tr[i][j]-3;
        st.push(std::make_pair(i, k));
        st.push(std::make_pair(k+1, j));
        break;
    }
  }
  return dp[0][L-1];
}

float
Nussinov::
decode(const VVF& p, VU& ss, std::string& str)
{
  if (!allow_isolated_) {
    const float score = decode_stacked_nussinov(p.size(), [&](uint i, uint j) {
      return nussinov_pair_score(1.0f, p[i][j], 0.0f, th_, false);
    }, ss, dp_threads_, used_workers_, parallel_steps_);
    make_brackets(ss, str);
    return score;
  }
  uint L=p.size();
  used_workers_ = 1;
  parallel_steps_ = 0;
  if (L == 0) { ss.clear(); str.clear(); return 0.0f; }
  assert(p[0].size()==L);

  // calculate scoring matrices for the current step
  VVF sm(L, VF(L, 0.0));
  for (uint i=0; i!=L-1; ++i)
    for (uint j=i+1; j!=L; ++j)
      sm[i][j] = p[i][j]-th_;

  VVF dp(L, VF(L, 0.0));
  VVU tr(L, VU(L, 0));
  const int workers = dp_worker_limit(dp_threads_, L);
#ifdef _OPENMP
#pragma omp parallel num_threads(workers) if(workers > 1)
#endif
  {
#ifdef _OPENMP
#pragma omp single
    used_workers_ = omp_get_num_threads();
#endif
    for (uint l=1; l<L; ++l)
    {
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
      for (uint i=0; i<L-l; ++i)
      {
        uint j=i+l;
        float v=0.0;
        int t=0;
        if (i+1<j)
        {
          v=dp[i+1][j];
          t=1;
        }
        if (i<j-1 && v<dp[i][j-1])
        {
          v=dp[i][j-1];
          t=2;
        }
        if (i+1<j-1 && v<dp[i+1][j-1]+sm[i][j])
        {
          v=dp[i+1][j-1]+sm[i][j];
          t=3;
        }
        for (uint k=i+1; k<j; ++k)
        {
          if (v<dp[i][k]+dp[k+1][j])
          {
            v=dp[i][k]+dp[k+1][j];
            t=k-i+3;
          }
        }
        dp[i][j]=v;
        tr[i][j]=t;
      }
    }
  }
  if (used_workers_ > 1) parallel_steps_ = L - 1;

  // trace back
  ss.resize(L);
  std::fill(ss.begin(), ss.end(), -1u);
  std::stack<std::pair<uint,uint> > st;
  st.push(std::make_pair(0, L-1));
  while (!st.empty())
  {
    const std::pair<uint,uint> p=st.top(); st.pop();
    const int i=p.first, j=p.second;
    switch (tr[i][j])
    {
      case 0:
        break;
      case 1:
        st.push(std::make_pair(i+1, j));
        break;
      case 2:
        st.push(std::make_pair(i, j-1));
        break;
      case 3:
        ss[i]=j;
        st.push(std::make_pair(i+1, j-1));
        break;
      default:
        const int k=i+tr[i][j]-3;
        st.push(std::make_pair(i, k));
        st.push(std::make_pair(k+1, j));
        break;
    }
  }

  make_brackets(ss, str);
  return dp[0][L-1];
}

float
Nussinov::
decode(const SparseFloatMatrix& p, const SparseFloatMatrix& pair_bonus,
       VU& ss, std::string& str)
{
  assert(p.rows() == pair_bonus.rows() &&
         p.columns() == pair_bonus.columns());
  if (!allow_isolated_) {
    const float score = decode_stacked_nussinov(p.rows(), [&](uint i, uint j) {
      return nussinov_pair_score(1.0f, p.get(i, j), -pair_bonus.get(i, j), th_, false);
    }, ss, dp_threads_, used_workers_, parallel_steps_);
    make_brackets(ss, str);
    return score;
  }
  VVF adjusted = p.dense();
  for (uint i = 0; i < pair_bonus.rows(); ++i)
    for (const auto [j, value] : pair_bonus.ordered_row(i))
      adjusted[i][j] += value;
  return decode(adjusted, ss, str);
}

void
Nussinov::
make_brackets(const VU& ss, std::string& str) const
{
  str=::make_brackets(ss);
}

namespace {
void
radix_sort_linear_nussinov_starts(std::vector<uint>& values,
                                  std::vector<uint>& scratch)
{
  if (values.size() < 2)
    return;

  scratch.resize(values.size());
  constexpr uint radix_bits = 8;
  constexpr uint radix_size = 1u << radix_bits;
  std::array<size_t, radix_size> counts;
  for (uint shift = 0; shift < sizeof(uint) * 8; shift += radix_bits) {
    counts.fill(0);
    for (const uint value : values)
      ++counts[(value >> shift) & (radix_size - 1)];

    size_t offset = 0;
    for (size_t& count : counts) {
      const size_t next = offset + count;
      count = offset;
      offset = next;
    }
    for (const uint value : values)
      scratch[counts[(value >> shift) & (radix_size - 1)]++] = value;
    values.swap(scratch);
  }
}

// A caller-supplied sparse support can contain span-two pairs.  The linear
// recurrence can decode such a pair because it leaves one position inside;
// adjacent pairs cannot be decoded.  Use one predicate in both the decoder
// and its certificate so they always optimize and bound exactly the same
// support without changing the historical DAFS search space.
bool
linear_nussinov_pair_is_usable(uint left, uint right, uint length)
{
  return left < length && right < length && left < right &&
         right - left >= LinearNussinov::cached_support_minimum_pair_span;
}

std::vector<std::vector<uint>>
dense_pair_support(uint length)
{
  std::vector<std::vector<uint>> pairs_by_right(length);
  for (uint right = 3; right < length; ++right)
    for (uint left = 0; left + 2 < right; ++left)
      pairs_by_right[right].push_back(left);
  return pairs_by_right;
}

void
add_sparse_pair_support(const SparseFloatMatrix& matrix,
                        std::vector<std::vector<uint>>& pairs_by_right)
{
  const uint length = pairs_by_right.size();
  assert(matrix.rows() == length && matrix.columns() == length);
  // This pass only collects support; sorting each sparse row here adds an
  // unnecessary O(m log m) cost before the right-endpoint grouping.
  matrix.for_each_nonzero([&](uint left, uint right, float value) {
    if (value != 0.0f && left + 2 < right && right < length)
      pairs_by_right[right].push_back(left);
  });
}

void
deduplicate_pair_support(std::vector<std::vector<uint>>& pairs_by_right)
{
  std::vector<uint> scratch;
  for (auto& lefts : pairs_by_right) {
    // Pair support is already grouped by right endpoint.  Radix sorting the
    // left indices keeps this union/deduplication pass O(M) for M sparse
    // entries (uint has fixed width), including a row with O(L) candidates.
    radix_sort_linear_nussinov_starts(lefts, scratch);
    lefts.erase(std::unique(lefts.begin(), lefts.end()), lefts.end());
  }
}

struct LinearNussinovCertificate
{
  std::vector<std::vector<double>> prefix_potentials;
  std::vector<double> left;
  std::vector<double> right;
  std::vector<double> incident;
  std::vector<double> half_incident;
  std::vector<double> potential;
  std::vector<double> tightened_forward;
  std::vector<double> tightened_reverse;
  std::vector<std::vector<std::pair<uint, double>>> adjacency;
  double additive_upper_bound = 0.0;
  double pruned_upper_bound = -std::numeric_limits<double>::infinity();
  size_t pruned_states = 0;

  double completion_bound(uint start, double score,
                          uint right, uint length) const
  {
    double outside = std::numeric_limits<double>::infinity();
    for (const auto& prefix : prefix_potentials) {
      outside = std::min(outside,
          prefix[start] + prefix[length] - prefix[right+1]);
    }
    return score + outside;
  }

  void record(uint start, double score, uint right, uint length)
  {
    pruned_upper_bound = std::max(
        pruned_upper_bound, completion_bound(start, score, right, length));
    ++pruned_states;
  }
};

void
tightened_vertex_cover(
    const std::vector<std::vector<std::pair<uint, double>>>& adjacency,
    const std::vector<double>& initial, bool reverse_first,
    std::vector<double>& potential, std::vector<double>& best)
{
  potential = initial;
  best = initial;
  double best_total = std::accumulate(best.begin(), best.end(), 0.0);
  constexpr uint sweeps = 4;
  for (uint sweep = 0; sweep < sweeps; ++sweep) {
    const bool reverse = reverse_first != ((sweep & 1u) != 0);
    for (uint step = 0; step < potential.size(); ++step) {
      const uint i = reverse
          ? static_cast<uint>(potential.size()) - 1 - step : step;
      double tightened = 0.0;
      for (const auto& [j, value] : adjacency[i])
        tightened = std::max(tightened, value - potential[j]);
      potential[i] = tightened;
    }
    const double total =
        std::accumulate(potential.begin(), potential.end(), 0.0);
    if (total < best_total) {
      best_total = total;
      best = potential;
    }
  }
}

void
make_prefix_potential(const std::vector<double>& cover,
                      std::vector<double>& prefix)
{
  prefix.resize(cover.size()+1);
  prefix[0] = 0.0;
  for (size_t i = 0; i < cover.size(); ++i)
    prefix[i+1] = prefix[i] + cover[i];
}

void
prepare_linear_nussinov_certificate(
    uint length, const std::vector<std::vector<uint>>& pairs_by_right,
    const std::vector<float>& pair_scores,
    const std::vector<size_t>& pair_score_offsets,
    LinearNussinovCertificate& certificate)
{
  certificate.left.assign(length, 0.0);
  certificate.right.assign(length, 0.0);
  certificate.incident.assign(length, 0.0);
  certificate.half_incident.resize(length);
  certificate.adjacency.resize(length);
  for (auto& row : certificate.adjacency)
    row.clear();
  certificate.prefix_potentials.resize(5);
  certificate.additive_upper_bound =
      std::numeric_limits<double>::infinity();
  certificate.pruned_upper_bound =
      -std::numeric_limits<double>::infinity();
  certificate.pruned_states = 0;

  for (uint j = 0; j < length; ++j) {
    for (size_t support_index = 0;
         support_index < pairs_by_right[j].size(); ++support_index) {
      const uint i = pairs_by_right[j][support_index];
      if (!linear_nussinov_pair_is_usable(i, j, length))
        continue;
      const double value = std::max(0.0,
          static_cast<double>(pair_scores[pair_score_offsets[j] +
                                         support_index]));
      if (!(value > 0.0))
        continue;
      certificate.left[i] = std::max(certificate.left[i], value);
      certificate.right[j] = std::max(certificate.right[j], value);
      certificate.incident[i] = std::max(certificate.incident[i], value);
      certificate.incident[j] = std::max(certificate.incident[j], value);
      certificate.adjacency[i].push_back({j, value});
      certificate.adjacency[j].push_back({i, value});
    }
  }

  for (uint i = 0; i < length; ++i)
    certificate.half_incident[i] = 0.5 * certificate.incident[i];

  make_prefix_potential(certificate.left,
                        certificate.prefix_potentials[0]);
  make_prefix_potential(certificate.right,
                        certificate.prefix_potentials[1]);
  make_prefix_potential(certificate.half_incident,
                        certificate.prefix_potentials[2]);
  tightened_vertex_cover(certificate.adjacency,
                          certificate.half_incident, false,
                          certificate.potential,
                          certificate.tightened_forward);
  make_prefix_potential(certificate.tightened_forward,
                        certificate.prefix_potentials[3]);
  tightened_vertex_cover(certificate.adjacency,
                          certificate.half_incident, true,
                          certificate.potential,
                          certificate.tightened_reverse);
  make_prefix_potential(certificate.tightened_reverse,
                        certificate.prefix_potentials[4]);

#ifndef NDEBUG
  const std::vector<const std::vector<double>*> covers = {
      &certificate.left, &certificate.right, &certificate.half_incident,
      &certificate.tightened_forward, &certificate.tightened_reverse};
  for (const auto* cover : covers)
    for (uint i = 0; i < length; ++i)
      for (const auto& [j, value] : certificate.adjacency[i])
        assert((*cover)[i] + (*cover)[j] + 1e-10 >= value);
#endif

  for (const auto& prefix : certificate.prefix_potentials) {
    certificate.additive_upper_bound = std::min(
        certificate.additive_upper_bound, prefix.back());
  }
  if (length == 0)
    certificate.additive_upper_bound = 0.0;
}
} // namespace

// Restored from d7f4b07, before IPknot Nussinov-kernel integration.
// Preserve float addition order, skip-left tie precedence, and the original
// candidate-list recurrence; this is the default nonlinear CLI backend.
namespace {
template <typename PairScore>
float decode_dafs_stacked_nussinov(uint length, PairScore pair_score, VU& ss)
{
  ss.assign(length, -1u);
  if (length == 0)
    return 0.0f;
  const double missing = -std::numeric_limits<double>::infinity();
  std::vector<std::vector<double>> complete(length, std::vector<double>(length));
  std::vector<std::vector<double>> paired(length, std::vector<double>(length, missing));
  VVU trace(length, VU(length, -1u));
  std::vector<std::vector<bool>> extend(length, std::vector<bool>(length));
  for (uint right = 0; right < length; ++right) {
    std::vector<std::pair<uint, double>> blocks;
    for (uint end = right + 1; end > 0; --end) {
      const uint left = end - 1;
      complete[left][right] = left < right ? complete[left][right-1] : 0.0;
      if (left + 2 < right) {
        const double local = pair_score(left, right);
        if (std::isfinite(local)) {
          const double inner_pair = paired[left+1][right-1];
          extend[left][right] = inner_pair > complete[left+1][right-1];
          paired[left][right] = local + std::max(
              inner_pair, complete[left+1][right-1]);
          const double block = local + inner_pair;
          if (block > 0.0)
            blocks.emplace_back(left, block);
        }
      }
      for (const auto [start, block] : blocks) {
        const double prefix = start > left ? complete[left][start-1] : 0.0;
        if (prefix + block > complete[left][right]) {
          complete[left][right] = prefix + block;
          trace[left][right] = start;
        }
      }
    }
  }
  std::vector<std::tuple<uint, uint, bool>> pending{{0, length-1, false}};
  while (!pending.empty()) {
    const auto [left, right, is_paired] = pending.back();
    pending.pop_back();
    if (is_paired) {
      ss[left] = right;
      pending.emplace_back(left+1, right-1, extend[left][right]);
    } else if (trace[left][right] != -1u) {
      const uint start = trace[left][right];
      ss[start] = right;
      if (left < start)
        pending.emplace_back(left, start-1, false);
      pending.emplace_back(start+1, right-1, true);
    } else if (left < right) {
      pending.emplace_back(left, right-1, false);
    }
  }
  return static_cast<float>(complete[0][length-1]);
}
} // namespace

DAFSNussinov::DAFSNussinov(float th, bool allow_isolated, int dp_threads)
    : th_(th), allow_isolated_(allow_isolated)
{
  // Retain option validation, while the original DAFS kernel stays serial.
  validate_dp_threads(dp_threads);
}

float
DAFSNussinov::
score(float w, const VVF& p, const VVF& q, const VU& ss)
{
  if (p.empty()) return 0.0f;
  uint L=p.size();
  assert(p[0].size()==L);
  float score = 0.0;

  for (uint i=0; i!=L; ++i)
  {
    if (ss[i] != -1u)
    {
      score += w * (p[i][ss[i]] - th_) - q[i][ss[i]];
    }
  }

  return score;
}

float
DAFSNussinov::
decode(float w, const VVF& p, const VVF& q, VU& ss)
{
  if (p.empty()) { ss.clear(); return 0.0f; }
  if (!allow_isolated_)
    return decode_dafs_stacked_nussinov(p.size(), [&](uint i, uint j) {
      return nussinov_pair_score(w, p[i][j], q[i][j], th_, false);
    }, ss);
  uint L=p.size();
  assert(p[0].size()==L);

  VVF dp(L, VF(L, 0.0));
  BP bp(L);
  VVU tr(L, VU(L, 0));
  for (uint l=1; l<L; ++l)
  {
    for (uint i=0; i+l<L; ++i)
    {
      uint j=i+l;
      float v=0.0;
      int t=0;
      if (i+1<j)
      {
        v=dp[i+1][j];
        t=1;
      }
      if (i<j-1 && v<dp[i][j-1])
      {
        v=dp[i][j-1];
        t=2;
      }
      if (i+1<j-1)
      {
        float s=w*(p[i][j]-th_)-q[i][j];
        if (s>0.0)
        {
          bp[j].push_back(std::make_pair(i,dp[i+1][j-1]+s));
          if (v<dp[i+1][j-1]+s)
          {
            v=dp[i+1][j-1]+s;
            t=3;
          }
        }
      }
      for (SV::const_iterator x=bp[j].begin(); x!=bp[j].end(); ++x)
      {
        const uint k=x->first;
        const float s=x->second;
        if (i<k)
        {
          if (v<dp[i][k-1]+s)
          {
            v=dp[i][k-1]+s;
            t=k-i+3;
          }
        }
      }
      dp[i][j]=v;
      tr[i][j]=t;
    }
  }

  // trace back
  ss.resize(L);
  std::fill(ss.begin(), ss.end(), -1u);
  std::stack<std::pair<uint,uint> > st;
  st.push(std::make_pair(0, L-1));
  while (!st.empty())
  {
    const std::pair<uint,uint> p=st.top(); st.pop();
    const int i=p.first, j=p.second;
    switch (tr[i][j])
    {
      case 0:
        break;
      case 1:
        st.push(std::make_pair(i+1, j));
        break;
      case 2:
        st.push(std::make_pair(i, j-1));
        break;
      case 3:
        ss[i]=j;
        st.push(std::make_pair(i+1, j-1));
        break;
      default:
        const int k=i+tr[i][j]-3;
        st.push(std::make_pair(i, k-1));
        ss[k]=j;
        st.push(std::make_pair(k+1, j-1));
        break;
    }
  }

  return dp[0][L-1];
}

float
DAFSNussinov::
decode(float w, const VVF& p, const SparseFloatMatrix& q, VU& ss)
{
  if (p.empty()) { ss.clear(); return 0.0f; }
  if (!allow_isolated_)
    return decode_dafs_stacked_nussinov(p.size(), [&](uint i, uint j) {
      return nussinov_pair_score(w, p[i][j], q.get(i, j), th_, false);
    }, ss);
  uint L=p.size();
  assert(p[0].size()==L);
  assert(q.rows()==L && q.columns()==L);

  VVF dp(L, VF(L, 0.0));
  BP bp(L);
  VVU tr(L, VU(L, 0));
  std::vector<const SV*> multiplier_rows(L);
  std::vector<size_t> multiplier_positions(L, 0);
  for (uint i=0; i<L; ++i)
    multiplier_rows[i] = &q.ordered_row(i);
  for (uint l=1; l<L; ++l)
  {
    for (uint i=0; i+l<L; ++i)
    {
      uint j=i+l;
      float v=0.0;
      int t=0;
      if (i+1<j)
      {
        v=dp[i+1][j];
        t=1;
      }
      if (i<j-1 && v<dp[i][j-1])
      {
        v=dp[i][j-1];
        t=2;
      }
      if (i+1<j-1)
      {
        // Multipliers can become nonzero only on a selected probability
        // candidate or a CBP projection, both of which have p>0.  Avoid a
        // multiplier lookup for the zero cells of a LinearFold matrix.
        float multiplier = 0.0f;
        if (p[i][j] > 0.0f)
        {
          const SV& row = *multiplier_rows[i];
          size_t& pos = multiplier_positions[i];
          while (pos < row.size() && row[pos].first < j) ++pos;
          if (pos < row.size() && row[pos].first == j)
            multiplier = row[pos].second;
        }
        float s=w*(p[i][j]-th_)-multiplier;
        if (s>0.0)
        {
          bp[j].push_back(std::make_pair(i,dp[i+1][j-1]+s));
          if (v<dp[i+1][j-1]+s)
          {
            v=dp[i+1][j-1]+s;
            t=3;
          }
        }
      }
      for (SV::const_iterator x=bp[j].begin(); x!=bp[j].end(); ++x)
      {
        const uint k=x->first;
        const float s=x->second;
        if (i<k)
        {
          if (v<dp[i][k-1]+s)
          {
            v=dp[i][k-1]+s;
            t=k-i+3;
          }
        }
      }
      dp[i][j]=v;
      tr[i][j]=t;
    }
  }

  ss.resize(L);
  std::fill(ss.begin(), ss.end(), -1u);
  std::stack<std::pair<uint,uint> > st;
  st.push(std::make_pair(0, L-1));
  while (!st.empty())
  {
    const std::pair<uint,uint> p=st.top(); st.pop();
    const int i=p.first, j=p.second;
    switch (tr[i][j])
    {
      case 0: break;
      case 1: st.push(std::make_pair(i+1, j)); break;
      case 2: st.push(std::make_pair(i, j-1)); break;
      case 3:
        ss[i]=j;
        st.push(std::make_pair(i+1, j-1));
        break;
      default:
        const int k=i+tr[i][j]-3;
        st.push(std::make_pair(i, k-1));
        ss[k]=j;
        st.push(std::make_pair(k+1, j-1));
        break;
    }
  }

  return dp[0][L-1];
}

float
DAFSNussinov::
decode(const VVF& p, VU& ss, std::string& str)
{
  if (p.empty()) { ss.clear(); str.clear(); return 0.0f; }
  if (!allow_isolated_) {
    const float score = decode_dafs_stacked_nussinov(p.size(), [&](uint i, uint j) {
      return nussinov_pair_score(1.0f, p[i][j], 0.0f, th_, false);
    }, ss);
    make_brackets(ss, str);
    return score;
  }
  uint L=p.size();
  assert(p[0].size()==L);

  VVF dp(L, VF(L, 0.0));
  BP bp(L);
  VVU tr(L, VU(L, 0));
  for (uint l=1; l<L; ++l)
  {
    for (uint i=0; i+l<L; ++i)
    {
      uint j=i+l;
      float v=0.0;
      int t=0;
      if (i+1<j)
      {
        v=dp[i+1][j];
        t=1;
      }
      if (i<j-1 && v<dp[i][j-1])
      {
        v=dp[i][j-1];
        t=2;
      }
      if (i+1<j-1)
      {
        float s=p[i][j]-th_;
        if (s>0.0)
        {
          bp[j].push_back(std::make_pair(i,dp[i+1][j-1]+s));
          if (v<dp[i+1][j-1]+s)
          {
            v=dp[i+1][j-1]+s;
            t=3;
          }
        }
      }
      for (SV::const_iterator x=bp[j].begin(); x!=bp[j].end(); ++x)
      {
        const uint k=x->first;
        const float s=x->second;
        if (i<k)
        {
          if (v<dp[i][k-1]+s)
          {
            v=dp[i][k-1]+s;
            t=k-i+3;
          }
        }
      }
      dp[i][j]=v;
      tr[i][j]=t;
    }
  }

  // trace back
  ss.resize(L);
  std::fill(ss.begin(), ss.end(), -1u);
  std::stack<std::pair<uint,uint> > st;
  st.push(std::make_pair(0, L-1));
  while (!st.empty())
  {
    const std::pair<uint,uint> p=st.top(); st.pop();
    const int i=p.first, j=p.second;
    switch (tr[i][j])
    {
      case 0:
        break;
      case 1:
        st.push(std::make_pair(i+1, j));
        break;
      case 2:
        st.push(std::make_pair(i, j-1));
        break;
      case 3:
        ss[i]=j;
        st.push(std::make_pair(i+1, j-1));
        break;
      default:
        const int k=i+tr[i][j]-3;
        st.push(std::make_pair(i, k-1));
        ss[k]=j;
        st.push(std::make_pair(k+1, j-1));
        break;
    }
  }

  make_brackets(ss, str);
  return dp[0][L-1];
}

float
DAFSNussinov::
decode(const SparseFloatMatrix& p, const SparseFloatMatrix& pair_bonus,
       VU& ss, std::string& str)
{
  assert(p.rows() == pair_bonus.rows() &&
         p.columns() == pair_bonus.columns());
  if (!allow_isolated_) {
    const float score = decode_dafs_stacked_nussinov(p.rows(), [&](uint i, uint j) {
      return nussinov_pair_score(1.0f, p.get(i, j), -pair_bonus.get(i, j), th_, false);
    }, ss);
    make_brackets(ss, str);
    return score;
  }
  VVF adjusted = p.dense();
  for (uint i = 0; i < pair_bonus.rows(); ++i)
    for (const auto [j, value] : pair_bonus.ordered_row(i))
      adjusted[i][j] += value;
  return decode(adjusted, ss, str);
}

void
DAFSNussinov::
make_brackets(const VU& ss, std::string& str) const
{
  str=::make_brackets(ss);
}

struct SparseNussinov::Workspace
{
  uint length = 0;
  bool dense_chart = false;
  // DDNussinov borrows this vector; reset_pairs() follows endpoint changes.
  std::vector<DDPair> pairs, next_pairs;
  std::vector<double> weights;
  std::vector<unsigned char> allowed;
  std::vector<int> selected;
  std::unique_ptr<DDNussinov> decoder;
};

SparseNussinov::SparseNussinov(float th, bool allow_isolated, int dp_threads)
    : th_(th), allow_isolated_(allow_isolated), dp_threads_(dp_threads),
      workspace_(std::make_unique<Workspace>())
{
  validate_dp_threads(dp_threads);
}

SparseNussinov::~SparseNussinov() = default;

int SparseNussinov::used_workers() const
{
  return workspace_->decoder ? workspace_->decoder->used_workers() : 1;
}

size_t SparseNussinov::parallel_steps() const
{
  return workspace_->decoder ? workspace_->decoder->parallel_steps() : 0;
}

size_t SparseNussinov::sparsified_candidates() const
{
  return workspace_->decoder ? workspace_->decoder->sparsified_candidates() : 0;
}

template <typename PairScore>
float SparseNussinov::decode_impl(
    uint length, const std::vector<std::vector<uint>>* support,
    PairScore pair_score, VU& ss)
{
  if (length > static_cast<uint>(std::numeric_limits<int>::max()))
    throw std::invalid_argument("Nussinov sequence is too long");
  auto& w = *workspace_;
  w.next_pairs.clear();
  w.weights.clear();
  const auto append = [&](uint left, uint right) {
    // Preserve DAFS's minimum span (three) rather than IPknot's span one.
    if (left >= right || right - left < 3) return;
    const float local = pair_score(left, right);
    if (local == -std::numeric_limits<float>::infinity()) return;
    if (!std::isfinite(local))
      throw std::invalid_argument("Nussinov pair scores must be finite");
    if (allow_isolated_ && !(local > 0.0f)) return;
    w.next_pairs.push_back({static_cast<int>(left), static_cast<int>(right), 0, 0.0});
    w.weights.push_back(local);
  };
  if (support) {
    for (uint right = 0; right < length; ++right)
      for (uint left : (*support)[right]) append(left, right);
  } else {
    // Scan dense input in row order, without allocating a quadratic endpoint
    // list. Only retained score candidates enter the reusable sparse oracle.
    for (uint left = 0; left < length; ++left)
      for (uint right = left+3; right < length; ++right) append(left, right);
  }
  if (w.next_pairs.size() > static_cast<size_t>(std::numeric_limits<int>::max()))
    throw std::invalid_argument("Too many Nussinov pairs");
  // A populated sparse graph can also benefit from the compact interval
  // chart. Bound its extra storage by 512 bytes per retained pair; extremely
  // sparse, long inputs keep the sparse chart without a quadratic allocation.
  const bool dense_chart = support == nullptr ||
      6.0 * length * (double(length) + 1.0) <= 512.0 * w.next_pairs.size();
  const bool same_graph = w.length == length && w.dense_chart == dense_chart &&
      w.pairs.size() == w.next_pairs.size() &&
      std::equal(w.pairs.begin(), w.pairs.end(), w.next_pairs.begin(),
          [](const DDPair& a, const DDPair& b) {
            return a.left == b.left && a.right == b.right;
          });
  if (!w.decoder || !same_graph) {
    auto decoder = std::move(w.decoder);
    if (w.length != length || w.dense_chart != dense_chart) decoder.reset();
    w.pairs.swap(w.next_pairs);
    w.length = length;
    w.dense_chart = dense_chart;
    // Reuse the public exact oracle, including closed-candidate sparsification,
    // suffix dominance, signed stacking states and lazy traceback allocation.
    // Finite-beam decoding below uses the same oracle with exact pruning
    // disabled and a pruning observer for DAFS completion certificates.
    // Dense inputs use IPknot's compact interval chart. Sparse inputs use
    // its sparse exact chart, including long sequences with few candidates.
    // Both charts apply the same exact closed-candidate sparsification.
    if (decoder) decoder->reset_pairs();
    else decoder = std::make_unique<DDNussinov>(
        length, w.pairs, 0, 0, !allow_isolated_, true, dp_threads_, dense_chart);
    w.decoder = std::move(decoder);
  }
  w.allowed.assign(w.pairs.size(), 1);
  const double score = w.decoder->decode(w.weights, w.allowed, w.selected);
  ss.assign(length, -1u);
  for (int id : w.selected)
    ss[w.pairs[id].left] = w.pairs[id].right;
  return static_cast<float>(score);
}

float SparseNussinov::decode(float w, const VVF& p, const VVF& q, VU& ss)
{
  return decode_impl(p.size(), nullptr, [&](uint i, uint j) {
    return nussinov_pair_score(w, p[i][j], q[i][j], th_, allow_isolated_);
  }, ss);
}

float SparseNussinov::decode(
    float w, const VVF& p, const SparseFloatMatrix& q, VU& ss)
{
  assert(q.rows() == p.size() && q.columns() == p.size());
  return decode_impl(p.size(), nullptr, [&](uint i, uint j) {
    return nussinov_pair_score(w, p[i][j], q.get(i, j), th_, allow_isolated_);
  }, ss);
}

float SparseNussinov::decode(
    float w, const SparseFloatMatrix& p, const VVF& q, VU& ss)
{
  return decode_impl(p.rows(), nullptr, [&](uint i, uint j) {
    return nussinov_pair_score(w, p.get(i, j), q[i][j], th_, allow_isolated_);
  }, ss);
}

float SparseNussinov::decode(
    float w, const SparseFloatMatrix& p, const SparseFloatMatrix& q, VU& ss)
{
  const uint length = p.rows();
  assert(p.columns() == length && q.rows() == length && q.columns() == length);
  std::vector<std::vector<uint>> support(length);
  // In this unusual scoring regime even absent entries are profitable.
  const bool dense = allow_isolated_ && -w * th_ > 0.0f;
  if (!dense) {
    add_sparse_pair_support(p, support);
    add_sparse_pair_support(q, support);
    deduplicate_pair_support(support);
  }
  return decode_impl(length, dense ? nullptr : &support, [&](uint i, uint j) {
    return nussinov_pair_score(w, p.get(i, j), q.get(i, j), th_, allow_isolated_);
  }, ss);
}

float SparseNussinov::decode(const VVF& p, VU& ss, std::string& str)
{
  const float score = decode_impl(p.size(), nullptr,
      [&](uint i, uint j) {
        return nussinov_pair_score(1.0f, p[i][j], 0.0f, th_, allow_isolated_);
      }, ss);
  make_brackets(ss, str);
  return score;
}

float SparseNussinov::decode(const SparseFloatMatrix& p, VU& ss, std::string& str)
{
  SparseFloatMatrix bonus;
  bonus.assign(p.rows(), p.columns());
  return decode(p, bonus, ss, str);
}

float SparseNussinov::decode(const SparseFloatMatrix& p,
    const SparseFloatMatrix& pair_bonus, VU& ss, std::string& str)
{
  const uint length = p.rows();
  assert(p.columns() == length && pair_bonus.rows() == length &&
         pair_bonus.columns() == length);
  std::vector<std::vector<uint>> support(length);
  const bool dense = allow_isolated_ && th_ < 0.0f;
  if (!dense) {
    add_sparse_pair_support(p, support);
    add_sparse_pair_support(pair_bonus, support);
    deduplicate_pair_support(support);
  }
  const float score = decode_impl(length, dense ? nullptr : &support, [&](uint i, uint j) {
    return nussinov_pair_score(1.0f, p.get(i, j), -pair_bonus.get(i, j),
                             th_, allow_isolated_);
  }, ss);
  make_brackets(ss, str);
  return score;
}

float SparseNussinov::score(float w, const VVF& p, const VVF& q, const VU& ss)
{
  float score = 0.0f;
  for (uint i = 0; i < ss.size(); ++i)
    if (ss[i] != -1u) score += w * (p[i][ss[i]] - th_) - q[i][ss[i]];
  return score;
}

void SparseNussinov::make_brackets(const VU& ss, std::string& str) const
{
  str = ::make_brackets(ss);
}

struct LinearNussinov::Workspace
{
  bool used_ipknot = false;
  uint length = 0;
  std::vector<std::vector<uint>> support;
  // The public IPknot oracle borrows this stable graph. Rebuild it only when
  // endpoints change; weights and masks change independently on every solve.
  std::vector<DDPair> pairs;
  std::vector<double> weights;
  std::vector<unsigned char> allowed;
  std::vector<int> selected;
  std::vector<float> pair_scores;
  std::vector<size_t> pair_score_offsets;
  std::unique_ptr<DDNussinov> decoder;
  int used_workers = 1;
  size_t parallel_steps = 0;
  LinearNussinovCertificate certificate;
};

LinearNussinov::LinearNussinov(float th, uint beam_size, bool allow_isolated,
                               int dp_threads)
    : th_(th), beam_size_(beam_size), allow_isolated_(allow_isolated),
      dp_threads_(dp_threads),
      workspace_(std::make_unique<Workspace>())
{
  validate_dp_threads(dp_threads);
}

LinearNussinov::~LinearNussinov() = default;

int LinearNussinov::used_workers() const { return workspace_->used_workers; }
size_t LinearNussinov::parallel_steps() const { return workspace_->parallel_steps; }
bool LinearNussinov::used_ipknot_kernel() const { return workspace_->used_ipknot; }

template <typename PairScore>
LinearNussinovResult
LinearNussinov::
decode_impl(uint L, const std::vector<std::vector<uint>>& pairs_by_right,
            PairScore pair_score, VU& ss, bool certify)
{
  Workspace& workspace = *workspace_;
  workspace.used_ipknot = false;
  workspace.used_workers = 1;
  workspace.parallel_steps = 0;
  if (L == 0) { ss.clear(); return {}; }
  if (L > static_cast<uint>(std::numeric_limits<int>::max()) ||
      beam_size_ > static_cast<uint>(std::numeric_limits<int>::max()))
    throw std::invalid_argument("LinearNussinov dimensions exceed the decoder limit");

  size_t score_count = 0;
  for (const auto& row : pairs_by_right) score_count += row.size();
  // Always use IPknot improved beam, including short and very sparse inputs.
  workspace.used_ipknot = true;
  ss.assign(L, -1u);
  workspace.pair_score_offsets.resize(L+1);
  size_t offset = 0;
  for (uint right = 0; right < L; ++right) {
    workspace.pair_score_offsets[right] = offset;
    offset += pairs_by_right[right].size();
  }
  workspace.pair_score_offsets[L] = offset;
  if (!workspace.decoder || workspace.length != L || workspace.support != pairs_by_right) {
    // Keep the chart buffers across endpoint changes. Moving ownership out
    // also leaves a clean rebuild path if graph construction throws.
    auto decoder = std::move(workspace.decoder);
    if (workspace.length != L) decoder.reset();
    workspace.length = L;
    workspace.support = pairs_by_right;
    workspace.pairs.clear();
    for (uint right = 0; right < L; ++right)
      for (uint left : pairs_by_right[right])
        if (linear_nussinov_pair_is_usable(left, right, L))
          workspace.pairs.push_back({int(left), int(right), 0, 0.0});
    if (workspace.pairs.size() > static_cast<size_t>(std::numeric_limits<int>::max()))
      throw std::invalid_argument("Too many LinearNussinov pairs");
    // IPknot supplies the faster beam, signed stacking and traceback kernels.
    // Native improved beam may use an empty interior when all inside states
    // were pruned. The observer retains DAFS completion certificates.
    if (decoder) decoder->reset_pairs();
    else decoder = std::make_unique<DDNussinov>(
        L, workspace.pairs, 0, std::max(1u, beam_size_), !allow_isolated_,
        true, dp_threads_, false);
    workspace.decoder = std::move(decoder);
  }
  workspace.weights.resize(workspace.pairs.size());
  workspace.allowed.resize(workspace.pairs.size());
  if (certify) workspace.pair_scores.assign(score_count, 0.0f);
  size_t id = 0;
  for (uint right = 0; right < L; ++right) {
    const auto& row = pairs_by_right[right];
    for (size_t k = 0; k < row.size(); ++k) {
      const uint left = row[k];
      if (!linear_nussinov_pair_is_usable(left, right, L)) continue;
      const float value = pair_score(left, right);
      if (certify)
        workspace.pair_scores[workspace.pair_score_offsets[right]+k] = value;
      workspace.weights[id] = std::isfinite(value) ? value : 0.0;
      workspace.allowed[id] = std::isfinite(value) && (!allow_isolated_ || value > 0);
      ++id;
    }
  }
  auto& certificate = workspace.certificate;
  if (certify)
    prepare_linear_nussinov_certificate(L, pairs_by_right, workspace.pair_scores,
                                       workspace.pair_score_offsets, certificate);
  DDNussinov::PruningObserver observer{nullptr, nullptr};
  if (certify) {
    observer.context = &workspace;
    observer.record = [](void* context, int start, int right, double score) {
      auto& w = *static_cast<Workspace*>(context);
      w.certificate.record(start, score, right, w.length);
    };
  }
  const double score = workspace.decoder->decode(
      workspace.weights, workspace.allowed, workspace.selected, observer);
  workspace.used_workers = workspace.decoder->used_workers();
  workspace.parallel_steps = workspace.decoder->parallel_steps();
  for (int selected : workspace.selected) {
    const auto& pair = workspace.pairs[selected];
    ss[pair.left] = pair.right;
  }
  LinearNussinovResult result;
  result.score = static_cast<float>(score);
  if (!certify) {
    result.upper_bound = result.score;
    result.pruned_upper_bound = result.score;
    result.additive_upper_bound = result.score;
    return result;
  }
  const double pruned_upper_bound = std::max(score, certificate.pruned_upper_bound);
  const double upper_bound = std::min(pruned_upper_bound, certificate.additive_upper_bound);
  result.pruned_upper_bound = std::max(result.score,
      RelaxedBounds::round_up_to_float(pruned_upper_bound));
  result.additive_upper_bound = std::max(result.score,
      RelaxedBounds::round_up_to_float(certificate.additive_upper_bound));
  result.upper_bound = std::max(result.score, RelaxedBounds::round_up_to_float(upper_bound));
  result.pruned_states = certificate.pruned_states;
  return result;
}

float
LinearNussinov::
decode(float w, const VVF& p, const VVF& q, VU& ss)
{
  const uint length = p.size();
  return decode_impl(length, dense_pair_support(length),
                     [&](uint i, uint j) {
                       return nussinov_pair_score(w, p[i][j], q[i][j], th_, allow_isolated_);
                     }, ss, false).score;
}

float
LinearNussinov::
decode(float w, const VVF& p, const SparseFloatMatrix& q, VU& ss)
{
  const uint length = p.size();
  return decode_impl(length, dense_pair_support(length),
                     [&](uint i, uint j) {
                       return nussinov_pair_score(w, p[i][j], q.get(i, j), th_, allow_isolated_);
                     }, ss, false).score;
}

float
LinearNussinov::
decode(float w, const SparseFloatMatrix& p, const VVF& q, VU& ss)
{
  const uint length = p.rows();
  return decode_impl(length, dense_pair_support(length),
                     [&](uint i, uint j) {
                       return nussinov_pair_score(w, p.get(i, j), q[i][j], th_, allow_isolated_);
                     }, ss, false).score;
}

float
LinearNussinov::
decode(float w, const SparseFloatMatrix& p, const SparseFloatMatrix& q, VU& ss)
{
  const uint length = p.rows();
  std::vector<std::vector<uint>> support(length);
  add_sparse_pair_support(p, support);
  add_sparse_pair_support(q, support);
  deduplicate_pair_support(support);
  return decode_impl(length, support,
                     [&](uint i, uint j) {
                       return nussinov_pair_score(w, p.get(i, j), q.get(i, j), th_, allow_isolated_);
                     }, ss, false).score;
}

LinearNussinovResult
LinearNussinov::
decode_certified(float w, const SparseFloatMatrix& p,
                 const VVF& q, VU& ss)
{
  const uint length = p.rows();
  return decode_impl(length, dense_pair_support(length),
                     [&](uint i, uint j) {
                       return nussinov_pair_score(w, p.get(i, j), q[i][j], th_, allow_isolated_);
                     }, ss, true);
}

LinearNussinovResult
LinearNussinov::
decode_certified(float w, const SparseFloatMatrix& p,
                 const SparseFloatMatrix& q, VU& ss)
{
  const uint length = p.rows();
  std::vector<std::vector<uint>> support(length);
  add_sparse_pair_support(p, support);
  add_sparse_pair_support(q, support);
  deduplicate_pair_support(support);
  return decode_impl(length, support,
                     [&](uint i, uint j) {
                       return nussinov_pair_score(w, p.get(i, j), q.get(i, j), th_, allow_isolated_);
                     }, ss, true);
}

LinearNussinovResult
LinearNussinov::
decode_certified(
    float w, const SparseFloatMatrix& p, const VVF& q,
    const std::vector<std::vector<uint>>& pairs_by_right, VU& ss)
{
  const uint length = p.rows();
  assert(pairs_by_right.size() == length);
  return decode_impl(length, pairs_by_right,
                     [&](uint i, uint j) {
                       return nussinov_pair_score(w, p.get(i, j), q[i][j], th_, allow_isolated_);
                     }, ss, true);
}

LinearNussinovResult
LinearNussinov::
decode_certified(
    float w, const SparseFloatMatrix& p, const SparseFloatMatrix& q,
    const std::vector<std::vector<uint>>& pairs_by_right, VU& ss)
{
  const uint length = p.rows();
  assert(pairs_by_right.size() == length);
  return decode_impl(length, pairs_by_right,
                     [&](uint i, uint j) {
                       return nussinov_pair_score(w, p.get(i, j), q.get(i, j), th_, allow_isolated_);
                     }, ss, true);
}

float
LinearNussinov::
decode(const VVF& p, VU& ss, std::string& str)
{
  const uint length = p.size();
  const float score = decode_impl(
      length, dense_pair_support(length),
      [&](uint i, uint j) { return nussinov_pair_score(1.0f, p[i][j], 0.0f, th_, allow_isolated_); }, ss, false).score;
  make_brackets(ss, str);
  return score;
}

float
LinearNussinov::
decode(const SparseFloatMatrix& p, VU& ss, std::string& str)
{
  const uint length = p.rows();
  std::vector<std::vector<uint>> support(length);
  add_sparse_pair_support(p, support);
  const float score = decode_impl(
      length, support,
      [&](uint i, uint j) { return nussinov_pair_score(1.0f, p.get(i, j), 0.0f, th_, allow_isolated_); }, ss, false).score;
  make_brackets(ss, str);
  return score;
}

float
LinearNussinov::
decode(const SparseFloatMatrix& p, const SparseFloatMatrix& pair_bonus,
       VU& ss, std::string& str)
{
  const uint length = p.rows();
  assert(p.columns() == length && pair_bonus.rows() == length &&
         pair_bonus.columns() == length);
  std::vector<std::vector<uint>> support(length);
  add_sparse_pair_support(p, support);
  add_sparse_pair_support(pair_bonus, support);
  deduplicate_pair_support(support);
  const float score = decode_impl(
      length, support,
      [&](uint i, uint j) {
        return nussinov_pair_score(1.0f, p.get(i, j), -pair_bonus.get(i, j), th_, allow_isolated_);
      }, ss, false).score;
  make_brackets(ss, str);
  return score;
}

void
LinearNussinov::
make_brackets(const VU& ss, std::string& str) const
{
  str=::make_brackets(ss);
}

static
std::string
make_brackets(const VU& ss)
{
  std::string s(ss.size(), '.');
  for (uint i=0; i!=ss.size(); ++i)
    if (ss[i]!=-1u)
    {
      s[i]=Fold::Decoder::left_brackets[0];
      s[ss[i]]=Fold::Decoder::right_brackets[0];
    }
  return s;
}
