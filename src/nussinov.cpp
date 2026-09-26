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
#include <array>
#include <algorithm>
#include <cassert>
#include <limits>
#include <numeric>
#include <stack>

static
std::string
make_brackets(const VU& ss);

float
Nussinov::
decode(float w, const VVF& p, const VVF& q, VU& ss)
{
  uint L=p.size();
  assert(p[0].size()==L);

  // calculate scoring matrices for the current step
  VVF sm(L, VF(L, 0.0));
  for (uint i=0; i!=L-1; ++i)
    for (uint j=i+1; j!=L; ++j)
      sm[i][j] = w*(p[i][j]-th_)-q[i][j];

  VVF dp(L, VF(L, 0.0));
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
  uint L=p.size();
  assert(p[0].size()==L);

  // calculate scoring matrices for the current step
  VVF sm(L, VF(L, 0.0));
  for (uint i=0; i!=L-1; ++i)
    for (uint j=i+1; j!=L; ++j)
      sm[i][j] = p[i][j]-th_;

  VVF dp(L, VF(L, 0.0));
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

float
SparseNussinov::
score(float w, const VVF& p, const VVF& q, const VU& ss)
{
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
SparseNussinov::
decode(float w, const VVF& p, const VVF& q, VU& ss)
{
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
SparseNussinov::
decode(float w, const VVF& p, const SparseFloatMatrix& q, VU& ss)
{
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
SparseNussinov::
decode(const VVF& p, VU& ss, std::string& str)
{
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
SparseNussinov::
decode(const SparseFloatMatrix& p, const SparseFloatMatrix& pair_bonus,
       VU& ss, std::string& str)
{
  assert(p.rows() == pair_bonus.rows() &&
         p.columns() == pair_bonus.columns());
  VVF adjusted = p.dense();
  for (uint i = 0; i < pair_bonus.rows(); ++i)
    for (const auto [j, value] : pair_bonus.ordered_row(i))
      adjusted[i][j] += value;
  return decode(adjusted, ss, str);
}

void
SparseNussinov::
make_brackets(const VU& ss, std::string& str) const
{
  str=::make_brackets(ss);
}

namespace {
struct LinearNussinovCell
{
  uint start;
  double score;
  uint pair_left;
  uint inside_start;
  bool paired;
};

using LinearNussinovBeam = std::vector<LinearNussinovCell>;

struct LinearNussinovRankedCell
{
  const LinearNussinovCell* cell;
  double priority;
};

bool
linear_nussinov_cell_is_better(const LinearNussinovRankedCell& lhs,
                               const LinearNussinovRankedCell& rhs)
{
  return lhs.priority != rhs.priority
      ? lhs.priority > rhs.priority
      : lhs.cell->start < rhs.cell->start;
}

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

const LinearNussinovCell*
find_linear_nussinov_cell(const LinearNussinovBeam& beam, uint start)
{
  const auto it = std::lower_bound(
      beam.begin(), beam.end(), start,
      [](const LinearNussinovCell& cell, uint value) {
        return cell.start < value;
      });
  return it != beam.end() && it->start == start ? &*it : nullptr;
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

  double completion_bound(const LinearNussinovCell& cell,
                          uint right, uint length) const
  {
    double outside = std::numeric_limits<double>::infinity();
    for (const auto& prefix : prefix_potentials) {
      outside = std::min(outside,
          prefix[cell.start] + prefix[length] - prefix[right+1]);
    }
    return cell.score + outside;
  }

  void record(const LinearNussinovCell& cell, uint right, uint length)
  {
    pruned_upper_bound = std::max(
        pruned_upper_bound, completion_bound(cell, right, length));
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

struct LinearNussinov::Workspace
{
  std::vector<LinearNussinovBeam> beams;
  std::vector<const LinearNussinovCell*> beam_roots;
  std::vector<LinearNussinovCell> candidates;
  std::vector<uint> candidate_generation;
  std::vector<uint> candidate_starts;
  std::vector<uint> radix_scratch;
  std::vector<uint> dominated_generation;
  std::vector<uint> retained_generation;
  std::vector<const LinearNussinovCell*> suffix_best;
  std::vector<LinearNussinovRankedCell> selected;
  std::vector<std::pair<uint, uint>> pending;
  std::vector<float> pair_scores;
  std::vector<size_t> pair_score_offsets;
  uint generation = 0;
  LinearNussinovCertificate certificate;
};

LinearNussinov::LinearNussinov(float th, uint beam_size)
    : th_(th), beam_size_(beam_size), workspace_(std::make_unique<Workspace>())
{
}

LinearNussinov::~LinearNussinov() = default;

template <typename PairScore>
LinearNussinovResult
LinearNussinov::
decode_impl(uint L, const std::vector<std::vector<uint>>& pairs_by_right,
            PairScore pair_score, VU& ss, bool certify)
{
  ss.assign(L, -1u);
  if (L == 0)
    return {};

  Workspace& workspace = *workspace_;
  LinearNussinovCertificate& certificate = workspace.certificate;
  if (certify) {
    workspace.pair_score_offsets.resize(L + 1);
    size_t score_count = 0;
    for (uint right = 0; right < L; ++right) {
      workspace.pair_score_offsets[right] = score_count;
      score_count += pairs_by_right[right].size();
    }
    workspace.pair_score_offsets[L] = score_count;
    workspace.pair_scores.assign(score_count, 0.0f);
    for (uint right = 0; right < L; ++right) {
      const size_t score_offset = workspace.pair_score_offsets[right];
      for (size_t support_index = 0;
           support_index < pairs_by_right[right].size(); ++support_index) {
        const uint left = pairs_by_right[right][support_index];
        if (linear_nussinov_pair_is_usable(left, right, L))
          workspace.pair_scores[score_offset + support_index] =
              pair_score(left, right);
      }
    }
    prepare_linear_nussinov_certificate(
        L, pairs_by_right, workspace.pair_scores,
        workspace.pair_score_offsets, certificate);
  }

  // D(i,j) is the best non-crossing matching on [i,j].  At each right
  // endpoint retain only the most promising interval starts.  For a caller-
  // supplied support with M entries, fixed beam b gives O(L + M*f(b)); the
  // dense overloads intentionally materialize their L^2 support before this
  // core is entered.
  workspace.beams.resize(L);
  for (auto& beam : workspace.beams)
    beam.clear();
  workspace.beam_roots.resize(L);
  workspace.candidates.resize(L);
  workspace.candidate_generation.resize(L, 0);
  workspace.candidate_starts.reserve(L);
  workspace.radix_scratch.reserve(L);
  workspace.dominated_generation.resize(L, 0);
  workspace.retained_generation.resize(L, 0);
  workspace.selected.reserve(std::min<size_t>(
      std::max(1u, beam_size_), L));
  auto& beams = workspace.beams;
  const uint beam_limit = std::max(1u, beam_size_);

  for (uint right = 0; right < L; ++right) {
    ++workspace.generation;
    if (workspace.generation == 0) {
      std::fill(workspace.candidate_generation.begin(),
                workspace.candidate_generation.end(), 0);
      std::fill(workspace.dominated_generation.begin(),
                workspace.dominated_generation.end(), 0);
      std::fill(workspace.retained_generation.begin(),
                workspace.retained_generation.end(), 0);
      workspace.generation = 1;
    }
    const uint generation = workspace.generation;
    auto& candidates = workspace.candidates;
    auto& candidate_starts = workspace.candidate_starts;
    candidate_starts.clear();

    const auto offer = [&](uint start, double score, bool paired, uint left,
                           uint inside_start) {
      if (workspace.candidate_generation[start] != generation) {
        workspace.candidate_generation[start] = generation;
        candidates[start] = LinearNussinovCell{
            start, score, left, inside_start, paired};
        candidate_starts.push_back(start);
      } else if (score > candidates[start].score) {
        candidates[start] = LinearNussinovCell{
            start, score, left, inside_start, paired};
      }
    };

    if (right > 0) {
      const auto& previous = beams[right-1];
      auto& suffix_best = workspace.suffix_best;
      suffix_best.resize(previous.size());
      const LinearNussinovCell* best = nullptr;
      for (size_t i = previous.size(); i-- > 0; ) {
        const LinearNussinovCell* cell = &previous[i];
        if (!best || cell->score > best->score ||
            (cell->score == best->score && cell->start > best->start))
          best = cell;
        suffix_best[i] = best;
      }

    }

    // Leave the new right endpoint unpaired, including the new singleton
    // interval.  Offering these first makes ties deterministic.
    offer(right, 0.0, false, -1u, -1u);
    if (right > 0)
      for (const auto& previous : beams[right-1])
        offer(previous.start, previous.score, false, -1u, -1u);

    for (size_t support_index = 0;
         support_index < pairs_by_right[right].size(); ++support_index) {
      const uint left = pairs_by_right[right][support_index];
      if (!linear_nussinov_pair_is_usable(left, right, L))
        continue;
      const double local_score = certify
          ? static_cast<double>(workspace.pair_scores[
                workspace.pair_score_offsets[right] + support_index])
          : pair_score(left, right);
      if (!(local_score > 0.0))
        continue;

      // A structure whose first used base is later than left+1 is a valid
      // inside structure with the skipped leading bases left unpaired.  This
      // suffix dominance removes the O(length) family of equivalent empty or
      // leading-unpaired states without losing an exact derivation.
      const LinearNussinovCell* inside = nullptr;
      if (right > 0) {
        const auto& previous = beams[right-1];
        const auto it = std::lower_bound(
            previous.begin(), previous.end(), left + 1,
            [](const LinearNussinovCell& cell, uint minimum_start) {
              return cell.start < minimum_start;
            });
        if (it != previous.end())
          inside = workspace.suffix_best[it - previous.begin()];
      }
      if (!inside)
        continue;

      // The pair can begin the interval (empty prefix), or follow any
      // retained interval ending immediately before its left endpoint.
      offer(left, inside->score + local_score, true, left, inside->start);
      if (left > 0) {
        for (const auto& prefix_interval : beams[left-1])
          offer(prefix_interval.start,
                prefix_interval.score + inside->score + local_score,
                true, left, inside->start);
      }
    }

    // The same start set drives suffix dominance, selection, and certificate
    // recording.  Certificate::record only accumulates a max and a count, so
    // preserving the original offer order would be identity-only overhead.
    constexpr size_t radix_cutoff = 128;
    if (candidate_starts.size() < radix_cutoff)
      std::sort(candidate_starts.begin(), candidate_starts.end());
    else
      radix_sort_linear_nussinov_starts(candidate_starts,
                                        workspace.radix_scratch);

    // If a later-starting state has at least the same score, an earlier state
    // represents the same structure plus unused leading bases and is
    // dominated for every future suffix query.  Root is retained because it
    // is the reported prefix optimum, but a dominated state is not a lost
    // derivation and therefore needs no pruning certificate.
    double best_suffix_score = -std::numeric_limits<double>::infinity();
    auto& selected = workspace.selected;
    selected.clear();
    selected.reserve(std::min<size_t>(beam_limit, candidate_starts.size()));
    // The heap algorithms put the element that is "largest" under their
    // comparator at the front.  The better-than comparator therefore makes
    // the worst retained candidate the front element, ready for replacement.
    const auto heap_compare = [](const LinearNussinovRankedCell& lhs,
                                 const LinearNussinovRankedCell& rhs) {
      return linear_nussinov_cell_is_better(lhs, rhs);
    };
    for (size_t position = candidate_starts.size(); position-- > 0; ) {
      const uint start = candidate_starts[position];
      const LinearNussinovCell& cell = candidates[start];
      if (cell.start != 0 && cell.score <= best_suffix_score) {
        workspace.dominated_generation[cell.start] = generation;
        continue;
      }
      best_suffix_score = std::max(best_suffix_score, cell.score);

      double priority = cell.score;
      if (cell.start != 0) {
        const LinearNussinovCell* prefix =
            workspace.beam_roots[cell.start-1];
        assert(prefix);
        priority = prefix->score + cell.score;
      }
      const LinearNussinovRankedCell ranked{&cell, priority};
      if (selected.size() < beam_limit) {
        selected.push_back(ranked);
        std::push_heap(selected.begin(), selected.end(), heap_compare);
      } else if (linear_nussinov_cell_is_better(
                     ranked, selected.front())) {
        std::pop_heap(selected.begin(), selected.end(), heap_compare);
        selected.back() = ranked;
        std::push_heap(selected.begin(), selected.end(), heap_compare);
      }
    }

    const LinearNussinovCell* root = &candidates[0];
    const auto has_root = [&]() {
      return std::any_of(selected.begin(), selected.end(),
                         [](const LinearNussinovRankedCell& ranked) {
                           return ranked.cell->start == 0;
                         });
    };
    if (!has_root()) {
      assert(selected.size() == beam_limit);
      const LinearNussinovRankedCell root_ranked{root, root->score};
      std::pop_heap(selected.begin(), selected.end(), heap_compare);
      selected.back() = root_ranked;
      std::push_heap(selected.begin(), selected.end(), heap_compare);
    }

    LinearNussinovBeam& beam = beams[right];
    beam.clear();
    beam.reserve(selected.size());
    for (const auto& ranked : selected)
      beam.push_back(*ranked.cell);
    if (certify) {
      for (const auto& cell : beam)
        workspace.retained_generation[cell.start] = generation;
      for (const uint start : candidate_starts)
        if (workspace.retained_generation[start] != generation &&
            workspace.dominated_generation[start] != generation)
          certificate.record(candidates[start], right, L);
    }
    std::sort(beam.begin(), beam.end(),
              [](const LinearNussinovCell& lhs,
                 const LinearNussinovCell& rhs) {
                return lhs.start < rhs.start;
              });
    const auto final_root = std::find_if(
        beam.begin(), beam.end(),
        [](const LinearNussinovCell& cell) { return cell.start == 0; });
    assert(final_root != beam.end());
    workspace.beam_roots[right] = &*final_root;
  }

  auto& pending = workspace.pending;
  pending.clear();
  pending.push_back({0, L-1});
  while (!pending.empty()) {
    const auto [start, right] = pending.back();
    pending.pop_back();
    if (start > right)
      continue;
    const LinearNussinovCell* cell =
        find_linear_nussinov_cell(beams[right], start);
    assert(cell);
    if (!cell->paired) {
      if (start < right)
        pending.push_back({start, right-1});
      continue;
    }

    const uint left = cell->pair_left;
    ss[left] = right;
    if (start < left)
      pending.push_back({start, left-1});
    if (cell->inside_start != -1u)
      pending.push_back({cell->inside_start, right-1});
  }

  const LinearNussinovCell* root = workspace.beam_roots[L-1];
  assert(root);
  LinearNussinovResult result;
  result.score = static_cast<float>(root->score);
  if (!certify) {
    result.upper_bound = result.score;
    result.pruned_upper_bound = result.score;
    result.additive_upper_bound = result.score;
    return result;
  }

  const double pruned_upper_bound = std::max(
      root->score, certificate.pruned_upper_bound);
  const double upper_bound = std::min(
      pruned_upper_bound, certificate.additive_upper_bound);
  result.pruned_upper_bound =
      std::max(result.score,
               RelaxedBounds::round_up_to_float(pruned_upper_bound));
  result.additive_upper_bound =
      std::max(result.score, RelaxedBounds::round_up_to_float(
                                 certificate.additive_upper_bound));
  result.upper_bound = std::max(
      result.score, RelaxedBounds::round_up_to_float(upper_bound));
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
                       return w * (p[i][j] - th_) - q[i][j];
                     }, ss, false).score;
}

float
LinearNussinov::
decode(float w, const VVF& p, const SparseFloatMatrix& q, VU& ss)
{
  const uint length = p.size();
  return decode_impl(length, dense_pair_support(length),
                     [&](uint i, uint j) {
                       return w * (p[i][j] - th_) - q.get(i, j);
                     }, ss, false).score;
}

float
LinearNussinov::
decode(float w, const SparseFloatMatrix& p, const VVF& q, VU& ss)
{
  const uint length = p.rows();
  return decode_impl(length, dense_pair_support(length),
                     [&](uint i, uint j) {
                       return w * (p.get(i, j) - th_) - q[i][j];
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
                       return w * (p.get(i, j) - th_) - q.get(i, j);
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
                       return w * (p.get(i, j) - th_) - q[i][j];
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
                       return w * (p.get(i, j) - th_) - q.get(i, j);
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
                       return w * (p.get(i, j) - th_) - q[i][j];
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
                       return w * (p.get(i, j) - th_) - q.get(i, j);
                     }, ss, true);
}

float
LinearNussinov::
decode(const VVF& p, VU& ss, std::string& str)
{
  const uint length = p.size();
  const float score = decode_impl(
      length, dense_pair_support(length),
      [&](uint i, uint j) { return p[i][j] - th_; }, ss, false).score;
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
      [&](uint i, uint j) { return p.get(i, j) - th_; }, ss, false).score;
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
        return p.get(i, j) - th_ + pair_bonus.get(i, j);
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
