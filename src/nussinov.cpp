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
#include <algorithm>
#include <cassert>
#include <limits>
#include <numeric>
#include <stack>
#include <unordered_map>
#include <unordered_set>

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

const LinearNussinovCell*
find_linear_nussinov_suffix_cell(
    const LinearNussinovBeam& beam, uint minimum_start)
{
  const LinearNussinovCell* best = nullptr;
  for (const auto& cell : beam) {
    if (cell.start < minimum_start)
      continue;
    if (!best || cell.score > best->score ||
        (cell.score == best->score && cell.start > best->start))
      best = &cell;
  }
  return best;
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
  for (uint left = 0; left < length; ++left)
    for (const auto [right, value] : matrix.ordered_row(left))
      if (value != 0.0f && left + 2 < right && right < length)
        pairs_by_right[right].push_back(left);
}

void
deduplicate_pair_support(std::vector<std::vector<uint>>& pairs_by_right)
{
  for (auto& lefts : pairs_by_right) {
    std::sort(lefts.begin(), lefts.end());
    lefts.erase(std::unique(lefts.begin(), lefts.end()), lefts.end());
  }
}

struct LinearNussinovCertificate
{
  std::vector<std::vector<double>> prefix_potentials;
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

std::vector<double>
tightened_vertex_cover(
    const std::vector<std::vector<std::pair<uint, double>>>& adjacency,
    const std::vector<double>& initial, bool reverse_first)
{
  std::vector<double> potential = initial;
  std::vector<double> best = initial;
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
  return best;
}

template <typename PairScore>
LinearNussinovCertificate
make_linear_nussinov_certificate(
    uint length, const std::vector<std::vector<uint>>& pairs_by_right,
    PairScore pair_score)
{
  std::vector<double> left(length, 0.0);
  std::vector<double> right(length, 0.0);
  std::vector<double> incident(length, 0.0);
  std::vector<std::vector<std::pair<uint, double>>> adjacency(length);

  for (uint j = 0; j < length; ++j) {
    for (const uint i : pairs_by_right[j]) {
      if (i >= length || j <= i + 2)
        continue;
      const double value = std::max(0.0,
          static_cast<double>(pair_score(i, j)));
      if (!(value > 0.0))
        continue;
      left[i] = std::max(left[i], value);
      right[j] = std::max(right[j], value);
      incident[i] = std::max(incident[i], value);
      incident[j] = std::max(incident[j], value);
      adjacency[i].push_back({j, value});
      adjacency[j].push_back({i, value});
    }
  }

  std::vector<double> half_incident(length, 0.0);
  for (uint i = 0; i < length; ++i)
    half_incident[i] = 0.5 * incident[i];
  std::vector<std::vector<double>> covers;
  covers.push_back(std::move(left));
  covers.push_back(std::move(right));
  covers.push_back(half_incident);
  covers.push_back(
      tightened_vertex_cover(adjacency, half_incident, false));
  covers.push_back(
      tightened_vertex_cover(adjacency, half_incident, true));

#ifndef NDEBUG
  for (const auto& cover : covers)
    for (uint i = 0; i < length; ++i)
      for (const auto& [j, value] : adjacency[i])
        assert(cover[i] + cover[j] + 1e-10 >= value);
#endif

  LinearNussinovCertificate certificate;
  certificate.additive_upper_bound =
      std::numeric_limits<double>::infinity();
  for (const auto& cover : covers) {
    std::vector<double> prefix(length+1, 0.0);
    for (uint i = 0; i < length; ++i)
      prefix[i+1] = prefix[i] + cover[i];
    certificate.additive_upper_bound = std::min(
        certificate.additive_upper_bound, prefix.back());
    certificate.prefix_potentials.push_back(std::move(prefix));
  }
  if (length == 0)
    certificate.additive_upper_bound = 0.0;
  return certificate;
}
} // namespace

template <typename PairScore>
LinearNussinovResult
LinearNussinov::
decode_impl(uint L, const std::vector<std::vector<uint>>& pairs_by_right,
            PairScore pair_score, VU& ss, bool certify)
{
  ss.assign(L, -1u);
  if (L == 0)
    return {};

  LinearNussinovCertificate certificate;
  if (certify)
    certificate = make_linear_nussinov_certificate(
        L, pairs_by_right, pair_score);

  // D(i,j) is the best non-crossing matching on [i,j].  At each right
  // endpoint retain only the most promising interval starts.  With a fixed
  // beam and threshold-sparse pair support this is linear in sequence length.
  std::vector<LinearNussinovBeam> beams(L);
  const uint beam_limit = std::max(1u, beam_size_);

  for (uint right = 0; right < L; ++right) {
    std::unordered_map<uint, LinearNussinovCell> candidates;
    candidates.reserve(beam_limit * 2 + pairs_by_right[right].size());

    const auto offer = [&](uint start, double score, bool paired, uint left,
                           uint inside_start) {
      const auto it = candidates.find(start);
      if (it == candidates.end()) {
        candidates.emplace(start, LinearNussinovCell{
            start, score, left, inside_start, paired});
      } else if (score > it->second.score) {
        it->second = LinearNussinovCell{
            start, score, left, inside_start, paired};
      }
    };

    // Leave the new right endpoint unpaired, including the new singleton
    // interval.  Offering these first makes ties deterministic.
    offer(right, 0.0, false, -1u, -1u);
    if (right > 0)
      for (const auto& previous : beams[right-1])
        offer(previous.start, previous.score, false, -1u, -1u);

    for (const uint left : pairs_by_right[right]) {
      const double local_score = pair_score(left, right);
      if (!(local_score > 0.0))
        continue;

      // A structure whose first used base is later than left+1 is a valid
      // inside structure with the skipped leading bases left unpaired.  This
      // suffix dominance removes the O(length) family of equivalent empty or
      // leading-unpaired states without losing an exact derivation.
      const LinearNussinovCell* inside =
          find_linear_nussinov_suffix_cell(beams[right-1], left+1);
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

    LinearNussinovBeam beam;
    beam.reserve(candidates.size());
    for (const auto& entry : candidates)
      beam.push_back(entry.second);

    // If a later-starting state has at least the same score, an earlier state
    // represents the same structure plus unused leading bases and is
    // dominated for every future suffix query.  Root is retained because it
    // is the reported prefix optimum, but a dominated state is not a lost
    // derivation and therefore needs no pruning certificate.
    std::sort(beam.begin(), beam.end(),
              [](const LinearNussinovCell& lhs,
                 const LinearNussinovCell& rhs) {
                return lhs.start > rhs.start;
              });
    std::unordered_set<uint> dominated;
    double best_suffix_score = -std::numeric_limits<double>::infinity();
    LinearNussinovBeam nondominated;
    nondominated.reserve(beam.size());
    for (const auto& cell : beam) {
      if (cell.start != 0 && cell.score <= best_suffix_score) {
        dominated.insert(cell.start);
        continue;
      }
      best_suffix_score = std::max(best_suffix_score, cell.score);
      nondominated.push_back(cell);
    }
    beam.swap(nondominated);

    const auto priority = [&](const LinearNussinovCell& cell) {
      if (cell.start == 0)
        return cell.score;
      const LinearNussinovCell* prefix =
          find_linear_nussinov_cell(beams[cell.start-1], 0);
      assert(prefix);
      return prefix->score + cell.score;
    };
    const auto better = [&](const LinearNussinovCell& lhs,
                            const LinearNussinovCell& rhs) {
      const double lhs_priority = priority(lhs);
      const double rhs_priority = priority(rhs);
      return lhs_priority != rhs_priority ? lhs_priority > rhs_priority
                                          : lhs.start < rhs.start;
    };

    const auto root_it = std::find_if(
        beam.begin(), beam.end(),
        [](const LinearNussinovCell& cell) { return cell.start == 0; });
    assert(root_it != beam.end());
    const LinearNussinovCell root = *root_it;
    std::sort(beam.begin(), beam.end(), better);
    if (beam.size() > beam_limit)
      beam.resize(beam_limit);
    if (std::none_of(beam.begin(), beam.end(),
                     [](const LinearNussinovCell& cell) {
                       return cell.start == 0;
                     }))
      beam.back() = root;
    if (certify) {
      std::unordered_set<uint> retained;
      retained.reserve(beam.size());
      for (const auto& cell : beam)
        retained.insert(cell.start);
      for (const auto& [start, cell] : candidates)
        if (retained.find(start) == retained.end() &&
            dominated.find(start) == dominated.end())
          certificate.record(cell, right, L);
    }
    std::sort(beam.begin(), beam.end(),
              [](const LinearNussinovCell& lhs,
                 const LinearNussinovCell& rhs) {
                return lhs.start < rhs.start;
              });
    beams[right] = std::move(beam);
  }

  std::stack<std::pair<uint, uint>> pending;
  pending.push({0, L-1});
  while (!pending.empty()) {
    const auto [start, right] = pending.top();
    pending.pop();
    if (start > right)
      continue;
    const LinearNussinovCell* cell =
        find_linear_nussinov_cell(beams[right], start);
    assert(cell);
    if (!cell->paired) {
      if (start < right)
        pending.push({start, right-1});
      continue;
    }

    const uint left = cell->pair_left;
    ss[left] = right;
    if (start < left)
      pending.push({start, left-1});
    if (cell->inside_start != -1u)
      pending.push({cell->inside_start, right-1});
  }

  const LinearNussinovCell* root =
      find_linear_nussinov_cell(beams[L-1], 0);
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
      RelaxedBounds::round_up_to_float(pruned_upper_bound);
  result.additive_upper_bound =
      RelaxedBounds::round_up_to_float(certificate.additive_upper_bound);
  result.upper_bound = RelaxedBounds::round_up_to_float(upper_bound);
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
