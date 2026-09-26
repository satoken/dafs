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

#include "ipknot.h"
#include <cassert>
#include <string>
#include <vector>
#include <algorithm>
#include "ip.h"

static
void
decompose_plevel(const VU& ss, VU& plevel);

static
std::string
make_brackets(const VU& ss, const VU& plevel);

IPknot::
IPknot(const VF& th, int n_th /*=1*/)
  : th_(th),
    alpha_(th_.size(), 1.0),
    levelwise_(true),
    stacking_constraints_(true),
    n_th_(n_th)
{
}

float
IPknot::
decode(float w, const VVF& p, const VVF& q, VU& ss)
{
  IP ip(IP::MAX, n_th_);
  make_objective(ip, w, p, q);
  make_constraints(ip);
  return solve(ip, ss);
}

float
IPknot::
decode(float w, const SparseFloatMatrix& p, const VVF& q, VU& ss)
{
  IP ip(IP::MAX, n_th_);
  make_objective(ip, w, p, q);
  make_constraints(ip);
  return solve(ip, ss);
}

float
IPknot::
decode(float w, const SparseFloatMatrix& p,
       const SparseFloatMatrix& q, VU& ss)
{
  IP ip(IP::MAX, n_th_);
  make_objective(ip, w, p, q);
  make_constraints(ip);
  return solve(ip, ss);
}

float
IPknot::
decode(const VVF& p, VU& ss, std::string& str)
{
  IP ip(IP::MAX, n_th_);
  make_objective(ip, p);
  make_constraints(ip);
  float s=solve(ip, ss);
  str=::make_brackets(ss, plevel_);
  return s;
}

float
IPknot::
decode(const SparseFloatMatrix& p, VU& ss, std::string& str)
{
  IP ip(IP::MAX, n_th_);
  make_objective(ip, p);
  make_constraints(ip);
  float s = solve(ip, ss);
  str = ::make_brackets(ss, plevel_);
  return s;
}

void
IPknot::
make_brackets(const VU& ss, std::string& str) const
{
  VU plevel;
  decompose_plevel(ss, plevel);
  str=::make_brackets(ss, plevel);
}

void
IPknot::
reset_candidates(uint length)
{
  v_left_.assign(th_.size(), SparseVariables(length));
  v_right_.assign(th_.size(), SparseVariables(length));
}

void
IPknot::
add_candidate(IP& ip, uint lv, uint i, uint j, float score)
{
  if (score <= 0.0f) return;
  const int var = ip.make_variable(score * alpha_[lv]);
  v_left_[lv][i].emplace_back(j, var);
  v_right_[lv][j].emplace_back(i, var);
}

void
IPknot::
make_objective(IP& ip, float w, const VVF& p, const VVF& q)
{
  const uint L = p.size();
  reset_candidates(L);
  for (uint j = 1; j < L; ++j)
    for (uint i = 0; i < j; ++i)
      for (uint lv = 0; lv < th_.size(); ++lv)
        add_candidate(ip, lv, i, j, w * (p[i][j] - th_[lv]) - q[i][j]);
  ip.update();
}

void
IPknot::
make_objective(IP& ip, float w, const SparseFloatMatrix& p,
               const VVF& q)
{
  const uint L = p.rows();
  reset_candidates(L);
  // Dense multipliers can make a pair absent from the BPP support profitable.
  for (uint j = 1; j < L; ++j)
    for (uint i = 0; i < j; ++i)
      for (uint lv = 0; lv < th_.size(); ++lv)
        add_candidate(ip, lv, i, j,
                      w * (p.get(i, j) - th_[lv]) - q[i][j]);
  ip.update();
}

void
IPknot::
make_objective(IP& ip, float w, const SparseFloatMatrix& p,
               const SparseFloatMatrix& q)
{
  const uint L = p.rows();
  reset_candidates(L);
  // Visit the union of probability and multiplier supports.  A negative
  // multiplier can introduce a pair absent from LinearPartition's output.
  for (uint i = 0; i < L; ++i) {
    for (const auto [j, probability] : p.ordered_row(i))
      if (i < j)
        for (uint lv = 0; lv < th_.size(); ++lv)
          add_candidate(ip, lv, i, j,
                        w * (probability - th_[lv]) - q.get(i, j));
    for (const auto [j, multiplier] : q.ordered_row(i))
      if (i < j && p.get(i, j) == 0.0f)
        for (uint lv = 0; lv < th_.size(); ++lv)
          add_candidate(ip, lv, i, j, -w * th_[lv] - multiplier);
  }
  ip.update();
}

void
IPknot::
make_objective(IP& ip, const VVF& p)
{
  const uint L = p.size();
  reset_candidates(L);
  for (uint j = 1; j < L; ++j)
    for (uint i = 0; i < j; ++i)
      for (uint lv = 0; lv < th_.size(); ++lv)
        add_candidate(ip, lv, i, j, p[i][j] - th_[lv]);
  ip.update();
}

void
IPknot::
make_objective(IP& ip, const SparseFloatMatrix& p)
{
  const uint L = p.rows();
  reset_candidates(L);
  for (uint i = 0; i < L; ++i)
    for (const auto [j, probability] : p.ordered_row(i))
      if (i < j)
        for (uint lv = 0; lv < th_.size(); ++lv)
          add_candidate(ip, lv, i, j, probability - th_[lv]);
  ip.update();
}

void
IPknot::
make_constraints(IP& ip)
{
  const uint L=v_left_[0].size();
  const uint P=th_.size();

  // constraint 1: each s_i is paired with at most one base
  for (uint i=0; i!=L; ++i)
  {
    int row = ip.make_constraint(IP::UP, 0, 1);
    for (uint lv=0; lv!=P; ++lv)
    {
      for (const auto [j, var] : v_right_[lv][i])
        ip.add_constraint(row, var, 1);
      for (const auto [j, var] : v_left_[lv][i])
        ip.add_constraint(row, var, 1);
    }
  }

  if (levelwise_)
  {
    // constraint 2: disallow pseudoknots in x[lv]
    for (uint lv=0; lv!=P; ++lv)
      for (uint i=0; i<L; ++i)
        for (const auto [j, first] : v_left_[lv][i])
        {
          for (uint k=i+1; k<j; ++k)
            for (const auto [l, second] : v_left_[lv][k])
            {
              if (j<l)
              {
                int row = ip.make_constraint(IP::UP, 0, 1);
                ip.add_constraint(row, first, 1);
                ip.add_constraint(row, second, 1);
              }
            }
        }

    // constraint 3: any x[t]_kl must be pseudoknotted with x[u]_ij for t>u
    for (uint lv=1; lv!=P; ++lv)
      for (uint k=0; k<L; ++k)
        for (const auto [l, var] : v_left_[lv][k])
        {
          for (uint plv=0; plv!=lv; ++plv)
          {
            int row = ip.make_constraint(IP::LO, 0, 0);
            ip.add_constraint(row, var, -1);
            for (uint i=0; i<k; ++i)
              for (const auto [j, crossing] : v_left_[plv][i])
              {
                if (k<j && j<l)
                  ip.add_constraint(row, crossing, 1);
              }
            for (uint i=k+1; i<l; ++i)
              for (const auto [j, crossing] : v_left_[plv][i])
              {
                if (l<j)
                  ip.add_constraint(row, crossing, 1);
              }
          }
        }
  }

  if (stacking_constraints_)
  {
    for (uint lv=0; lv!=P; ++lv)
    {
      // upstream
      for (uint i=0; i<L; ++i)
      {
        int row = ip.make_constraint(IP::LO, 0, 0);
        for (const auto [j, var] : v_right_[lv][i])
          ip.add_constraint(row, var, -1);
        if (i>0)
          for (const auto [j, var] : v_right_[lv][i-1])
            ip.add_constraint(row, var, 1);
        if (i+1<L)
          for (const auto [j, var] : v_right_[lv][i+1])
            ip.add_constraint(row, var, 1);
      }

      // downstream
      for (uint i=0; i<L; ++i)
      {
        int row = ip.make_constraint(IP::LO, 0, 0);
        for (const auto [j, var] : v_left_[lv][i])
          ip.add_constraint(row, var, -1);
        if (i>0)
          for (const auto [j, var] : v_left_[lv][i-1])
            ip.add_constraint(row, var, 1);
        if (i+1<L)
          for (const auto [j, var] : v_left_[lv][i+1])
            ip.add_constraint(row, var, 1);
      }
    }
  }
}

float
IPknot::
solve(IP& ip, VU& ss)
{
  const uint L=v_left_[0].size();
  const uint P=th_.size();

  // execute optimization
  float s=ip.solve();

  // build the result
  ss.resize(L);
  std::fill(ss.begin(), ss.end(), -1u);
  plevel_.resize(L);
  std::fill(plevel_.begin(), plevel_.end(), -1u);
  for (uint lv=0; lv!=P; ++lv)
  {
    for (uint i=0; i<L; ++i)
      for (const auto [j, var] : v_left_[lv][i])
        if (ip.get_value(var)>0.5)
        {
          ss[i]=j; //ss[j]=i;
          plevel_[i]=plevel_[j]=lv;
        }
  }

  if (!levelwise_) decompose_plevel(ss, plevel_);

  return s;
}

LinearIPknot::
LinearIPknot(const VF& th, uint beam_size)
  : th_(th), beam_size_(std::max(1u, beam_size))
{
}

void
LinearIPknot::
dense_to_sparse(const VVF& dense, SparseFloatMatrix& sparse)
{
  const uint length = dense.size();
  sparse.assign(length, length);
  for (uint i = 0; i < length; ++i) {
    assert(dense[i].size() == length);
    for (uint j = 0; j < length; ++j)
      if (dense[i][j] != 0.0f)
        sparse.set(i, j, dense[i][j]);
  }
}

float
LinearIPknot::
decode(float w, const VVF& p, const VVF& q, VU& ss)
{
  SparseFloatMatrix sparse_p, sparse_q;
  dense_to_sparse(p, sparse_p);
  dense_to_sparse(q, sparse_q);
  return decode_sparse(w, sparse_p, sparse_q, ss);
}

float
LinearIPknot::
decode(float w, const SparseFloatMatrix& p, const VVF& q, VU& ss)
{
  SparseFloatMatrix sparse_q;
  dense_to_sparse(q, sparse_q);
  return decode_sparse(w, p, sparse_q, ss);
}

float
LinearIPknot::
decode(float w, const SparseFloatMatrix& p,
       const SparseFloatMatrix& q, VU& ss)
{
  return decode_sparse(w, p, q, ss);
}

float
LinearIPknot::
decode(const VVF& p, VU& ss, std::string& str)
{
  SparseFloatMatrix sparse_p;
  dense_to_sparse(p, sparse_p);
  SparseFloatMatrix empty_q;
  empty_q.assign(p.size(), p.size());
  const float score = decode_sparse(1.0f, sparse_p, empty_q, ss);
  str = brackets_for(ss);
  return score;
}

float
LinearIPknot::
decode(const SparseFloatMatrix& p, VU& ss, std::string& str)
{
  SparseFloatMatrix empty_q;
  empty_q.assign(p.rows(), p.columns());
  const float score = decode_sparse(1.0f, p, empty_q, ss);
  str = brackets_for(ss);
  return score;
}

float
LinearIPknot::
decode_sparse(float w, const SparseFloatMatrix& p,
             const SparseFloatMatrix& q, VU& ss)
{
  const uint length = p.rows();
  assert(p.columns() == length && q.rows() == length &&
         q.columns() == length && !th_.empty());
  ss.assign(length, -1u);
  plevel_.assign(length, -1u);

  std::vector<std::pair<uint, uint>> selected_pairs;
  std::vector<char> used(length, 0);
  float total_score = 0.0f;
  const uint levels = std::min<uint>(th_.size(), n_support_brackets);

  // The support is visited once per level.  With a fixed LinearPartition beam
  // and fixed level count this is linear in the retained sparse support.
  for (uint level = 0; level < levels; ++level) {
    SparseFloatMatrix level_p, level_q;
    level_p.assign(length, length);
    level_q.assign(length, length);

    // Build O(1) crossing witnesses for the sparse candidate scan.  A pair
    // (i,j) crosses a previously selected pair (k,l) iff either
    //   k < i < l < j, or i < k < j < l.
    // The first case needs an endpoint l in (i,j) among intervals that start
    // before i.  A prefix maximum is insufficient: a larger endpoint can
    // hide a smaller valid endpoint.  Every selected level is noncrossing,
    // so a stack per level supplies both the smallest active right endpoint
    // and the largest active start.  Each selected pair is pushed and popped
    // once, giving a strict O(levels * L + selected_pairs) sweep.
    std::vector<uint> min_right_before(length, -1u);
    std::vector<uint> max_active_start(length, 0);
    if (level > 0 && length != 0) {
      std::vector<uint> selected_right(length, -1u);
      std::vector<uint> selected_level(length, -1u);
      for (const auto& [left, right] : selected_pairs)
        if (left < length && right < length && plevel_[left] < levels) {
          selected_right[left] = right;
          selected_level[left] = plevel_[left];
        }
      std::vector<std::vector<std::pair<uint, uint>>> active(levels);
      for (uint position = 0; position < length; ++position) {
        uint minimum_right = length;
        uint maximum_start = 0;
        bool found_active = false;
        for (uint selected_level_index = 0;
             selected_level_index < levels; ++selected_level_index) {
          auto& stack = active[selected_level_index];
          while (!stack.empty() && stack.back().second <= position)
            stack.pop_back();
          if (!stack.empty()) {
            found_active = true;
            minimum_right = std::min(minimum_right, stack.back().second);
            maximum_start = std::max(maximum_start,
                                     stack.back().first);
          }
        }
        min_right_before[position] = found_active ? minimum_right : -1u;
        max_active_start[position] = found_active ? maximum_start : 0;

        const uint right = selected_right[position];
        const uint selected_level_index = selected_level[position];
        if (right > position && right < length &&
            selected_level_index < levels)
          active[selected_level_index].emplace_back(position, right);
      }
    }
    const auto consider = [&](uint i, uint j) {
      if (i >= j || used[i] || used[j])
        return;
      if (level > 0) {
        const bool left_crossing = min_right_before[i] < j;
        const bool right_crossing = max_active_start[j] > i;
        if (!left_crossing && !right_crossing)
          return;
      }
      const float probability = p.get(i, j);
      const float multiplier = q.get(i, j);
      const float score = w * (probability - th_[level]) - multiplier;
      if (!(score > 0.0f))
        return;
      if (probability != 0.0f)
        level_p.set(i, j, probability);
      if (multiplier != 0.0f)
        level_q.set(i, j, multiplier);
    };
    // These scans only need coverage.  ordered_row() sorts every raw row,
    // which is unnecessary here and can cost O(m log m) for a large row.  The
    // q scan includes q-only entries because a negative multiplier can make a
    // pair absent from the posterior support profitable.
    p.for_each_nonzero([&](uint i, uint j, float value) {
      if (value != 0.0f)
        consider(i, j);
    });
    q.for_each_nonzero([&](uint i, uint j, float value) {
      if (value != 0.0f && p.get(i, j) == 0.0f)
        consider(i, j);
    });
    LinearNussinov decoder(th_[level], beam_size_);
    VU level_structure;
    const float level_score = decoder.decode(w, level_p, level_q,
                                             level_structure);
    bool selected_any = false;
    for (uint i = 0; i < length; ++i) {
      const uint j = level_structure[i];
      if (j == -1u || j <= i || used[i] || used[j])
        continue;
      selected_any = true;
      ss[i] = j;
      used[i] = used[j] = 1;
      plevel_[i] = plevel_[j] = level;
      selected_pairs.emplace_back(i, j);
    }
    if (!selected_any)
      continue;
    total_score += level_score;
  }
  last_structure_ = ss;
  return total_score;
}

void
LinearIPknot::
linear_decompose_plevel(const VU& ss, VU& plevel)
{
  // Process pairs by increasing left endpoint.  For a pair (i,j), all
  // already processed pairs start before i, so the only possible crossing is
  // k < i < l < j.  A noncrossing level has nested active arcs, making the
  // top of a stack its smallest active right endpoint.  With the fixed
  // bracket alphabet this is O(n * n_support_brackets) time and O(n) memory.
  const uint length = ss.size();
  const uint level_count = Fold::Decoder::n_support_brackets;
  plevel.assign(length, -1u);
  if (length == 0 || level_count == 0)
    return;

  std::vector<std::vector<std::pair<uint, uint>>> active(level_count);
  for (uint left = 0; left < length; ++left) {
    const uint right = ss[left];
    if (right == -1u || right <= left || right >= length)
      continue;
    uint level = level_count;
    for (uint candidate = 0; candidate < level_count; ++candidate) {
      auto& stack = active[candidate];
      while (!stack.empty() && stack.back().second <= left)
        stack.pop_back();
      if (stack.empty() || stack.back().second >= right) {
        level = candidate;
        break;
      }
    }
    if (level == level_count)
      continue;
    plevel[left] = plevel[right] = level;
    active[level].emplace_back(left, right);
  }
}

std::string
LinearIPknot::
brackets_for(const VU& ss) const
{
  VU levels;
  if (last_structure_ == ss && plevel_.size() == ss.size())
    levels = plevel_;
  else
    linear_decompose_plevel(ss, levels);
  return ::make_brackets(ss, levels);
}

void
LinearIPknot::
make_brackets(const VU& ss, std::string& str) const
{
  VU levels;
  if (last_structure_ == ss && plevel_.size() == ss.size())
    levels = plevel_;
  else
    linear_decompose_plevel(ss, levels);
  str = ::make_brackets(ss, levels);
}

struct cmp_by_degree : public std::less<int>
{
  cmp_by_degree(const VVU& g) : g_(g) {}
  bool operator()(int x, int y) const { return g_[y].size()<g_[x].size(); }
  const VVU& g_;
};

struct cmp_by_count : public std::less<int>
{
  cmp_by_count(const VU& count) : count_(count) { }
  bool operator()(int x, int y) const { return count_[y]<count_[x]; }
  const VU& count_;
};

static void
decompose_plevel(const VU& ss, VU& plevel)
{
  // resolve the symbol of brackets by the graph coloring problem
  uint L=ss.size();
    
  // make an adjacent graph, in which pseudoknotted base-pairs are connected.
  VVU g(L);
  for (uint i=0; i!=L; ++i)
  {
    if (ss[i]==-1u || ss[i]<=i) continue;
    uint j=ss[i];
    for (uint k=i+1; k!=L; ++k)
    {
      if (ss[k]==-1u || ss[k]<=k) continue;
      uint l=ss[k];
      if (k<j && j<l)
      {
        g[i].push_back(k);
        g[k].push_back(i);
      }
    }
  }
  // vertices are indexed by the position of the left base
  VU v;
  for (uint i=0; i!=ss.size(); ++i)
    if (ss[i]!=-1u && i<ss[i]) 
      v.push_back(i);
  // sort vertices by degree
  std::sort(v.begin(), v.end(), cmp_by_degree(g));

  // determine colors
  VU c(L, -1u);
  uint max_color=0;
  for (uint i=0; i!=v.size(); ++i)
  {
    // find the smallest color that is unused
    VU used;
    for (uint j=0; j!=g[v[i]].size(); ++j)
      if (c[g[v[i]][j]]!=-1u) used.push_back(c[g[v[i]][j]]);
    std::sort(used.begin(), used.end());
    used.erase(std::unique(used.begin(), used.end()), used.end());
    uint j=0;
    for (j=0; j!=used.size(); ++j)
      if (used[j]!=j) break;
    c[v[i]]=j;
    max_color=std::max(max_color, j);
  }

  // renumber colors in decentant order by the number of base-pairs for each color
  VU count(max_color+1, 0);
  for (uint i=0; i!=c.size(); ++i)
    if (c[i]!=-1u) count[c[i]]++;
  VU idx(count.size());
  for (uint i=0; i!=idx.size(); ++i) idx[i]=i;
  sort(idx.begin(), idx.end(), cmp_by_count(count));
  VU rev(idx.size());
  for (uint i=0; i!=rev.size(); ++i) rev[idx[i]]=i;
  plevel.resize(L);
  for (uint i=0; i!=c.size(); ++i)
    plevel[i]= c[i]!=-1u ? rev[c[i]] : -1u;
}

static
std::string
make_brackets(const VU& ss, const VU& plevel) 
{
  std::string r(ss.size(), '.');
  for (uint i=0; i!=ss.size(); ++i)
  {
    if (ss[i]!=-1u && i<ss[i])
    {
      uint j=ss[i];
      assert(plevel[i]!=-1u);
      if (plevel[i]<Fold::Decoder::n_support_brackets)
      {
        r[i]=Fold::Decoder::left_brackets[plevel[i]];
        r[j]=Fold::Decoder::right_brackets[plevel[i]];
      }
    }
  }
  return r;
}
