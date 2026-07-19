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
#include <cstring>
#include <cassert>
#include <cmath>
#include <unistd.h>
#include <vector>
#include <queue>
#include <stack>
#include <algorithm>
#include <iostream>
#include <memory>
#include <system_error>
#include <random>
#include <set>
#include <cstdint>
#include <limits>
#include <array>
//#include <fstream>
#include "fa.h"
#include "fold.h"
#include "nussinov.h"
#include "ipknot.h"
#include "align.h"
#include "linearalign.h"
#include "needleman_wunsch.h"
#include "alifold.h"
#include "ip.h"
#include "typedefs.h"
#include "gradient_manager.h"
#include "dafs.h"
#include "linfold_wrapper.h"

namespace Vienna
{
  extern "C"
  {
#include <ViennaRNA/fold.h>
#include <ViennaRNA/fold_vars.h>
#include <ViennaRNA/part_func.h>
#include <ViennaRNA/alifold.h>
#include <ViennaRNA/aln_util.h>
#include <ViennaRNA/utils.h>
#include <ViennaRNA/PS_dot.h>
  };
};

#include "cxxopts.hpp"
#include "spdlog/spdlog.h"
#include "spdlog/sinks/basic_file_sink.h"
#include "spdlog/stopwatch.h"
#define CUTOFF 0.01

DAFS::DAFS()
    : use_dynamic_cbp_(false),
      use_sparse_structure_lagrangian_(false),
      use_sparse_alignment_lagrangian_(false),
      use_linear_structure_decoder_(false),
      use_linear_alignment_decoder_(false),
      use_alifold_(false),
      use_alifold1_(true),
      use_linear_profile_folding_(false),
      g_(42)
{
}

DAFS::~DAFS()
{
}

float DAFS::solve(VU &x, VU &y, VU &z,
                  const SparseFloatMatrix &p_x,
                  const SparseFloatMatrix &p_y,
                  const SparseFloatMatrix &p_z,
                  const ALN &aln1, const ALN &aln2)
{
#if defined(WITH_GLPK) || defined(WITH_CPLEX) || defined(WITH_GUROBI)
  return t_max_ != 0 ? solve_by_dd(x, y, z, p_x, p_y, p_z, aln1, aln2) : solve_by_ip(x, y, z, p_x, p_y, p_z, aln1, aln2);
#else
  return solve_by_dd(x, y, z, p_x, p_y, p_z, aln1, aln2);
#endif
}

static void
transpose_mp(const MP &mp, MP &mp_trans, uint x, uint y)
{
  assert(mp.size() == x);
  mp_trans.resize(y);
  for (uint i = 0; i != mp.size(); ++i)
    for (const auto [j, p] : mp[i])
    {
      mp_trans[j].push_back(std::make_pair(i, p));
    }
  for (uint j = 0; j != mp_trans.size(); ++j)
    sort(mp_trans[j].begin(), mp_trans[j].end());
}

#if 0
static void
print_mp(std::ostream &os, const MP &mp)
{
  for (uint i = 0; i != mp.size(); ++i)
  {
    os << i << ":";
    for (const auto& v : mp[i])
    os << " " << v.first << ":" << v.second;
    os << std::endl;
  }
}

static void
print_bp(std::ostream &os, const BP &bp)
{
  for (uint i = 0; i != bp.size(); ++i)
  {
    os << i << ":";
    for (const auto& v : bp[i])
    os << " " << v.first << ":" << v.second;
    os << std::endl;
  }
}

static void
print_matching_probability(std::ostream &os, const VVF &p)
{
  for (uint i = 0; i != p.size(); ++i)
  {
    std::cout << i << ":";
    for (uint k = 0; k != p[i].size(); ++k)
      if (p[i][k] > CUTOFF)
        os << " " << k << ":" << p[i][k];
    os << std::endl;
  }
}

static void
print_basepairing_probability(std::ostream &os, const VVF &p)
{
  for (uint i = 0; i != p.size(); ++i)
  {
    std::cout << i << ":";
    for (uint j = i + 1; j != p.size(); ++j)
      if (p[i][j] > CUTOFF)
        os << " " << j << ":" << p[i][j];
    os << std::endl;
  }
}
#endif

#if 0
static void
save_bp(std::ostream &os, const std::vector<BP> &bp)
{
  for (uint x = 0; x != bp.size(); ++x)
  {
    os << "> " << x << std::endl;
    for (uint i = 0; i != bp[x].size(); ++i)
    {
      os << i;
      for (const auto& j : bp[x][i])
      os << " " << j.first << ":" << j.second;
      os << std::endl;
    }
  }
}

static void
save_mp(std::ostream &os, const std::vector<std::vector<MP>> &mp)
{
  for (uint x = 0; x != mp.size() - 1; ++x)
  {
    for (uint y = x + 1; y != mp[x].size(); ++y)
    {
      os << "> " << x << " " << y << std::endl;
      for (uint i = 0; i != mp[x][y].size(); ++i)
      {
        os << i;
        for (const auto& k : mp[x][y][i])
        os << " " << k.first << ":" << k.second;
        os << std::endl;
      }
    }
  }
}
#endif

void DAFS::
    relax_matching_probability()
{
  const uint N = fa_.size();
  std::vector<std::vector<MP>> mp(N, std::vector<MP>(N));
  assert(mp_.size() == N);
  assert(mp_[0].size() == N);
  for (uint x = 0; x != N - 1; ++x)
  {
    const uint L1 = fa_[x].size();
    for (uint y = x + 1; y != N; ++y)
    {
      const uint L2 = fa_[y].size();
      SparseFloatMatrix posterior;
      posterior.assign(L1, L2);
      assert(L1 == mp_[x][y].size());

      float sum_w = 0.0;
      for (uint z = 0; z != N; ++z)
      {
        const uint L3 = fa_[z].size();
        assert(L3 == mp_[z][x].size());
        assert(L3 == mp_[z][y].size());
        float w = sim_[z][x] * sim_[z][y];
        ;
        if (w_pct_a_ < 0.0)
          w *= 1.0 / N;
        else if (z == x || z == y)
          w *= (1.0 - w_pct_a_) / 2;
        else
          w *= w_pct_a_ / (N - 2);
        sum_w += w;
        for (uint k = 0; k != L3; ++k)
        {
          for (const auto [i, p_ik] : mp_[z][x][k])
          {
            for (const auto [j, p_jk] : mp_[z][y][k])
            {
              assert(i < L1);
              assert(j < L2);
              posterior.add(i, j, p_ik * p_jk * w);
            }
          }
        }
      }

      mp[x][y].resize(L1);
      for (uint i = 0; i != L1; ++i)
        for (const auto [j, value] : posterior.ordered_row(i)) {
          const float v = value / sum_w;
          if (v > CUTOFF)
            mp[x][y][i].push_back(std::make_pair(j, v));
        }
      transpose_mp(mp[x][y], mp[y][x], L1, L2);
    }
  }

  for (uint x = 0; x != N; ++x)
  {
    mp[x][x].resize(fa_[x].size());
    for (uint i = 0; i != fa_[x].size(); ++i)
      mp[x][x][i].push_back(std::make_pair(i, 1.0f));
  }
  std::swap(mp_, mp);
}

void DAFS::
    relax_basepairing_probability()
{
  const uint N = bp_.size();
  std::vector<BP> bp(N);
  for (uint x = 0; x != N; ++x)
  {
    const uint L1 = bp_[x].size();
    SparseFloatMatrix p;
    p.assign(L1, L1);

    float sum_w = 0.0;
    for (uint y = 0; y != N; ++y)
    {
      const uint L2 = bp_[y].size();
      assert(L2 == mp_[y][x].size());
      float w = sim_[y][x];
      if (w_pct_s_ < 0.0)
        w *= 1.0 / N;
      else if (y == x)
        w *= 1.0 - w_pct_s_;
      else
        w *= w_pct_s_ / (N - 1);
      sum_w += w;
      for (uint k = 0; k != L2; ++k)
      {
        for (const auto [l, p_kl] : bp_[y][k])
        {
          for (const auto [i, p_ik] : mp_[y][x][k])
          {
            for (const auto [j, p_jl] : mp_[y][x][l])
            {
              if (i < j)
                p.add(i, j, p_kl * p_ik * p_jl * w);
            }
          }
        }
      }
    }

    bp[x].resize(L1);
    for (uint i = 0; i + 1 < L1; ++i)
      for (const auto [j, value] : p.ordered_row(i))
      {
        const float v = value / sum_w;
        if (v > CUTOFF)
          bp[x][i].push_back(std::make_pair(j, v));
      }
  }
  std::swap(bp_, bp);
}

void DAFS::
    relax_fourway_consistency()
{
  const uint N = fa_.size();
  std::vector<std::vector<MP>> mp(N, std::vector<MP>(N));
  assert(mp_.size() == N);
  assert(mp_[0].size() == N);
  for (uint x = 0; x != N - 1; ++x)
  {
    const uint L1 = fa_[x].size();
    for (uint y = x + 1; y != N; ++y)
    {
      const uint L2 = fa_[y].size();
      SparseFloatMatrix posterior;
      posterior.assign(L1, L2);
      assert(L1 == mp_[x][y].size());

      for (uint i = 0; i != L1; ++i)
      {
        for (const auto [k, p_ik] : mp_[x][y][i])
        {
          posterior.add(i, k, p_ik * (1.0 - w_pct_f_));
          for (const auto [j, p_ij] : bp_[x][i])
          {
            auto ll1 = mp_[x][y][j].begin();
            auto ll2 = bp_[y][k].begin();
            while (ll1 != mp_[x][y][j].end() && ll2 != bp_[y][k].end())
            {
              if (ll1->first < ll2->first)
                ++ll1;
              else if (ll1->first > ll2->first)
                ++ll2;
              else /* if (ll1->first==ll2->first) */
              {
                const uint l = ll1->first;
                const float p_jl = ll1->second;
                const float p_kl = ll2->second;
                posterior.add(i, k, p_ij * p_kl * p_jl * w_pct_f_);
                posterior.add(j, l, p_ij * p_kl * p_ik * w_pct_f_);
                ++ll1;
                ++ll2;
              }
            }
          }
        }
      }

      mp[x][y].resize(L1);
      for (uint i = 0; i != L1; ++i)
        for (const auto [j, v] : posterior.ordered_row(i))
          if (v > CUTOFF)
            mp[x][y][i].push_back(std::make_pair(j, v));
      transpose_mp(mp[x][y], mp[y][x], L1, L2);
    }
  }

  for (uint x = 0; x != N; ++x)
  {
    mp[x][x].resize(fa_[x].size());
    for (uint i = 0; i != fa_[x].size(); ++i)
      mp[x][x][i].push_back(std::make_pair(i, 1.0f));
  }
  std::swap(mp_, mp);
}

void DAFS::
    build_tree()
{
  uint n = fa_.size();
  tree_.resize(2 * n - 1);
  std::fill(tree_.begin(), tree_.end(), std::make_pair(0.0, std::make_pair(-1u, -1u)));

  VVF d(n, VF(n, 0.0));
  VU idx(2 * n - 1, -1u);
  for (uint i = 0; i != n; ++i)
    idx[i] = i;

  std::priority_queue<node_t> pq;
  for (uint i = 0; i != n - 1; ++i)
  {
    for (uint j = i + 1; j != n; ++j)
    {
      d[i][j] = d[j][i] = sim_[i][j];
      pq.push(std::make_pair(sim_[i][j], std::make_pair(i, j)));
    }
  }

  while (!pq.empty())
  {
    node_t t = pq.top();
    pq.pop();
    if (idx[t.second.first] != -1u && idx[t.second.second] != -1u)
    {
      assert(n < tree_.size());
      const uint l = idx[t.second.first];
      const uint r = idx[t.second.second];
      idx[t.second.first] = idx[t.second.second] = -1u;
      for (uint i = 0; i != n; ++i)
      {
        if (idx[i] != -1u)
        {
          uint ii = idx[i];
          d[ii][l] = d[l][ii] = (d[ii][l] + d[ii][r]) * t.first / 2;
          pq.push(std::make_pair(d[ii][l], std::make_pair(i, n)));
        }
      }
      tree_[n] = t;
      idx[n++] = l;
    }
  }
  assert(n == tree_.size());
}

// print the guide tree
void DAFS::
    print_tree(std::ostream &os, int i) const
{
  if (tree_[i].second.first == -1u)
  {
    assert(tree_[i].second.second == -1u);
    os << fa_[i].name();
  }
  else
  {
    os << "[ " << tree_[i].first << " ";
    print_tree(os, tree_[i].second.first);
    os << " ";
    print_tree(os, tree_[i].second.second);
    os << " ]";
  }
}

void DAFS::
    calculate_profile_basepairing_probability(const ALN &aln, BP &bp) const
{
  if (use_linear_profile_folding_)
    s_model_->calculate(profile_consensus_sequence(aln), bp);
  else
    Alifold(0.0 /*CUTOFF*/).fold(aln, fa_, bp);
}

void DAFS::
    calculate_profile_basepairing_probability(
        const ALN &aln, const std::string &constraint, BP &bp) const
{
  if (use_linear_profile_folding_)
    s_model_->calculate(profile_consensus_sequence(aln), constraint, bp);
  else
    Alifold(0.0 /*CUTOFF*/).fold(aln, fa_, constraint, bp);
}

std::string DAFS::
    profile_consensus_sequence(const ALN &aln) const
{
  const uint L = aln.front().second.size();
  std::vector<std::array<uint, 4>> counts(L);
  for (auto& count : counts)
    count.fill(0);

  for (const auto& entry : aln) {
    assert(entry.second.size() == L);
    const std::string& sequence = fa_[entry.first].seq();
    uint sequence_position = 0;
    for (uint column = 0; column < L; ++column) {
      if (!entry.second[column])
        continue;
      assert(sequence_position < sequence.size());
      switch (sequence[sequence_position++]) {
      case 'A': case 'a': ++counts[column][0]; break;
      case 'C': case 'c': ++counts[column][1]; break;
      case 'G': case 'g': ++counts[column][2]; break;
      case 'U': case 'u': case 'T': case 't': ++counts[column][3]; break;
      default: break;
      }
    }
    assert(sequence_position == sequence.size());
  }

  static constexpr char bases[] = {'A', 'C', 'G', 'U'};
  std::string consensus(L, 'A');
  for (uint column = 0; column < L; ++column) {
    uint best = 0;
    for (uint base = 1; base < 4; ++base)
      if (counts[column][base] > counts[column][best])
        best = base;
    consensus[column] = bases[best];
  }
  return consensus;
}

void DAFS::
    average_matching_probability(SparseFloatMatrix &posterior,
                                 const ALN &aln1, const ALN &aln2) const
{
  const uint L1 = aln1.front().second.size();
  const uint L2 = aln2.front().second.size();
  const uint N1 = aln1.size();
  const uint N2 = aln2.size();
  SparseFloatMatrix p;
  p.assign(L1, L2);
  for (const auto& it1 : aln1)
  {
    assert(L1 == it1.second.size());
    VU idx1(fa_[it1.first].size());
    for (uint i = 0, ii = 0; i != L1; ++i)
      if (it1.second[i])
        idx1[ii++] = i;
    for (const auto& it2 : aln2)
    {
      assert(L2 == it2.second.size());
      VU idx2(fa_[it2.first].size());
      for (uint j = 0, jj = 0; j != L2; ++j)
        if (it2.second[j])
          idx2[jj++] = j;
      const MP &m = mp_[it1.first][it2.first];
      for (uint ii = 0; ii != m.size(); ++ii) {
        assert(ii < m.size());
        for (const auto [jj, value] : m[ii])
          p.add(idx1[ii], idx2[jj], value / (N1 * N2));
      }
    }
  }
  p.prune(CUTOFF, 1.0f);
  posterior = std::move(p);
}

void DAFS::
    average_basepairing_probability(SparseFloatMatrix &posterior,
                                    const ALN &aln, bool use_alifold) const
{
  // calculate an averaged base-pairing probabilities
  const uint L = aln.front().second.size();
  const uint N = aln.size();
  SparseFloatMatrix p;
  p.assign(L, L);
  for (const auto& it : aln)
  {
    assert(L == it.second.size());
    uint s = it.first;
    VU idx(fa_[s].size());
    for (uint i = 0, j = 0; i != L; ++i)
      if (it.second[i])
        idx[j++] = i;
    const BP &bp = bp_[s];
    for (uint i = 0; i != bp.size(); ++i)
      for (uint j = 0; j != bp[i].size(); ++j)
        p.add(idx[i], idx[bp[i][j].first], bp[i][j].second / N);
  }

  // mix base-pairing probabilities by alifold and averaged base-pairing probabilities
  if (use_alifold)
  {
    BP bp;
    calculate_profile_basepairing_probability(aln, bp);
    assert(L == bp.size());
    for (uint i = 0; i != bp.size(); ++i)
      for (uint j = 0; j != bp[i].size(); ++j)
        p.add(i, bp[i][j].first, bp[i][j].second);

    p.scale(0.5f);
  }

  p.prune(CUTOFF);
  for (uint i = 0; i < L; ++i)
    for (const auto [j, value] : p.ordered_row(i))
      assert(j > i && value <= 1.0f);
  posterior = std::move(p);
}

void DAFS::
    update_basepairing_probability(SparseFloatMatrix &posterior,
                                   const VU &ss, const std::string &str,
                                   const ALN &aln, bool use_alifold) const
{
  const uint L = aln.front().second.size();
  const uint N = aln.size();
  const uint plevel = th_s_.size();
  SparseFloatMatrix p;
  p.assign(L, L);

  // calculate an averaged base-pairing probabilities
  //   which are constrained by the previous prediction
  for (const auto& it : aln)
  {
    assert(L == it.second.size());
    // calculate the mapping
    uint s = it.first;
    VU idx(fa_[s].size()); // from the sequence to the alignment
    VU rev(L, -1u);        // from the alignment to the sequence
    for (uint i = 0, j = 0; i != L; ++i)
      if (it.second[i])
      {
        idx[j] = i;
        rev[i] = j;
        j++;
      }

    for (uint plv = 0; plv != plevel; ++plv)
    {
      // make the constraint from the prediction
      std::string con(fa_[s].size(), '?');
      for (uint i = 0; i != L; ++i)
      {
        if (ss[i] != -1u && rev[i] != -1u && rev[ss[i]] != -1u)
        {
          if (str[i] == Fold::Decoder::left_brackets[plv])
          {
            con[rev[i]] = '(';
            con[rev[ss[i]]] = ')';
          }
          else
          {
            con[rev[i]] = con[rev[ss[i]]] = '.';
          }
        }
      }

      // calculate base-pairing probabilities under the constraint
      BP bp;
      s_model_->calculate(fa_[s].seq(), con, bp);
      for (uint i = 0; i != bp.size(); ++i)
        for (uint j = 0; j != bp[i].size(); ++j)
          p.add(idx[i], idx[bp[i][j].first], bp[i][j].second / N);
    }
  }

  // mix base-pairing probabilities by alifold and averaged base-pairing probabilities
  if (use_alifold)
  {
    for (uint plv = 0; plv != plevel; ++plv)
    {
      // make the constraint from the prediction
      std::string con(L, '?');
      for (uint i = 0; i != L; ++i)
      {
        if (ss[i] != -1u)
        {
          if (str[i] == Fold::Decoder::left_brackets[plv])
          {
            con[i] = '(';
            con[ss[i]] = ')';
          }
          else
          {
            con[i] = con[ss[i]] = '.';
          }
        }
      }

      // calculate base-pairing probabilities under the constraint
      BP bp;
      calculate_profile_basepairing_probability(aln, con, bp);
      assert(L == bp.size());
      for (uint i = 0; i != bp.size(); ++i)
        for (uint j = 0; j != bp[i].size(); ++j)
          p.add(i, bp[i][j].first, bp[i][j].second);
    }

    p.scale(0.5f);
  }

  p.prune(CUTOFF);
  for (uint i = 0; i < L; ++i)
    for (const auto [j, value] : p.ordered_row(i))
      assert(j > i && value <= 1.0f);
  posterior = std::move(p);
}

static float
calculate_similarity_score(const MP &mp, uint L1, uint L2)
{
  assert(mp.size() == L1);

  VVF dp(L1 + 1, VF(L2 + 1, 0.0));
  VVI tr(L1 + 1, VI(L2 + 1, 0));
  for (uint i = 1; i != L1 + 1; ++i)
  {
    uint j = 1;
    for (const auto& jj : mp[i - 1])
    {
      for (; j - 1 < jj.first; ++j)
      {
        dp[i][j] = dp[i][j - 1];
        tr[i][j] = tr[i][j - 1] + 1;
        if (dp[i][j] < dp[i - 1][j])
        {
          dp[i][j] = dp[i - 1][j];
          tr[i][j] = tr[i - 1][j] + 1;
        }
      }

      dp[i][j] = dp[i - 1][j - 1] + jj.second;
      tr[i][j] = tr[i - 1][j - 1] + 1;
      if (dp[i][j] < dp[i][j - 1])
      {
        dp[i][j] = dp[i][j - 1];
        tr[i][j] = tr[i][j - 1] + 1;
      }
      if (dp[i][j] < dp[i - 1][j])
      {
        dp[i][j] = dp[i - 1][j];
        tr[i][j] = tr[i - 1][j] + 1;
      }
      ++j;
    }

    for (; j < L2 + 1; ++j)
    {
      dp[i][j] = dp[i][j - 1];
      tr[i][j] = tr[i][j - 1] + 1;
      if (dp[i][j] < dp[i - 1][j])
      {
        dp[i][j] = dp[i - 1][j];
        tr[i][j] = tr[i - 1][j] + 1;
      }
    }
  }
  //return dp[L1][L2]/(std::min(L1, L2)); // simple version
  return dp[L1][L2] / tr[L1][L2];
}

void DAFS::
    project_alignment(ALN &aln, const ALN &aln1, const ALN &aln2, const VU &z) const
{
  const uint L1 = aln1[0].second.size();
  const uint L2 = aln2[0].second.size();
  aln.resize(aln1.size() + aln2.size());
  uint c = 0;
  for (uint i = 0; i != z.size(); ++i)
    if (z[i] != -1u)
      c++;
  const uint L = L1 + L2 - c;
  ALN::iterator p = aln.begin();
  for (const auto& q : aln1)
  {
    p->first = q.first;
    p->second.resize(L, false);
    uint r = 0, k = 0;
    for (uint i = 0; i != z.size(); ++i)
    {
      if (z[i] != -1u)
      {
        while (k < z[i])
        {
          p->second[r++] = false;
          k++;
        }
        p->second[r++] = q.second[i];
        ++k;
      }
      else
        p->second[r++] = q.second[i];
    }
    while (k < L2)
    {
      p->second[r++] = false;
      k++;
    }
    ++p;
  }
  for (const auto& q : aln2)
  {
    p->first = q.first;
    p->second.resize(L, false);
    uint k = 0, r = 0;
    for (uint i = 0; i != z.size(); ++i)
    {
      if (z[i] != -1u)
      {
        while (k < z[i])
          p->second[r++] = q.second[k++];
        p->second[r++] = q.second[k++];
      }
      else
        p->second[r++] = false;
    }
    while (k < L2)
      p->second[r++] = q.second[k++];
    ++p;
  }
}

void DAFS::
    project_secondary_structure(VU &xx, VU &yy, const VU &x, const VU &y, const VU &z) const
{
  const uint L1 = x.size();
  const uint L2 = y.size();
  VU idx1(L1, -1u), idx2(L2, -1u);
  uint r = 0, k = 0;
  for (uint i = 0; i != z.size(); ++i)
  {
    if (z[i] != -1u)
    {
      while (k < z[i])
      {
        idx2[k] = r;
        ++r;
        ++k;
      }
      idx1[i] = r;
      idx2[k] = r;
      ++r;
      ++k;
    }
    else
    {
      idx1[i] = r;
      ++r;
    }
  }
  while (k < L2)
  {
    idx2[k] = r;
    ++r;
    ++k;
  }
  const uint L = r;

  xx.resize(L);
  std::fill(xx.begin(), xx.end(), -1u);
  yy.resize(L);
  std::fill(yy.begin(), yy.end(), -1u);
  for (uint i = 0; i != L1; ++i)
    if (x[i] != -1u)
      xx[idx1[i]] = idx1[x[i]];
  for (uint k = 0; k != L2; ++k)
    if (y[k] != -1u)
      yy[idx2[k]] = idx2[y[k]];
}

void DAFS::
    output_verbose(const VU &x, const VU &y, const VU &z, const ALN &aln1, const ALN &aln2) const
{
  ALN aln;
  project_alignment(aln, aln1, aln2, z);
  VU xx, yy;
  project_secondary_structure(xx, yy, x, y, z);

  std::string x_str, y_str;
  s_decoder_->make_brackets(xx, x_str);
  s_decoder_->make_brackets(yy, y_str);

  output(std::cout, aln.begin(), aln.begin() + aln1.size());
  std::cout << x_str << std::endl;

  output(std::cout, aln.begin() + aln1.size(), aln.end());
  std::cout << y_str << std::endl;

  std::cout << std::endl;
}

void DAFS::
    align_alignments(ALN &aln, const ALN &aln1, const ALN &aln2)
{
  // calculate posteriors
  SparseFloatMatrix p_x, p_y, p_z;
  average_basepairing_probability(p_x, aln1, use_alifold_);
  average_basepairing_probability(p_y, aln2, use_alifold_);
  average_matching_probability(p_z, aln1, aln2);

  // solve the problem
  VU x, y, z;
  solve(x, y, z, p_x, p_y, p_z, aln1, aln2);

  // build the result
  project_alignment(aln, aln1, aln2, z);
}

float DAFS::
    align_alignments(VU &ss, ALN &aln, const ALN &aln1, const ALN &aln2)
{
  // calculate posteriors
  SparseFloatMatrix p_x, p_y, p_z;
  average_basepairing_probability(p_x, aln1, use_alifold_);
  if (use_bp_update_)
  {
    std::string str;
    VU ss;
    s_decoder_->decode(p_x, ss, str);
    update_basepairing_probability(p_x, ss, str, aln1, use_alifold_);
  }

  average_basepairing_probability(p_y, aln2, use_alifold_);
  if (use_bp_update_)
  {
    std::string str;
    VU ss;
    s_decoder_->decode(p_y, ss, str);
    update_basepairing_probability(p_y, ss, str, aln2, use_alifold_);
  }

  average_matching_probability(p_z, aln1, aln2);

  // solve the problem
  VU x, y, z;
  float s = solve(x, y, z, p_x, p_y, p_z, aln1, aln2);

  // build the result
  project_alignment(aln, aln1, aln2, z);
  VU xx, yy;
  project_secondary_structure(xx, yy, x, y, z);
  assert(xx.size() == yy.size());
  ss.resize(xx.size());
  std::fill(ss.begin(), ss.end(), -1u);
  for (uint i = 0; i != ss.size(); ++i)
    if (xx[i] == yy[i])
      ss[i] = xx[i];

#if 0 // calculate the score exactly
  // calculate alignment score
  float a_score=0.0;
  for (uint i=0; i!=z.size(); ++i)
    if (z[i]!=-1u) a_score += p_z[i][z[i]]-th_a_;

  // calculate folding score
  assert(x.size()==z.size());
  float s_score=0.0;
  for (uint i=0; i!=x.size(); ++i)
  {
    if (x[i]!=-1u && z[i]!=-1u)
    {
      const uint j=x[i];
      const uint k=z[i];
      if (z[j]!=-1u && z[j]==y[k])
      {
        const uint l=y[k];
        s_score += p_x[i][j]-th_s_;
        s_score += p_y[k][l]-th_s_;
      }
    }
  }

  return w_*s_score + a_score;
#else
  return s; // this is the minimized upper bound of the exact score.
#endif
}

float DAFS::
    calculate_alignment_only_score(VU &ss, ALN &aln, const ALN &aln1, const ALN &aln2) const
{
  // 1. Use alignment decoder only with zero Lagrange multipliers
  SparseFloatMatrix p_z;
  average_matching_probability(p_z, aln1, aln2);
  
  // Initialize alignment decoder
  a_decoder_->initialize(p_z);
  
  VU z;
  float alignment_score = a_decoder_->decode(p_z, z);
  
  // Build the combined alignment from the decoded alignment
  project_alignment(aln, aln1, aln2, z);
  
  // 2. Compute common secondary structure like the final step in run method
  SparseFloatMatrix p;
  average_basepairing_probability(p, aln, false /*use_alifold1_*/);
#if 0 
  if (use_bp_update1_ && s_decoder1_)
  {
    std::string str;
    VU ss_temp;
    s_decoder1_->decode(p, ss_temp, str);
    update_basepairing_probability(p, ss_temp, str, aln, use_alifold1_);
  }
#endif
  std::string str;
  float structure_score = w_ * 2 * s_decoder_->decode(p, ss, str);

  // Return total score (alignment + secondary structure)
#if 0
  spdlog::info("Initial lower bound for adaptive method: {} (={}+{})", alignment_score + structure_score, alignment_score, structure_score);
  if (verbose_ >= 1)
    output(std::cout, aln);
#endif
  return alignment_score + structure_score;
}

float DAFS::
    repair_feasible_solution(VU& repaired_x, VU& repaired_y, VU& repaired_z,
                             const VU& x, const VU& y, const VU& z,
                             const SparseFloatMatrix& p_x,
                             const SparseFloatMatrix& p_y,
                             const SparseFloatMatrix& p_z,
                             uint N1, uint N2, float min_th_s) const
{
  const uint L1 = p_x.rows();
  const uint L2 = p_y.rows();
  repaired_x.assign(L1, -1u);
  repaired_y.assign(L2, -1u);
  repaired_z = z;

  // Every decoded z is a feasible monotone alignment.  Its original (not
  // Lagrangian-shifted) contribution is therefore a valid primal score.
  float score = 0.0f;
  for (uint i = 0; i < L1; ++i) {
    const uint k = z[i];
    if (k != -1u)
      score += p_z.get(i, k) - th_a_;
  }

  const float x_weight = w_ * 2.0f * N1 / (N1 + N2);
  const float y_weight = w_ * 2.0f * N2 / (N1 + N2);
  // With multiple IPknot levels the level assignment is not retained in VU.
  // Using the largest threshold is conservative and keeps this a lower bound.
  const float repair_th = *std::max_element(th_s_.begin(), th_s_.end());

  // Keep only base pairs on which both folding subproblems and both alignment
  // endpoints agree.  Removing pairs preserves the structural constraints,
  // so this intersection is feasible for the original coupled problem.
  for (uint i = 0; i + 1 < L1; ++i) {
    const uint j = x[i];
    if (j == -1u || j <= i || j >= L1)
      continue;
    const uint k = z[i];
    const uint l = z[j];
    if (k == -1u || l == -1u || k >= l || l >= L2 || y[k] != l)
      continue;
    if (!is_valid_cbp(i, j, k, l, p_x, p_y, p_z,
                      N1, N2, min_th_s))
      continue;

    repaired_x[i] = j;
    repaired_y[k] = l;
    score += x_weight * (p_x.get(i, j) - repair_th);
    score += y_weight * (p_y.get(k, l) - repair_th);
  }

  return score;
}

float DAFS::
    solve_by_dd(VU &x, VU &y, VU &z,
                const SparseFloatMatrix &p_x,
                const SparseFloatMatrix &p_y,
                const SparseFloatMatrix &p_z,
                const ALN &aln1, const ALN &aln2)
{
  const uint L1 = p_x.rows();
  const uint L2 = p_y.rows();
  const uint N1 = aln1.size();
  const uint N2 = aln2.size();

  std::vector<CBP> cbp;
  VU w_cbp;
  VVU c_x(L1), c_y(L2), c_z(L1);
  std::vector<uint> cbp_inactive_steps;
  std::vector<std::pair<uint, uint>> q_x_support, q_y_support;
  std::unordered_set<uint64_t> q_x_support_set, q_y_support_set;
  const auto support_key = [](uint i, uint j) {
    return (static_cast<uint64_t>(i) << 32) | static_cast<uint64_t>(j);
  };
  const auto add_support = [&](std::vector<std::pair<uint, uint>>& support,
                               std::unordered_set<uint64_t>& support_set,
                               uint i, uint j) {
    if (support_set.insert(support_key(i, j)).second)
      support.emplace_back(i, j);
  };

  // Sparse alignment support used by the exact pricing oracle.  With a fixed
  // alignment beam, every row has O(beam) entries.
  VVU p_z_forward(L1), p_z_reverse(L2);
  std::vector<std::pair<uint, uint>> p_x_support, p_y_support, p_z_support;
  for (uint i = 0; i < L1; ++i)
    for (const auto [j, value] : p_x.ordered_row(i))
      if (j > i && value != 0.0f)
        p_x_support.emplace_back(i, j);
  for (uint k = 0; k < L2; ++k)
    for (const auto [l, value] : p_y.ordered_row(k))
      if (l > k && value != 0.0f)
        p_y_support.emplace_back(k, l);
  for (uint i = 0; i < L1; ++i) {
    for (const auto [k, value] : p_z.ordered_row(i)) {
      if (value != 0.0f) {
        p_z_support.emplace_back(i, k);
      }
      if (value > CUTOFF) {
        p_z_forward[i].push_back(k);
        p_z_reverse[k].push_back(i);
      }
    }
  }
  
  // Clear CBP set for dynamic generation
  cbp_set_.clear();
  
  float min_th_s = *std::min_element(th_s_.begin(), th_s_.end());

  // precalculate the range for alignment, i.e. alignment envelope
  // (moved before dynamic CBP generation to enable decoder-based initial solution)
  VVF dense_p_x, dense_p_y, dense_p_z;
  if (!use_linear_structure_decoder_) {
    dense_p_x = p_x.dense();
    dense_p_y = p_y.dense();
  }
  if (!use_linear_alignment_decoder_)
    dense_p_z = p_z.dense();
  if (use_linear_alignment_decoder_)
    a_decoder_->initialize(p_z);
  else
    a_decoder_->initialize(dense_p_z);

  if (!use_dynamic_cbp_) {
    // Original static pre-enumeration
    for (const auto& [i, j] : p_x_support)
      if (p_x.get(i, j) > CUTOFF)
        for (const uint k : p_z_forward[i])
          for (const uint l : p_z_forward[j])
            if (k < l && p_y.get(k, l) > CUTOFF)
                {
                  assert(p_x.get(i, j) <= 1.0);
                  assert(p_y.get(k, l) <= 1.0);
                  float p = (N1 * p_x.get(i, j) + N2 * p_y.get(k, l)) / (N1 + N2);
                  float q = (p_z.get(i, k) + p_z.get(j, l)) / 2;
                  if (p - min_th_s > 0.0 && w_ * (p - min_th_s) + (q - th_a_) > 0.0)
                  {
                    cbp.push_back(std::make_pair(std::make_pair(i, j), std::make_pair(k, l)));
                    c_x[i].push_back(j);
                    c_y[k].push_back(l);
                    c_z[i].push_back(k);
                    c_z[j].push_back(l);
                  }
                }
              
    for (uint i = 0; i != c_x.size(); ++i)
    {
      std::sort(c_x[i].begin(), c_x[i].end());
      c_x[i].erase(std::unique(c_x[i].begin(), c_x[i].end()), c_x[i].end());
    }
    for (uint k = 0; k != c_y.size(); ++k)
    {
      std::sort(c_y[k].begin(), c_y[k].end());
      c_y[k].erase(std::unique(c_y[k].begin(), c_y[k].end()), c_y[k].end());
    }
    for (uint i = 0; i != c_z.size(); ++i)
    {
      std::sort(c_z[i].begin(), c_z[i].end());
      c_z[i].erase(std::unique(c_z[i].begin(), c_z[i].end()), c_z[i].end());
    }
  }  // end of else block for static pre-enumeration

  // Initialize gradient manager
  GradientManager gm(eta0_, 0.0f, 1.0f,
                     use_sparse_structure_lagrangian_,
                     use_sparse_alignment_lagrangian_);
  gm.initialize(L1, L2);
  
  // Certified bounds and the best feasible primal solution recovered so far.
  // The former alignment-only bound could contain a consensus pair across a
  // gap in one profile; use only explicitly verified coupled solutions here.
  float lb = std::numeric_limits<float>::lowest();
  float best_ub = std::numeric_limits<float>::infinity();
  float best_feasible_score = std::numeric_limits<float>::lowest();
  VU best_x, best_y, best_z;
  
  float s_prev = 0.0;
  uint violated = 0;
  uint t;
  for (t = 0; t != t_max_; ++t)
  {
    // solve the subproblems
    float s_x = 0.0f, s_y = 0.0f, s_z = 0.0f;
    if (gm.uses_sparse_structure_storage()) {
      if (use_linear_structure_decoder_) {
        s_x = s_decoder_->decode(w_ * 2 * N1 / (N1 + N2),
                                 p_x, gm.sparse_q_x(), x);
        s_y = s_decoder_->decode(w_ * 2 * N2 / (N1 + N2),
                                 p_y, gm.sparse_q_y(), y);
      } else {
        s_x = s_decoder_->decode(w_ * 2 * N1 / (N1 + N2),
                                 dense_p_x, gm.sparse_q_x(), x);
        s_y = s_decoder_->decode(w_ * 2 * N2 / (N1 + N2),
                                 dense_p_y, gm.sparse_q_y(), y);
      }
    } else {
      if (use_linear_structure_decoder_) {
        s_x = s_decoder_->decode(w_ * 2 * N1 / (N1 + N2),
                                 p_x, gm.dense_q_x(), x);
        s_y = s_decoder_->decode(w_ * 2 * N2 / (N1 + N2),
                                 p_y, gm.dense_q_y(), y);
      } else {
        s_x = s_decoder_->decode(w_ * 2 * N1 / (N1 + N2),
                                 dense_p_x, gm.dense_q_x(), x);
        s_y = s_decoder_->decode(w_ * 2 * N2 / (N1 + N2),
                                 dense_p_y, gm.dense_q_y(), y);
      }
    }
    if (gm.uses_sparse_alignment_storage()) {
      s_z = use_linear_alignment_decoder_
          ? a_decoder_->decode(p_z, gm.sparse_q_z(), z)
          : a_decoder_->decode(dense_p_z, gm.sparse_q_z(), z);
    } else {
      s_z = use_linear_alignment_decoder_
          ? a_decoder_->decode(p_z, gm.dense_q_z(), z)
          : a_decoder_->decode(dense_p_z, gm.dense_q_z(), z);
    }
    float s = s_x + s_y + s_z;

    // Beam maxima are feasible subproblem values, not certified maxima.
    // For each Linear component, independently selecting every positive
    // local item is a relaxation and therefore a safe upper bound.
    const float x_weight = w_ * 2 * N1 / (N1 + N2);
    const float y_weight = w_ * 2 * N2 / (N1 + N2);
    const auto relaxed_structure_bound = [](const auto& support,
                                            const SparseFloatMatrix& p, float weight,
                                            float threshold,
                                            const auto& multiplier) {
      float bound = 0.0f;
      for (const auto& [i, j] : support)
        bound += std::max(0.0f, weight * (p.get(i, j) - threshold)
                                - multiplier(i, j));
      return bound;
    };
    const auto relaxed_alignment_bound = [](const auto& support,
                                            const SparseFloatMatrix& p, float threshold,
                                            const auto& multiplier) {
      float bound = 0.0f;
      for (const auto& [i, k] : support)
        bound += std::max(0.0f, p.get(i, k) - threshold + multiplier(i, k));
      return bound;
    };
    float certified_s = s;
    if (use_linear_structure_decoder_) {
      certified_s -= s_x + s_y;
      certified_s += relaxed_structure_bound(
          p_x_support, p_x, x_weight, th_s_[0],
          [&](uint i, uint j) { return gm.q_x(i, j); });
      certified_s += relaxed_structure_bound(
          p_y_support, p_y, y_weight, th_s_[0],
          [&](uint k, uint l) { return gm.q_y(k, l); });
    }
    if (use_linear_alignment_decoder_) {
      certified_s -= s_z;
      certified_s += relaxed_alignment_bound(
          p_z_support, p_z, th_a_,
          [&](uint i, uint k) { return gm.q_z(i, k); });
    }

    if (verbose_ >= 2)
      output_verbose(x, y, z, aln1, aln2);

    // Dynamic CBP generation: add CBPs based on current solution
    if (use_dynamic_cbp_) {
      size_t cbp_before = cbp.size();

      // These are the only new structure coordinates whose multipliers can
      // leave zero in this iteration.  Retaining their history gives an exact
      // support for all potentially positive q_x and q_y entries.
      for (uint i = 0; i < L1; ++i)
        if (x[i] != -1u && x[i] > i)
          add_support(q_x_support, q_x_support_set, i, x[i]);
      for (uint k = 0; k < L2; ++k)
        if (y[k] != -1u && y[k] > k)
          add_support(q_y_support, q_y_support_set, k, y[k]);

      generate_cbp_from_solution(x, y, z, p_x, p_y, p_z, N1, N2, min_th_s, cbp, c_x, c_y, c_z);

      // Exact pricing is required for the restricted Lagrangian to remain an
      // upper bound for the original problem.  Materialize every missing
      // column with positive reduced cost before evaluating L(q).
      const size_t heuristic_cbp_end = cbp.size();
      generate_positive_reduced_cost_cbp(gm, p_x, p_y, p_z,
                                         p_z_forward, p_z_reverse,
                                         q_x_support, q_y_support,
                                         N1, N2, min_th_s,
                                         cbp, c_x, c_y, c_z);
      for (size_t u = cbp_before; u < cbp.size(); ++u) {
        add_support(q_x_support, q_x_support_set,
                    cbp[u].first.first, cbp[u].first.second);
        add_support(q_y_support, q_y_support_set,
                    cbp[u].second.first, cbp[u].second.second);
      }
      cbp_inactive_steps.resize(cbp.size(), 0);
      if (cbp.size() > cbp_before) {
        spdlog::debug("Dynamic CBP: Added {} new CBPs at iteration {} (total: {})",
                      cbp.size() - cbp_before, t, cbp.size());
      }
      if (cbp.size() > heuristic_cbp_end) {
        spdlog::debug("Dynamic CBP pricing: Added {} positive-reduced-cost CBPs at iteration {}",
                      cbp.size() - heuristic_cbp_end, t);
      }
      
      // Sort and remove duplicates in projection arrays after adding new CBPs
      for (uint i = 0; i != c_x.size(); ++i) {
        std::sort(c_x[i].begin(), c_x[i].end());
        c_x[i].erase(std::unique(c_x[i].begin(), c_x[i].end()), c_x[i].end());
      }
      for (uint k = 0; k != c_y.size(); ++k) {
        std::sort(c_y[k].begin(), c_y[k].end());
        c_y[k].erase(std::unique(c_y[k].begin(), c_y[k].end()), c_y[k].end());
      }
      for (uint i = 0; i != c_z.size(); ++i) {
        std::sort(c_z[i].begin(), c_z[i].end());
        c_z[i].erase(std::unique(c_z[i].begin(), c_z[i].end()), c_z[i].end());
      }
    }

    // Calculate Lagrangian value
    w_cbp.clear();
    VF cbp_reduced_cost(cbp.size(), 0.0f);
    for (uint u = 0; u != cbp.size(); ++u)
    {
      const auto &[i, j] = cbp[u].first;
      const auto &[k, l] = cbp[u].second;
      const float s_w = gm.q_x(i, j) + gm.q_y(k, l)
                      - gm.q_z(i, k) - gm.q_z(j, l);
      cbp_reduced_cost[u] = s_w;
      if (s_w > 0.0f)
      {
        s += s_w;
        certified_s += s_w;
        w_cbp.push_back(u);
      }
    }

    best_ub = std::min(best_ub, certified_s);

    VU repaired_x, repaired_y, repaired_z;
    const float feasible_score = repair_feasible_solution(
        repaired_x, repaired_y, repaired_z, x, y, z,
        p_x, p_y, p_z, N1, N2, min_th_s);
    if (feasible_score > best_feasible_score) {
      best_feasible_score = feasible_score;
      best_x.swap(repaired_x);
      best_y.swap(repaired_y);
      best_z.swap(repaired_z);
    }
    if (feasible_score > lb) {
      lb = feasible_score;
      gm.set_lower_bound(lb);
      spdlog::debug("Step: {}, improved feasible LB to {}", t, lb);
    }

    // Update gradients
    violated = gm.update_gradients(cbp, x, y, z, w_cbp, c_x, c_y, c_z,
                                   t, certified_s);

    spdlog::debug("Step: {}, Polyak alpha: {}, BeamL: {}, CertifiedL: {}, BestUB: {}, LB: {}, Violated: {}",
                  t, gm.get_step_size(), s, certified_s, best_ub, lb, violated);

    // Preserve the Lagrangian value from the current iteration, including
    // when the algorithm converges on the first iteration.
    s_prev = s;

    if (use_dynamic_cbp_ && !cbp.empty()) {
      // Keep the dynamic CBP set as a compact working set.  A column is stale
      // only when it has not participated in the auxiliary solution or either
      // folding solution and its reduced cost is non-positive.  Alignment
      // endpoints alone do not activate a consensus base pair.
      // Eight iterations avoids deleting columns during short oscillations.
      constexpr uint CBP_INACTIVE_PATIENCE = 8;
      std::vector<char> selected_w(cbp.size(), false);
      for (const uint u : w_cbp)
        selected_w[u] = true;

      size_t write = 0;
      size_t removed = 0;
      for (size_t u = 0; u < cbp.size(); ++u) {
        const auto& [ij, kl] = cbp[u];
        const auto& [i, j] = ij;
        const auto& [k, l] = kl;
        const bool participates = selected_w[u] || x[i] == j || y[k] == l;
        if (participates || cbp_reduced_cost[u] > 0.0f)
          cbp_inactive_steps[u] = 0;
        else
          ++cbp_inactive_steps[u];

        if (cbp_inactive_steps[u] >= CBP_INACTIVE_PATIENCE &&
            cbp_reduced_cost[u] <= 0.0f) {
          ++removed;
          continue;
        }
        if (write != u) {
          cbp[write] = cbp[u];
          cbp_inactive_steps[write] = cbp_inactive_steps[u];
        }
        ++write;
      }

      if (removed != 0) {
        cbp.resize(write);
        cbp_inactive_steps.resize(write);
        cbp_set_.clear();
        c_x.assign(L1, VU());
        c_y.assign(L2, VU());
        c_z.assign(L1, VU());
        for (const CBP& candidate : cbp) {
          cbp_set_.insert(candidate);
          const auto& [ij, kl] = candidate;
          const auto& [i, j] = ij;
          const auto& [k, l] = kl;
          c_x[i].push_back(j);
          c_y[k].push_back(l);
          c_z[i].push_back(k);
          c_z[j].push_back(l);
        }
        const auto sort_unique = [](VU& projection) {
          std::sort(projection.begin(), projection.end());
          projection.erase(std::unique(projection.begin(), projection.end()),
                           projection.end());
        };
        for (VU& projection : c_x) sort_unique(projection);
        for (VU& projection : c_y) sort_unique(projection);
        for (VU& projection : c_z) sort_unique(projection);
        spdlog::debug("Dynamic CBP: Removed {} stale CBPs at iteration {} (total: {})",
                      removed, t, cbp.size());
      }
    }

    const float certified_gap = best_ub - lb;
    const float gap_tolerance = 1e-4f * std::max(1.0f, std::abs(lb));
    if (violated == 0 ||
        (certified_gap >= 0.0f && certified_gap <= gap_tolerance))
      break; // exact agreement or a sufficiently small certified gap
  }

  if (!best_z.empty()) {
    x = std::move(best_x);
    y = std::move(best_y);
    z = std::move(best_z);
  }
  spdlog::info("Step: {}, BestUB: {}, Violated: {}, LB: {}, Gap: {}",
               t, best_ub, violated, lb, best_ub - lb);

  return std::isfinite(best_ub) ? best_ub : s_prev;
}

float DAFS::
    solve_by_ip(VU &x, VU &y, VU &z,
                const SparseFloatMatrix &p_x,
                const SparseFloatMatrix &p_y,
                const SparseFloatMatrix &p_z,
                const ALN &aln1, const ALN &aln2) const
{
  const uint L1 = p_x.rows();
  const uint L2 = p_y.rows();

  // integer programming
  IP ip(IP::MAX, 1);

  // variables
  VVI v_x(L1, VI(L1, -1));
  VVI v_y(L2, VI(L2, -1));
  VVI v_z(L1, VI(L2, -1));
  VI v_w;

  float min_th_s = th_s_[0];
  for (VF::const_iterator it = th_s_.begin(); it != th_s_.end(); ++it)
    min_th_s = std::min(min_th_s, *it);

  // enumerate the candidates of aligned bases
  for (uint i = 0; i != L1; ++i)
    for (uint k = 0; k != L2; ++k)
      if (p_z.get(i, k) > CUTOFF)
        v_z[i][k] = ip.make_variable(p_z.get(i, k) - th_a_);

  // enumerate the candidates of consensus base-pairs
  std::vector<CBP> cbp;
  for (uint i = 0; i != L1 - 1; ++i)
    for (uint j = i + 1; j != L1; ++j)
      if (p_x.get(i, j) > CUTOFF)
        for (uint k = 0; k != L2 - 1; ++k)
          if (p_z.get(i, k) > CUTOFF)
            for (uint l = k + 1; l != L2; ++l)
              if (p_y.get(k, l) > CUTOFF && p_z.get(j, l) > CUTOFF)
              {
                assert(p_x.get(i, j) <= 1.0);
                assert(p_y.get(k, l) <= 1.0);
                float p = (p_x.get(i, j) + p_y.get(k, l)) / 2;
                float q = (p_z.get(i, k) + p_z.get(j, l)) / 2;
                if (p - min_th_s > 0.0 && w_ * (p - min_th_s) + (q - th_a_) > 0.0)
                {
                  cbp.push_back(std::make_pair(std::make_pair(i, j), std::make_pair(k, l)));
                  v_w.push_back(ip.make_variable(0.0));
                  if (v_x[i][j] < 0)
                    v_x[i][j] = ip.make_variable(w_ * (p_x.get(i, j) - min_th_s));
                  if (v_y[k][l] < 0)
                    v_y[k][l] = ip.make_variable(w_ * (p_y.get(k, l) - min_th_s));
                }
              }
  ip.update();

  // constraints: each base is paired with at most one base (in a)
  for (uint i = 0; i < L1; ++i)
  {
    int row = ip.make_constraint(IP::UP, 0, 1);
    for (uint j = 0; j < i; ++j)
      if (v_x[j][i] >= 0)
        ip.add_constraint(row, v_x[j][i], 1);
    for (uint j = i + 1; j < L1; ++j)
      if (v_x[i][j] >= 0)
        ip.add_constraint(row, v_x[i][j], 1);
  }

  // constraints: no pseudoknots are allowed (in a)
  for (uint i = 0; i < L1 - 1; ++i)
    for (uint j = i + 1; j < L1; ++j)
      if (v_x[i][j] >= 0)
        for (uint k = i + 1; k < j; ++k)
          for (uint l = j + 1; l < L1; ++l)
            if (v_x[k][l] >= 0)
            {
              int row = ip.make_constraint(IP::UP, 0, 1);
              ip.add_constraint(row, v_x[i][j], 1);
              ip.add_constraint(row, v_x[k][l], 1);
            }

  // constraints: each base is paired with at most one base (in b)
  for (uint i = 0; i < L2; ++i)
  {
    int row = ip.make_constraint(IP::UP, 0, 1);
    for (uint j = 0; j < i; ++j)
      if (v_y[j][i] >= 0)
        ip.add_constraint(row, v_y[j][i], 1);
    for (uint j = i + 1; j < L2; ++j)
      if (v_y[i][j] >= 0)
        ip.add_constraint(row, v_y[i][j], 1);
  }

  // constraints: no pseudoknots are allowed (in b)
  for (uint i = 0; i < L2 - 1; ++i)
    for (uint j = i + 1; j < L2; ++j)
      if (v_y[i][j] >= 0)
        for (uint k = i + 1; k < j; ++k)
          for (uint l = j + 1; l < L2; ++l)
            if (v_y[k][l] >= 0)
            {
              int row = ip.make_constraint(IP::UP, 0, 1);
              ip.add_constraint(row, v_y[i][j], 1);
              ip.add_constraint(row, v_y[k][l], 1);
            }

  // constraints: each base is aligned with at most one base
  for (uint i = 0; i < L1; ++i)
  {
    int row = ip.make_constraint(IP::UP, 0, 1);
    for (uint k = 0; k < L2; ++k)
      if (v_z[i][k] >= 0)
        ip.add_constraint(row, v_z[i][k], 1);
  }
  for (uint k = 0; k < L2; ++k)
  {
    int row = ip.make_constraint(IP::UP, 0, 1);
    for (uint i = 0; i < L1; ++i)
      if (v_z[i][k] >= 0)
        ip.add_constraint(row, v_z[i][k], 1);
  }

  // constraints: no crossing matches are allowed
  for (uint i = 0; i < L1; ++i)
    for (uint k = 0; k < L2; ++k)
      if (v_z[i][k] >= 0)
        for (uint j = i + 1; j < L1; ++j)
          for (uint l = 0; l < k; ++l)
            if (v_z[j][l] >= 0)
            {
              int row = ip.make_constraint(IP::UP, 0, 1);
              ip.add_constraint(row, v_z[i][k], 1);
              ip.add_constraint(row, v_z[j][l], 1);
            }

  // constraints for consensus base pairs
  VVI r_x(L1, VI(L1, -1));
  for (uint i = 0; i < L1 - 1; ++i)
    for (uint j = i + 1; j < L1; ++j)
      if (v_x[i][j] >= 0)
      {
        r_x[i][j] = ip.make_constraint(IP::FX, 0, 0);
        ip.add_constraint(r_x[i][j], v_x[i][j], 1);
      }

  VVI r_y(L2, VI(L2, -1));
  for (uint i = 0; i < L2 - 1; ++i)
    for (uint j = i + 1; j < L2; ++j)
      if (v_y[i][j] >= 0)
      {
        r_y[i][j] = ip.make_constraint(IP::FX, 0, 0);
        ip.add_constraint(r_y[i][j], v_y[i][j], 1);
      }

  VVI r_z(L1, VI(L2, -1));
  for (uint i = 0; i < L1; ++i)
    for (uint k = 0; k < L2; ++k)
      if (v_z[i][k] >= 0)
      {
        r_z[i][k] = ip.make_constraint(IP::LO, 0, 0);
        ip.add_constraint(r_z[i][k], v_z[i][k], 1);
      }

  for (uint u = 0; u != cbp.size(); ++u)
  {
    const uint i = cbp[u].first.first, j = cbp[u].first.second;
    const uint k = cbp[u].second.first, l = cbp[u].second.second;
    assert(r_x[i][j] >= 0 && v_x[i][j] >= 0);
    ip.add_constraint(r_x[i][j], v_w[u], -1);
    assert(r_y[k][l] >= 0 && v_y[k][l] >= 0);
    ip.add_constraint(r_y[k][l], v_w[u], -1);
    assert(r_z[i][k] >= 0 && v_z[i][k] >= 0);
    ip.add_constraint(r_z[i][k], v_w[u], -1);
    assert(r_z[j][l] >= 0 && v_z[j][l] >= 0);
    ip.add_constraint(r_z[j][l], v_w[u], -1);
  }

  // execute optimization
  float s = ip.solve();

  // build the result
  x.resize(L1);
  std::fill(x.begin(), x.end(), -1u);
  for (uint i = 0; i < L1 - 1; ++i)
    for (uint j = i + 1; j < L1; ++j)
      if (v_x[i][j] >= 0 && ip.get_value(v_x[i][j]) > 0.5)
        x[i] = j;

  y.resize(L2);
  std::fill(y.begin(), y.end(), -1u);
  for (uint i = 0; i < L2 - 1; ++i)
    for (uint j = i + 1; j < L2; ++j)
      if (v_y[i][j] >= 0 && ip.get_value(v_y[i][j]) > 0.5)
        y[i] = j;

  z.resize(L1);
  std::fill(z.begin(), z.end(), -1u);
  for (uint i = 0; i < L1; ++i)
    for (uint k = 0; k < L2; ++k)
      if (v_z[i][k] >= 0 && ip.get_value(v_z[i][k]) > 0.5)
        z[i] = k;

  return s;
}

void DAFS::
    align(ALN &aln, int ch)
{
  if (tree_[ch].second.first == -1u)
  {
    assert(tree_[ch].second.second == -1u);
    aln.resize(1);
    aln[0].first = ch;
    aln[0].second = std::vector<bool>(fa_[ch].size(), true);
  }
  else
  {
    ALN aln1, aln2;
    align(aln1, tree_[ch].second.first);
    align(aln2, tree_[ch].second.second);
    align_alignments(aln, aln1, aln2);
  }
}

float DAFS::
    align(VU &ss, ALN &aln, int ch)
{
  float s = 0.0;
  if (tree_[ch].second.first == -1u)
  {
    assert(tree_[ch].second.second == -1u);
    aln.resize(1);
    aln[0].first = ch;
    aln[0].second = std::vector<bool>(fa_[ch].size(), true);
    spdlog::info("Aligning #{} ({})", ch, fa_[ch].name());
  }
  else
  {
    ALN aln1, aln2;
    VU ss1, ss2;
    align(ss1, aln1, tree_[ch].second.first);
    align(ss2, aln2, tree_[ch].second.second);
    spdlog::info("Aligning #{} and #{} into #{}", tree_[ch].second.first, tree_[ch].second.second, ch);
    s = align_alignments(ss, aln, aln1, aln2);

    // iterative refinement
    if (aln.size() >= 3 && n_refinement_ > 0) 
    {
      spdlog::info("Refining alignment for {} sequences", aln.size());
      auto n_itr = std::min(n_refinement_, (uint)aln.size());
      for (uint i = 0; i != n_itr; ++i)
      {
        VU ss_temp = ss;
        ALN aln_temp = aln;
        float s_temp;

        s_temp = refine(ss_temp, aln_temp);
        //std::cout << s << " " << s_temp << std::endl;
        if (s_temp > s)
        {
          s = s_temp;
          std::swap(ss, ss_temp);
          std::swap(aln, aln_temp);
        }
      }
    }
  }
  return s;
}

float DAFS::
    refine(VU &ss, ALN &aln)
{
  VU group[2];
  VU idx(aln.size(), -1u);
  for (auto i=0u; i!=idx.size(); ++i)
    idx[i] = i;
  std::shuffle(idx.begin(), idx.end(), g_);
  auto mid = idx.size() / 2;
  for (auto i=0; i!=mid; ++i)
    group[0].push_back(idx[i]);
  for (auto i=mid; i!=idx.size(); ++i)
    group[1].push_back(idx[i]);

  ALN a[2];
  for (uint i = 0; i != 2; ++i)
  {
    uint N = group[i].size();
    a[i].resize(N);
    uint L = aln[group[i][0]].second.size();
    for (uint j = 0; j != N; ++j)
      a[i][j].first = aln[group[i][j]].first;
    for (uint k = 0; k != L; ++k)
    {
      bool gap = true;
      for (uint j = 0; j != N; ++j)
        gap &= !aln[group[i][j]].second[k];
      if (!gap)
      {
        for (uint j = 0; j != N; ++j)
          a[i][j].second.push_back(aln[group[i][j]].second[k]);
      }
    }
  }

  ALN r;
  float s = align_alignments(ss, r, a[0], a[1]);
  std::swap(r, aln);
  return s;
}

void DAFS::
    output(std::ostream &os, const ALN &aln) const
{
  output(os, aln.begin(), aln.end());
}

void DAFS::
    output(std::ostream &os, ALN::const_iterator b, ALN::const_iterator e) const
{
  for (ALN::const_iterator a = b; a != e; ++a)
  {
    uint s = a->first;
    os << ">"
       << " " << fa_[s].name() << std::endl;
    for (uint j = 0, k = 0; j != a->second.size(); ++j)
    {
      if (a->second[j])
        os << fa_[s].seq()[k++];
      else
        os << '-';
    }
    os << std::endl;
  }
}

DAFS&
DAFS::
parse_options(int& argc, char**& argv)
{
  cxxopts::Options options(argv[0], "DAFS: dual decomposition for simultaneous aligning and folding RNA sequences.");
  options.add_options()
    ("h,help", "Print usage")
    ("version", "Print version")
    ("input", "Input file", cxxopts::value<std::string>(), "FILE")
    ("r,refinement", "The number of iteration of the iterative refinment", cxxopts::value<int>()->default_value("0"), "N")
    ("seed", "Random seed for shuffling", cxxopts::value<int>()->default_value("42"), "N")
    ("w,weight", "Weight of the expected accuracy score for secondary structures", cxxopts::value<float>()->default_value("4.0"))
    ("eta", "Initial step width for the subgradient optimization", cxxopts::value<float>()->default_value("0.5"))
    ("m,max-iter", "The maximum number of iteration of the subgradient optimization", cxxopts::value<int>()->default_value("600"), "T")
    ("f,fourway-pct", "Weight of four-way PCT", cxxopts::value<float>()->default_value("0.0"))
    ("v,verbose", "The level of verbose outputs", cxxopts::value<int>()->default_value("0"))
    ("dynamic-cbp", "Use dynamic CBP generation instead of pre-enumeration")
    ("dense-lagrangian", "Force dense Lagrange multiplier storage")
    ;

  options.add_options("Aligning")
    ("a,align-model", "Alignment model for calcualating matching probablities (value=CONTRAlign, ProbCons, LinearAlign)", 
      cxxopts::value<std::string>()->default_value("ProbCons"))
    ("align-beam", "Beam size for LinearAlign model", cxxopts::value<int>()->default_value("100"))
    ("p,align-pct", "Weight of PCT for matching probabilities", cxxopts::value<float>()->default_value("0.25"))
    ("u,align-th", "Threshold for matching probabilities", cxxopts::value<float>()->default_value("0.01"))
    ("align-aux", "Load matching probability matrices from FILENAME", cxxopts::value<std::string>(), "FILENAME");

  options.add_options("Folding")
    ("s,fold-model", "Folding model for calculating base-pairing probablities (value=Boltzmann, Vienna, CONTRAfold, lpv, lpc)",
      cxxopts::value<std::string>()->default_value("Boltzmann"))
    ("fold-decoder", "Decoder for common secondary structure prediction (value=Nussinov, IPknot)",
      cxxopts::value<std::string>()->default_value("Nussinov"))
    ("q,fold-pct", "Weight of PCT for base-pairing probabilities", cxxopts::value<float>()->default_value("0.25"))
    ("t,fold-th", "Threshold for base-pairing probabilities", cxxopts::value<std::vector<float>>()->default_value("0.2"))
    ("g,gamma", "Specify the threshold for base-pairing probabilities by 1/(gamma+1))", cxxopts::value<std::vector<float>>())
    ("no-alifold", "No use of RNAalifold for calculating base-pairing probabilities")
    ("T,fold-th1", "Threshold for base-pairing probabilities of the conclusive common secondary structures", cxxopts::value<std::vector<float>>())
    ("G,gamma1", "Specify the threshold for base-pairing probabilities of the conclusive common secondary structuresby 1/(gamma+1))", cxxopts::value<std::vector<float>>())
    ("ipknot", "Set optimized parameters for IPknot decoding (--fold-decoder=IPknot -g4,8 -G2,4 --bp-update1)")
    ("bp-update", "Use the iterative update of BPs")
    ("bp-update1", "Use the iterative update of BPs for the final prediction")
    ("fold-aux", "Load base-pairing probability matrices from FILENAME", cxxopts::value<std::string>(), "FILENAME")
    ("linfold-beam", "Beam size for LinFold algorithm", cxxopts::value<int>()->default_value("100"), "SIZE");

  options.parse_positional({"input"});
  options.positional_help("FILE").show_positional_help();

  try
  {
    auto res = options.parse(argc, argv);
    if (res.count("help"))
    {
      std::cout << options.help({"", "Aligning", "Folding"}) << std::endl;
      exit(0);
    }
    if (res.count("version"))
    {
      std::cout << "DAFS version " << VERSION << std::endl;
      exit(0);
    }
    // general options
    n_refinement_ = res["refinement"].as<int>();
    g_ = std::mt19937(res["seed"].as<int>());
    w_ = res["weight"].as<float>();
    eta0_ = res["eta"].as<float>();
    t_max_ = res["max-iter"].as<int>();
    w_pct_f_ = res["fourway-pct"].as<float>();
    verbose_ = res["verbose"].as<int>();
    use_dynamic_cbp_ = res.count("dynamic-cbp") > 0;
    switch (verbose_)
    {
    default:
    case 0:
      spdlog::set_level(spdlog::level::warn);
      break;
    case 1:
      spdlog::set_level(spdlog::level::info);
      break;
    case 2:
      spdlog::set_level(spdlog::level::debug);
      break;
    }

    // options for alignments
    w_pct_a_ = res["align-pct"].as<float>();
    th_a_ = res["align-th"].as<float>();
    const std::string align_model = res["align-model"].as<std::string>();
    const int requested_align_beam = res["align-beam"].as<int>();
    const uint align_beam = std::max(1, requested_align_beam);
    if (requested_align_beam <= 0 && align_model == "LinearAlign")
      spdlog::warn("--align-beam must be positive for linear complexity; using 1");

    if (res["align-aux"].count())
      a_model_ = std::make_unique<AUXAlign>(res["align-aux"].as<std::string>(), CUTOFF);
    else if (align_model == "CONTRAlign")
      a_model_ = std::make_unique<CONTRAlign>(th_a_);
    else if (align_model == "ProbCons")
      a_model_ = std::make_unique<ProbCons>(th_a_);
    else if (align_model == "LinearAlign")
      a_model_ = std::make_unique<LinearAlign>(th_a_, align_beam);
    else
      throw "Unknown alignment model: " + align_model;
    assert(a_model_);
    use_linear_alignment_decoder_ = !res["align-aux"].count() &&
                                    align_model == "LinearAlign";
    if (use_linear_alignment_decoder_)
      a_decoder_ = std::make_unique<LinearNeedlemanWunsch>(
          th_a_, align_beam);
    else
      a_decoder_ = std::make_unique<SparseNeedlemanWunsch>(th_a_);

    // options for folding
    w_pct_s_ = res["fold-pct"].as<float>();
    use_alifold_ = res["no-alifold"].count() == 0;
    const std::string fold_model = res["fold-model"].as<std::string>();
    const int requested_fold_beam = res["linfold-beam"].as<int>();
    const uint fold_beam = std::max(1, requested_fold_beam);
    if (requested_fold_beam <= 0 &&
        (fold_model == "lpv" || fold_model == "lpc" || fold_model == "LinFold"))
      spdlog::warn("--linfold-beam must be positive for linear complexity; using 1");
    if (res["fold-aux"].count())
      s_model_ = std::make_unique<AUXFold>(res["fold-aux"].as<std::string>(), CUTOFF);
    else if (fold_model == "Boltzmann")
      s_model_ = std::make_unique<RNAfold>(true, nullptr, CUTOFF);
    else if (fold_model == "Vienna")
      s_model_ = std::make_unique<RNAfold>(false, nullptr, CUTOFF);
    else if (fold_model == "CONTRAfold")
      s_model_ = std::make_unique<CONTRAfold>(CUTOFF);
    else if (fold_model == "lpv" || fold_model == "LinFold")
    {
      s_model_ = std::make_unique<LinFoldWrapper>(CUTOFF, LinFoldWrapper::ModelType::LPV, fold_beam);
    }
    else if (fold_model == "lpc")
    {
      s_model_ = std::make_unique<LinFoldWrapper>(CUTOFF, LinFoldWrapper::ModelType::LPC, fold_beam);
    }
    else
      throw "Unknown folding model: " + fold_model;
    assert(s_model_);
    const bool linear_probability_folding = !res["fold-aux"].count() &&
        (fold_model == "lpv" || fold_model == "lpc" ||
         fold_model == "LinFold");
    use_linear_profile_folding_ = linear_probability_folding;
    // --no-alifold applies both during progressive alignment and during the
    // final consensus fold.  In LinearPartition mode "alifold" denotes the
    // linear consensus-profile surrogate, not Vienna RNAalifold.
    use_alifold1_ = use_alifold_;

    if (res["fold-th"].count())
    {
      th_s_ = res["fold-th"].as<std::vector<float>>();
    }
    else if (res["gamma"].count())
    {
      th_s_ = res["gamma"].as<std::vector<float>>();
      for (uint i = 0; i != th_s_.size(); ++i)
        th_s_[i] = 1.0 / (1.0 + th_s_[i]);
    }
    else if (res["ipknot"].count())
    {
      th_s_.resize(2);
      th_s_[0] = 1.0 / (1.0 + 4.0);
      th_s_[1] = 1.0 / (1.0 + 8.0);
    }
    else
    {
      th_s_ = res["fold-th"].as<std::vector<float>>();
    }

    VF th_s1;
    if (res["fold-th1"].count())
    {
      th_s1 = res["fold-th1"].as<std::vector<float>>();
    }
    else if (res["gamma1"].count())
    {
      th_s1 = res["gamma1"].as<std::vector<float>>();
      for (uint i = 0; i != th_s1.size(); ++i)
        th_s1[i] = 1.0 / (1.0 + th_s1[i]);
    }
    else if (res["ipknot"].count())
    {
      th_s1.resize(2);
      th_s1[0] = 1.0 / (1.0 + 2.0);
      th_s1[1] = 1.0 / (1.0 + 4.0);
    }
    else
    {
      th_s1 = th_s_;
    }

    const std::string fold_decoder = res["fold-decoder"].as<std::string>();
    if (fold_decoder == "IPknot" || res["ipknot"].count())
    {
      s_decoder_ = std::make_unique<IPknot>(th_s_);
      s_decoder1_ = std::make_unique<IPknot>(th_s1);
    }
    else if (fold_decoder == "Nussinov")
    {
      const bool linear_folding = fold_model == "lpv" ||
                                  fold_model == "lpc" ||
                                  fold_model == "LinFold";
      use_linear_structure_decoder_ = !res["fold-aux"].count() && linear_folding;
      if (use_linear_structure_decoder_) {
        s_decoder_ = std::make_unique<LinearNussinov>(th_s_[0], fold_beam);
        s_decoder1_ = std::make_unique<LinearNussinov>(th_s1[0], fold_beam);
      } else {
        s_decoder_ = std::make_unique<SparseNussinov>(th_s_[0]);
        s_decoder1_ = std::make_unique<SparseNussinov>(th_s1[0]);
      }
    }
    else
      throw "Unknown folding decoder: " + res["fold-decoder"].as<std::string>();
    assert(s_decoder_);

    // Sparse multiplier lookup has expected O(1) cost but a larger constant
    // than direct dense indexing.  Select it independently for probability
    // components produced by Linear models; other components retain the
    // original dense decoder path unchanged.
    const bool linear_folding = fold_model == "lpv" ||
                                fold_model == "lpc" ||
                                fold_model == "LinFold";
    const bool force_dense = res["dense-lagrangian"].count();
    use_sparse_structure_lagrangian_ = !force_dense &&
                                        !res["fold-aux"].count() &&
                                        use_linear_structure_decoder_ &&
                                        fold_decoder == "Nussinov" &&
                                        !res["ipknot"].count();
    use_sparse_alignment_lagrangian_ = !force_dense &&
                                        !res["align-aux"].count() &&
                                        use_linear_alignment_decoder_;
    spdlog::info("DAFS decoder: structure={}, alignment={}",
                 use_linear_structure_decoder_ ? "beam-max" : "exact-DP",
                 use_linear_alignment_decoder_ ? "beam-max" : "exact-DP");
    spdlog::info("Lagrangian storage: structure={}, alignment={}",
                 use_sparse_structure_lagrangian_ ? "sparse" : "dense",
                 use_sparse_alignment_lagrangian_ ? "sparse" : "dense");
    spdlog::info("Profile folding: {}",
                 !use_alifold_ ? "disabled" :
                 (use_linear_profile_folding_ ? "linear-consensus" :
                                                "Vienna-RNAalifold"));

    use_bp_update_ = res["bp-update"].count() > 0;
    use_bp_update1_ = res["bp-update1"].count() > 0 ^ res["ipknot"].count() > 0;

    // read sequences
    Fasta::load(fa_, res["input"].as<std::string>().c_str());
  }
  catch (const cxxopts::exceptions::exception& e)
  {
    std::cout << options.help() << std::endl;
    exit(0);
  }

  return *this;
}

int DAFS::
    run()
{
  const uint N = fa_.size();

  // calculate base-pairing probabilities
  s_model_->calculate(fa_, bp_);
#if 0
  {
    std::ofstream os("bp");
    save_bp(os, bp_);
  }
#endif

  // calculate matching probabilities
  a_model_->calculate(fa_, mp_);
  for (uint i = 0; i != N; ++i)
    for (uint j = i + 1; j != N; ++j)
      transpose_mp(mp_[i][j], mp_[j][i], fa_[i].size(), fa_[j].size());
#if 0
  {
    std::ofstream os("mp");
    save_mp(os, mp_);
  }
#endif

  // four-way probabilistic consistency tranformation
  if (w_pct_f_ != 0.0)
    relax_fourway_consistency();

  // calculate probabilistic similarity scores
  // which are used for building guide trees and PCTs
  sim_.resize(N, VF(N));
  const auto* linear_similarity_decoder = use_linear_alignment_decoder_
      ? dynamic_cast<const LinearNeedlemanWunsch*>(a_decoder_.get())
      : nullptr;
  assert(!use_linear_alignment_decoder_ || linear_similarity_decoder != nullptr);
  for (uint i = 0; i != N; ++i)
  {
    sim_[i][i] = 1.0;
    for (uint j = i + 1; j != N; ++j)
      sim_[i][j] = sim_[j][i] = linear_similarity_decoder
          ? linear_similarity_decoder->similarity_score(
                mp_[i][j], fa_[i].size(), fa_[j].size())
          : calculate_similarity_score(
                mp_[i][j], fa_[i].size(), fa_[j].size());
  }

  // probabilistic consistency tranformation for base-pairing probabilitiy matrix
  if (w_pct_s_ != 0.0)
    relax_basepairing_probability();

  // probabilistic consistency tranformation for matching probability matrix
  if (w_pct_a_ != 0.0)
    relax_matching_probability();

  // compute the guide tree
  build_tree();
  print_tree(std::cout, tree_.size() - 1);
  std::cout << std::endl;

  // compute progressive alignments along with the guide tree
  VU ss;
  ALN aln;
  float s;
  s = align(ss, aln, tree_.size() - 1);

#if 0
  // iterative refinement
  for (uint i = 0; i != n_refinement_; ++i)
  {
    VU ss_temp = ss;
    ALN aln_temp = aln;
    float s_temp;

    s_temp = refine(ss_temp, aln_temp);
    //std::cout << s << " " << s_temp << std::endl;
    if (s_temp > s)
    {
      s = s_temp;
      std::swap(ss, ss_temp);
      std::swap(aln, aln_temp);
    }
  }
#endif

  std::string str;
  if (s_decoder1_)
  {
    // compute the common secondary structures from the averaged base-pairing matrix
    SparseFloatMatrix p;
    average_basepairing_probability(p, aln, use_alifold1_);
    if (use_bp_update1_)
    {
      std::string str;
      VU ss;
      s_decoder1_->decode(p, ss, str);
      update_basepairing_probability(p, ss, str, aln, use_alifold1_);
    }
    s_decoder1_->decode(p, ss, str);
  }
  else
    s_decoder_->make_brackets(ss, str);

  // output the alignment
  std::sort(aln.begin(), aln.end());
  std::cout << ">SS_cons" << std::endl
            << str << std::endl;
  output(std::cout, aln);

#if 0
  Alifold ali(0.0);
  float cv;
  std::cout << ali.energy_of_struct(aln, fa_, str, cv) << std::endl
            << cv << std::endl;
#endif

  return 0;
}

int main(int argc, char *argv[])
{
  try
  {
    DAFS dafs;
    return dafs.parse_options(argc, argv).run();
  }
  catch (const char *str)
  {
    std::cerr << str << std::endl;
  }
  catch (std::string str)
  {
    std::cerr << str << std::endl;
  }
  catch (std::system_error e)
  {
    std::cerr << e.what() << std::endl;
  }
  return EXIT_FAILURE;
}

// Dynamic CBP generation methods implementation

bool DAFS::is_valid_cbp(uint i, uint j, uint k, uint l, 
                        const SparseFloatMatrix& p_x,
                        const SparseFloatMatrix& p_y,
                        const SparseFloatMatrix& p_z,
                        uint N1, uint N2, float min_th_s) const {
    // Use the same logic as the original dafs.cpp:995-1001
    if (p_x.get(i, j) > CUTOFF && p_z.get(i, k) > CUTOFF &&
        p_y.get(k, l) > CUTOFF && p_z.get(j, l) > CUTOFF) {
        
        assert(p_x.get(i, j) <= 1.0);
        assert(p_y.get(k, l) <= 1.0);
        float p = (N1 * p_x.get(i, j) + N2 * p_y.get(k, l)) / (N1 + N2);
        float q = (p_z.get(i, k) + p_z.get(j, l)) / 2;
        return (p - min_th_s > 0.0 && w_ * (p - min_th_s) + (q - th_a_) > 0.0);
    }
    return false;
}

void DAFS::add_cbp_if_new(const CBP& candidate, std::vector<CBP>& cbp, 
                          VVU& c_x, VVU& c_y, VVU& c_z) {
    if (cbp_set_.insert(candidate).second) {  // Only add if it's new
        cbp.push_back(candidate);
        auto [i, j] = candidate.first;
        auto [k, l] = candidate.second;
        
        c_x[i].push_back(j);
        c_y[k].push_back(l);
        c_z[i].push_back(k);
        c_z[j].push_back(l);
    }
}

void DAFS::generate_cbp_from_solution(const VU& x, const VU& y, const VU& z,
                                      const SparseFloatMatrix& p_x,
                                      const SparseFloatMatrix& p_y,
                                      const SparseFloatMatrix& p_z,
                                      uint N1, uint N2, float min_th_s,
                                      std::vector<CBP>& cbp, VVU& c_x, VVU& c_y, VVU& c_z) {
    const uint L1 = p_x.rows();
    const uint L2 = p_y.rows();

    // Forward separation: project each base pair selected in x through the
    // current alignment.  Do not require y to have selected the projected base
    // pair; that disagreement is what the multipliers need to resolve.
    for (uint i = 0; i < L1-1; ++i) {
        const uint j = x[i];
        if (j == -1u || j <= i)
            continue;
        const uint k = z[i];
        const uint l = z[j];
        if (k != -1u && l != -1u && k < l &&
            is_valid_cbp(i, j, k, l, p_x, p_y, p_z,
                         N1, N2, min_th_s)) {
            add_cbp_if_new({{i, j}, {k, l}}, cbp, c_x, c_y, c_z);
        }
    }

    // Reverse separation catches the symmetric case: a base pair selected in
    // y whose aligned endpoints are not selected as a base pair in x.
    VU inverse_z(L2, -1u);
    for (uint i = 0; i < L1; ++i)
        if (z[i] != -1u && z[i] < L2)
            inverse_z[z[i]] = i;

    for (uint k = 0; k < L2-1; ++k) {
        const uint l = y[k];
        if (l == -1u || l <= k)
            continue;
        const uint i = inverse_z[k];
        const uint j = inverse_z[l];
        if (i != -1u && j != -1u && i < j &&
            is_valid_cbp(i, j, k, l, p_x, p_y, p_z,
                         N1, N2, min_th_s)) {
            add_cbp_if_new({{i, j}, {k, l}}, cbp, c_x, c_y, c_z);
        }
    }
}

void DAFS::generate_positive_reduced_cost_cbp(
    const GradientManager& gm,
    const SparseFloatMatrix& p_x,
    const SparseFloatMatrix& p_y,
    const SparseFloatMatrix& p_z,
    const VVU& p_z_forward, const VVU& p_z_reverse,
    const std::vector<std::pair<uint, uint>>& q_x_support,
    const std::vector<std::pair<uint, uint>>& q_y_support,
    uint N1, uint N2, float min_th_s,
    std::vector<CBP>& cbp, VVU& c_x, VVU& c_y, VVU& c_z) {
    const uint L1 = p_x.rows();
    const uint L2 = p_y.rows();
    if (L1 < 2 || L2 < 2)
        return;

    const auto price_candidate = [&](uint i, uint j, uint k, uint l) {
        if (i >= L1 || j >= L1 || k >= L2 || l >= L2 || i >= j || k >= l)
            return;
        const CBP candidate = {{i, j}, {k, l}};
        if (cbp_set_.find(candidate) != cbp_set_.end())
            return;
        if (!is_valid_cbp(i, j, k, l, p_x, p_y, p_z,
                          N1, N2, min_th_s))
            return;
        const float reduced_cost = gm.q_x(i, j) + gm.q_y(k, l)
                                 - gm.q_z(i, k) - gm.q_z(j, l);
        if (reduced_cost > 0.0f)
            add_cbp_if_new(candidate, cbp, c_x, c_y, c_z);
    };

    // Since q_z is projected onto the non-negative orthant, a positive
    // reduced cost implies q_x(i,j)>0 or q_y(k,l)>0.  Searching from both
    // positive supports is therefore exhaustive, while the two alignment
    // adjacency lists avoid the four-dimensional dense scan.
    for (const auto& [i, j] : q_x_support) {
        if (gm.q_x(i, j) <= 0.0f || p_x.get(i, j) <= CUTOFF)
            continue;
        for (const uint k : p_z_forward[i])
            for (const uint l : p_z_forward[j])
                if (k < l && p_y.get(k, l) > CUTOFF)
                    price_candidate(i, j, k, l);
    }

    for (const auto& [k, l] : q_y_support) {
        if (gm.q_y(k, l) <= 0.0f || p_y.get(k, l) <= CUTOFF)
            continue;
        for (const uint i : p_z_reverse[k])
            for (const uint j : p_z_reverse[l])
                if (i < j && p_x.get(i, j) > CUTOFF)
                    price_candidate(i, j, k, l);
    }
}
