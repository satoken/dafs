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
#include <unordered_map>
#include <cstdint>
#include <limits>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <sstream>
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
#include "relaxed_bounds.h"
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

namespace
{
using Clock = std::chrono::steady_clock;

double elapsed_seconds(const Clock::time_point& start)
{
  return std::chrono::duration<double>(Clock::now() - start).count();
}

std::string json_escape(const std::string& value)
{
  std::ostringstream out;
  for (const unsigned char ch : value) {
    switch (ch) {
    case '"': out << "\\\""; break;
    case '\\': out << "\\\\"; break;
    case '\b': out << "\\b"; break;
    case '\f': out << "\\f"; break;
    case '\n': out << "\\n"; break;
    case '\r': out << "\\r"; break;
    case '\t': out << "\\t"; break;
    default:
      if (ch < 0x20)
        out << "\\u" << std::hex << std::setw(4) << std::setfill('0')
            << static_cast<unsigned>(ch) << std::dec;
      else
        out << ch;
    }
  }
  return out.str();
}

std::string json_number(double value)
{
  if (!std::isfinite(value))
    return "null";
  std::ostringstream out;
  out << std::setprecision(17) << value;
  return out.str();
}

size_t sparse_entries(const std::vector<SV>& matrix)
{
  size_t count = 0;
  for (const SV& row : matrix)
    count += row.size();
  return count;
}
} // namespace

DAFS::DAFS()
    : w_ribosum_(0.075f),
      align_probability_beam_(100),
      align_dd_beam_(100),
      fold_probability_beam_(100),
      fold_dd_beam_(100),
      fold_final_beam_(100),
      metrics_merge_id_(0),
      metrics_dd_iterations_(0),
      metrics_cbp_peak_(0),
      metrics_cbp_added_(0),
      metrics_cbp_priced_(0),
      metrics_cbp_removed_(0),
      metrics_dd_seconds_(0.0),
      use_dynamic_cbp_(false),
      use_sparse_structure_lagrangian_(false),
      use_sparse_alignment_lagrangian_(false),
      use_linear_structure_decoder_(false),
      use_linear_alignment_decoder_(false),
      use_alifold_(false),
      use_alifold1_(false),
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
#if defined(WITH_GLPK) || defined(WITH_CPLEX) || defined(WITH_GUROBI) || \
    defined(WITH_SCIP) || defined(WITH_HIGHS)
  return t_max_ != 0 ? solve_by_dd(x, y, z, p_x, p_y, p_z, aln1, aln2) : solve_by_ip(x, y, z, p_x, p_y, p_z, aln1, aln2);
#else
  if (t_max_ == 0)
    throw std::runtime_error(
        "--max-iter=0 requires an integer-programming backend "
        "(GLPK, CPLEX, Gurobi, SCIP, or HiGHS)");
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
  Alifold(0.0 /*CUTOFF*/).fold(aln, fa_, bp);
}

void DAFS::
    calculate_profile_basepairing_probability(
        const ALN &aln, const std::string &constraint, BP &bp) const
{
  Alifold(0.0 /*CUTOFF*/).fold(aln, fa_, constraint, bp);
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
                             uint N1, uint N2, float min_th_s,
                             const std::function<float(
                                 uint, uint, uint, uint)>& pair_match_score,
                             float& intersection_score_out,
                             float& consensus_score_out) const
{
  const uint L1 = p_x.rows();
  const uint L2 = p_y.rows();
  repaired_x.assign(L1, -1u);
  repaired_y.assign(L2, -1u);
  repaired_z = z;

  // Every decoded z is a feasible monotone alignment.  Its original (not
  // Lagrangian-shifted) contribution is therefore a valid primal score.
  float alignment_score = 0.0f;
  for (uint i = 0; i < L1; ++i) {
    const uint k = z[i];
    if (k != -1u)
      alignment_score += p_z.get(i, k) - th_a_;
  }
  float intersection_score = alignment_score;

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
    const float pair_match = pair_match_score(i, j, k, l);
    if (!is_valid_cbp(i, j, k, l, p_x, p_y, p_z,
                      N1, N2, min_th_s, pair_match))
      continue;

    repaired_x[i] = j;
    repaired_y[k] = l;
    intersection_score += x_weight * (p_x.get(i, j) - repair_th);
    intersection_score += y_weight * (p_y.get(k, l) - repair_th);
    intersection_score += pair_match;
  }

  intersection_score_out = intersection_score;
  consensus_score_out = intersection_score;
  if (!use_linear_structure_decoder_)
    return intersection_score;

  // The decoded x/y intersection above is safe but unnecessarily
  // conservative.  For a fixed monotone z, any non-crossing matching on the
  // first profile maps to a non-crossing matching on the second profile.
  // Re-optimize those mapped consensus-pair contributions with the linear
  // Nussinov decoder.  Threshold sparsity bounds the candidate scan, and a
  // fixed folding beam keeps the repair linear in profile length.
  SparseFloatMatrix consensus_pair_score;
  consensus_pair_score.assign(L1, L1);
  for (uint i = 0; i < L1; ++i) {
    const uint k = z[i];
    if (k == -1u || k >= L2)
      continue;
    for (const auto& [j, probability] : p_x.ordered_row(i)) {
      if (j <= i + 2 || j >= L1)
        continue;
      const uint l = z[j];
      if (l == -1u || l >= L2 || l <= k + 2)
        continue;
      const float pair_match = pair_match_score(i, j, k, l);
      if (!is_valid_cbp(i, j, k, l, p_x, p_y, p_z,
                        N1, N2, min_th_s, pair_match))
        continue;
      const float contribution =
          x_weight * (probability - repair_th) +
          y_weight * (p_y.get(k, l) - repair_th) + pair_match;
      if (contribution > 0.0f)
        consensus_pair_score.set(i, j, contribution);
    }
  }

  LinearNussinov repair_decoder(0.0f, fold_dd_beam_);
  VU consensus_x;
  std::string brackets;
  const float consensus_structure_score =
      repair_decoder.decode(consensus_pair_score, consensus_x, brackets);
  const float consensus_score = alignment_score + consensus_structure_score;
  consensus_score_out = consensus_score;
  if (consensus_score > intersection_score) {
    VU consensus_y(L2, -1u);
    for (uint i = 0; i < L1; ++i) {
      const uint j = consensus_x[i];
      if (j == -1u || j <= i)
        continue;
      const uint k = z[i];
      const uint l = z[j];
      assert(k != -1u && l != -1u && k < l);
      consensus_y[k] = l;
    }
    repaired_x.swap(consensus_x);
    repaired_y.swap(consensus_y);
    return consensus_score;
  }

  return intersection_score;
}

float DAFS::
    solve_by_dd(VU &x, VU &y, VU &z,
                const SparseFloatMatrix &p_x,
                const SparseFloatMatrix &p_y,
                const SparseFloatMatrix &p_z,
                const ALN &aln1, const ALN &aln2)
{
  const auto dd_start = Clock::now();
  const uint merge_id = ++metrics_merge_id_;
  const uint L1 = p_x.rows();
  const uint L2 = p_y.rows();
  const uint N1 = aln1.size();
  const uint N2 = aln2.size();
  const RibosumProfile ribosum_x(aln1, fa_);
  const RibosumProfile ribosum_y(aln2, fa_);
  std::unordered_map<CBP, float, CBPHash> pair_match_cache;
  const auto pair_match_score = [&](uint i, uint j, uint k, uint l) {
    if (w_ribosum_ == 0.0f)
      return 0.0f;
    const CBP candidate = {{i, j}, {k, l}};
    const auto found = pair_match_cache.find(candidate);
    if (found != pair_match_cache.end())
      return found->second;
    const float score = w_ribosum_ *
        ribosum_x.pair_score(i, j, ribosum_y, k, l);
    pair_match_cache.emplace(candidate, score);
    return score;
  };

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
    if (support_set.insert(support_key(i, j)).second) {
      support.emplace_back(i, j);
      return true;
    }
    return false;
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

  // Linear Nussinov needs the union of probability and multiplier supports.
  // Probability support is fixed and multiplier support only grows.  Keep the
  // right-indexed union incrementally instead of scanning two sparse matrices
  // and sorting/deduplicating their entries in every DD iteration.
  VVU x_pairs_by_right(L1), y_pairs_by_right(L2);
  for (const auto& [i, j] : p_x_support)
    x_pairs_by_right[j].push_back(i);
  for (const auto& [k, l] : p_y_support)
    y_pairs_by_right[l].push_back(k);
  const auto add_decoder_support = [](VVU& pairs_by_right,
                                      uint left, uint right) {
    assert(right < pairs_by_right.size());
    VU& lefts = pairs_by_right[right];
    const auto position = std::lower_bound(lefts.begin(), lefts.end(), left);
    if (position == lefts.end() || *position != left)
      lefts.insert(position, left);
  };
  
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
                  const float ribosum = pair_match_score(i, j, k, l);
                  if (is_valid_cbp(i, j, k, l, p_x, p_y, p_z,
                                   N1, N2, min_th_s, ribosum))
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
  } else if (w_ribosum_ > 0.0f) {
    // A positive intrinsic w coefficient can have positive reduced cost even
    // while every multiplier is zero.  Materialize all such columns once and
    // retain them, preserving exact pricing without a dense four-dimensional
    // scan in every iteration.  Sparse probability degrees bound this scan by
    // O(|p_x| / th_a^2), hence it is linear in sequence length for fixed
    // thresholds.
    for (const auto& [i, j] : p_x_support) {
      for (const uint k : p_z_forward[i]) {
        for (const uint l : p_z_forward[j]) {
          if (k >= l || p_y.get(k, l) <= CUTOFF)
            continue;
          const float ribosum = pair_match_score(i, j, k, l);
          if (ribosum > 0.0f &&
              is_valid_cbp(i, j, k, l, p_x, p_y, p_z,
                           N1, N2, min_th_s, ribosum))
            add_cbp_if_new({{i, j}, {k, l}}, cbp, c_x, c_y, c_z);
        }
      }
    }
    for (VU& projection : c_x) {
      std::sort(projection.begin(), projection.end());
      projection.erase(std::unique(projection.begin(), projection.end()),
                       projection.end());
    }
    for (VU& projection : c_y) {
      std::sort(projection.begin(), projection.end());
      projection.erase(std::unique(projection.begin(), projection.end()),
                       projection.end());
    }
    for (VU& projection : c_z) {
      std::sort(projection.begin(), projection.end());
      projection.erase(std::unique(projection.begin(), projection.end()),
                       projection.end());
    }
  }

  const size_t initial_cbp = cbp.size();
  metrics_cbp_peak_ = std::max(metrics_cbp_peak_, initial_cbp);
  if (use_dynamic_cbp_)
    metrics_cbp_added_ += initial_cbp;

  // Initialize gradient manager
  GradientManager gm(eta0_, 0.0f, 1.0f,
                     use_sparse_structure_lagrangian_,
                     use_sparse_alignment_lagrangian_);
  gm.initialize(L1, L2);
  // Beam subgradients are useful for primal recovery but do not minimize the
  // certified matching relaxation.  Maintain an independent, always-sparse
  // dual track whose subproblems are the convex left/row relaxations.  A fixed
  // number of tracks preserves linear complexity.
  const bool use_certified_dual_track =
      use_linear_structure_decoder_ || use_linear_alignment_decoder_;
  std::unique_ptr<GradientManager> certified_gm;
  if (use_certified_dual_track) {
    certified_gm = std::make_unique<GradientManager>(
        eta0_, 0.0f, 1.0f, true, true);
    certified_gm->initialize(L1, L2);
  }
  
  // Certified bounds and the best feasible primal solution recovered so far.
  // The former alignment-only bound could contain a consensus pair across a
  // gap in one profile; use only explicitly verified coupled solutions here.
  float lb = std::numeric_limits<float>::lowest();
  float best_ub = std::numeric_limits<float>::infinity();
  float best_feasible_score = std::numeric_limits<float>::lowest();
  VU best_x, best_y, best_z;
  
  float s_prev = 0.0;
  uint violated = 0;
  uint certified_violated = 0;
  uint t;
  uint iterations_completed = 0;
  std::string stop_reason = "max_iterations";
  for (t = 0; t != t_max_; ++t)
  {
    size_t iteration_added = 0;
    size_t iteration_priced = 0;
    size_t iteration_removed = 0;
    // solve the subproblems
    float s_x = 0.0f, s_y = 0.0f, s_z = 0.0f;
    const float x_weight = w_ * 2 * N1 / (N1 + N2);
    const float y_weight = w_ * 2 * N2 / (N1 + N2);
    LinearNussinovResult x_linear_result, y_linear_result;
    if (use_linear_structure_decoder_) {
      auto* linear_decoder = dynamic_cast<LinearNussinov*>(s_decoder_.get());
      assert(linear_decoder);
      if (gm.uses_sparse_structure_storage()) {
        x_linear_result = linear_decoder->decode_certified(
            x_weight, p_x, gm.sparse_q_x(), x_pairs_by_right, x);
        y_linear_result = linear_decoder->decode_certified(
            y_weight, p_y, gm.sparse_q_y(), y_pairs_by_right, y);
      } else {
        x_linear_result = linear_decoder->decode_certified(
            x_weight, p_x, gm.dense_q_x(), x_pairs_by_right, x);
        y_linear_result = linear_decoder->decode_certified(
            y_weight, p_y, gm.dense_q_y(), y_pairs_by_right, y);
      }
      s_x = x_linear_result.score;
      s_y = y_linear_result.score;
    } else if (gm.uses_sparse_structure_storage()) {
      s_x = s_decoder_->decode(x_weight, dense_p_x,
                               gm.sparse_q_x(), x);
      s_y = s_decoder_->decode(y_weight, dense_p_y,
                               gm.sparse_q_y(), y);
    } else {
      s_x = s_decoder_->decode(x_weight, dense_p_x, gm.dense_q_x(), x);
      s_y = s_decoder_->decode(y_weight, dense_p_y, gm.dense_q_y(), y);
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

    // Beam maxima are feasible subproblem values, not certified maxima.  The
    // linear Nussinov decoder records every genuinely pruned state and bounds
    // its best possible completion with feasible vertex-cover potentials.
    // Its upper_bound is therefore already a rigorous bound for the complete
    // sparse subproblem.  Do not rebuild the same edge lists and vertex cover
    // here: that duplicated O(support) work without tightening any long-chain
    // case in the factorial benchmark.
    float certified_s = s;
    RelaxedBounds::AlignmentBound z_bound_details;
    float x_certified_bound = 0.0f;
    float y_certified_bound = 0.0f;
    if (use_linear_structure_decoder_) {
      certified_s -= s_x + s_y;
      x_certified_bound = x_linear_result.upper_bound;
      y_certified_bound = y_linear_result.upper_bound;
      // The beam result is feasible.  This maximum is normally redundant,
      // and protects certification against implementation or rounding drift.
      certified_s += std::max(x_certified_bound, s_x) +
                     std::max(y_certified_bound, s_y);
    }
    if (use_linear_alignment_decoder_) {
      certified_s -= s_z;
      z_bound_details = RelaxedBounds::alignment_bound(
          p_z_support, L1, L2, [&](uint i, uint k) {
            return p_z.get(i, k) - th_a_ + gm.q_z(i, k);
          });
      certified_s += std::max(z_bound_details.best, s_z);
    }

    VU certified_x, certified_y, certified_z, certified_w_cbp;
    float certified_track_s = std::numeric_limits<float>::infinity();
    double certified_track_accumulator =
        std::numeric_limits<double>::infinity();
    if (use_certified_dual_track) {
      certified_track_accumulator = RelaxedBounds::structure_left_solution(
          p_x_support, L1, [&](uint i, uint j) {
            return x_weight * (p_x.get(i, j) - th_s_[0]) -
                   certified_gm->q_x(i, j);
          }, certified_x);
      certified_track_accumulator += RelaxedBounds::structure_left_solution(
          p_y_support, L2, [&](uint k, uint l) {
            return y_weight * (p_y.get(k, l) - th_s_[0]) -
                   certified_gm->q_y(k, l);
          }, certified_y);
      certified_track_accumulator += RelaxedBounds::alignment_row_solution(
          p_z_support, L1, L2, [&](uint i, uint k) {
            return p_z.get(i, k) - th_a_ + certified_gm->q_z(i, k);
          }, certified_z);
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
        if (x[i] != -1u && x[i] > i &&
            add_support(q_x_support, q_x_support_set, i, x[i]))
          add_decoder_support(x_pairs_by_right, i, x[i]);
      for (uint k = 0; k < L2; ++k)
        if (y[k] != -1u && y[k] > k &&
            add_support(q_y_support, q_y_support_set, k, y[k]))
          add_decoder_support(y_pairs_by_right, k, y[k]);
      if (use_certified_dual_track) {
        for (uint i = 0; i < L1; ++i)
          if (certified_x[i] != -1u && certified_x[i] > i &&
              add_support(q_x_support, q_x_support_set,
                          i, certified_x[i]))
            add_decoder_support(x_pairs_by_right, i, certified_x[i]);
        for (uint k = 0; k < L2; ++k)
          if (certified_y[k] != -1u && certified_y[k] > k &&
              add_support(q_y_support, q_y_support_set,
                          k, certified_y[k]))
            add_decoder_support(y_pairs_by_right, k, certified_y[k]);
      }

      generate_cbp_from_solution(x, y, z, p_x, p_y, p_z,
                                 N1, N2, min_th_s,
                                 ribosum_x, ribosum_y,
                                 cbp, c_x, c_y, c_z);

      // Exact pricing is required for the restricted Lagrangian to remain an
      // upper bound for the original problem.  Materialize every missing
      // column with positive reduced cost under either dual track before
      // evaluating L(q).  A single union scan is exhaustive for both tracks.
      const size_t heuristic_cbp_end = cbp.size();
      generate_positive_reduced_cost_cbp(
          gm, certified_gm.get(), p_x, p_y, p_z,
          p_z_forward, p_z_reverse,
          q_x_support, q_y_support,
          N1, N2, min_th_s,
          ribosum_x, ribosum_y,
          cbp, c_x, c_y, c_z);
      for (size_t u = cbp_before; u < cbp.size(); ++u) {
        const auto [i, j] = cbp[u].first;
        const auto [k, l] = cbp[u].second;
        if (add_support(q_x_support, q_x_support_set, i, j))
          add_decoder_support(x_pairs_by_right, i, j);
        if (add_support(q_y_support, q_y_support_set, k, l))
          add_decoder_support(y_pairs_by_right, k, l);
      }
      cbp_inactive_steps.resize(cbp.size(), 0);
      iteration_added = heuristic_cbp_end - cbp_before;
      iteration_priced = cbp.size() - heuristic_cbp_end;
      metrics_cbp_added_ += iteration_added + iteration_priced;
      metrics_cbp_priced_ += iteration_priced;
      metrics_cbp_peak_ = std::max(metrics_cbp_peak_, cbp.size());
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
    VF certified_cbp_reduced_cost(cbp.size(), 0.0f);
    for (uint u = 0; u != cbp.size(); ++u)
    {
      const auto &[i, j] = cbp[u].first;
      const auto &[k, l] = cbp[u].second;
      const float s_w = pair_match_score(i, j, k, l)
                      + gm.q_x(i, j) + gm.q_y(k, l)
                      - gm.q_z(i, k) - gm.q_z(j, l);
      cbp_reduced_cost[u] = s_w;
      if (s_w > 0.0f)
      {
        s += s_w;
        certified_s += s_w;
        w_cbp.push_back(u);
      }
      if (use_certified_dual_track) {
        const float certified_s_w = pair_match_score(i, j, k, l)
            + certified_gm->q_x(i, j) + certified_gm->q_y(k, l)
            - certified_gm->q_z(i, k) - certified_gm->q_z(j, l);
        certified_cbp_reduced_cost[u] = certified_s_w;
        if (certified_s_w > 0.0f) {
          certified_track_accumulator += static_cast<double>(certified_s_w);
          certified_w_cbp.push_back(u);
        }
      }
    }

    const float update_s = s;
    if (use_certified_dual_track)
      certified_track_s =
          RelaxedBounds::round_up_to_float(certified_track_accumulator);
    best_ub = std::min({best_ub, certified_s, certified_track_s});

    VU repaired_x, repaired_y, repaired_z;
    float intersection_feasible_score = 0.0f;
    float consensus_feasible_score = 0.0f;
    const float feasible_score = repair_feasible_solution(
        repaired_x, repaired_y, repaired_z, x, y, z,
        p_x, p_y, p_z, N1, N2, min_th_s, pair_match_score,
        intersection_feasible_score, consensus_feasible_score);
    if (feasible_score > best_feasible_score) {
      best_feasible_score = feasible_score;
      best_x.swap(repaired_x);
      best_y.swap(repaired_y);
      best_z.swap(repaired_z);
    }
    if (feasible_score > lb) {
      lb = feasible_score;
      gm.set_lower_bound(lb);
      if (use_certified_dual_track)
        certified_gm->set_lower_bound(lb);
      spdlog::debug("Step: {}, improved feasible LB to {}", t, lb);
    }

    // Update gradients
    // The subgradient comes from the beam solutions x/y/z/w, so its Polyak
    // numerator must use their Lagrangian value.  The relaxed certified value
    // remains reserved for UB tracking and stopping guarantees.
    violated = gm.update_gradients(cbp, x, y, z, w_cbp, c_x, c_y, c_z,
                                   t, update_s);
    certified_violated = use_certified_dual_track
        ? certified_gm->update_gradients(
              cbp, certified_x, certified_y, certified_z, certified_w_cbp,
              c_x, c_y, c_z, t, certified_track_s)
        : 0;

    spdlog::debug("Step: {}, Polyak alpha: {}, BeamL: {}, CertifiedL: {}, BestUB: {}, LB: {}, Violated: {}",
                  t, gm.get_step_size(), update_s, certified_s, best_ub, lb, violated);

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
      for (const uint u : certified_w_cbp)
        selected_w[u] = true;

      size_t write = 0;
      size_t removed = 0;
      for (size_t u = 0; u < cbp.size(); ++u) {
        const auto& [ij, kl] = cbp[u];
        const auto& [i, j] = ij;
        const auto& [k, l] = kl;
        const float ribosum = pair_match_score(i, j, k, l);
        const bool participates = selected_w[u] || x[i] == j || y[k] == l ||
            (use_certified_dual_track &&
             (certified_x[i] == j || certified_y[k] == l));
        if (participates || cbp_reduced_cost[u] > 0.0f ||
            (use_certified_dual_track &&
             certified_cbp_reduced_cost[u] > 0.0f))
          cbp_inactive_steps[u] = 0;
        else
          ++cbp_inactive_steps[u];

        if (ribosum <= 0.0f &&
            cbp_inactive_steps[u] >= CBP_INACTIVE_PATIENCE &&
            cbp_reduced_cost[u] <= 0.0f &&
            (!use_certified_dual_track ||
             certified_cbp_reduced_cost[u] <= 0.0f)) {
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
      iteration_removed = removed;
      metrics_cbp_removed_ += removed;
    }

    const float certified_gap = best_ub - lb;
    const float gap_tolerance = 1e-4f * std::max(1.0f, std::abs(lb));
    ++iterations_completed;
    ++metrics_dd_iterations_;
    write_metric(
        "dd_iteration",
        {{"merge_id", std::to_string(merge_id)},
         {"iteration", std::to_string(t)},
         {"beam_lagrangian", json_number(s)},
         {"polyak_lagrangian", json_number(update_s)},
         {"certified_lagrangian", json_number(certified_s)},
         {"certified_track_lagrangian", json_number(certified_track_s)},
         {"certified_track_violated", std::to_string(certified_violated)},
         {"certified_track_polyak_update", use_certified_dual_track
              ? json_number(certified_gm->get_last_update_size()) : "null"},
         {"x_beam", json_number(s_x)},
         {"x_bound", use_linear_structure_decoder_
                         ? json_number(x_certified_bound) : "null"},
         {"x_bound_beam_certificate", use_linear_structure_decoder_
              ? json_number(x_linear_result.upper_bound) : "null"},
         {"x_bound_pruned_certificate", use_linear_structure_decoder_
              ? json_number(x_linear_result.pruned_upper_bound) : "null"},
         {"x_bound_additive_certificate", use_linear_structure_decoder_
              ? json_number(x_linear_result.additive_upper_bound) : "null"},
         {"x_pruned_states", use_linear_structure_decoder_
              ? std::to_string(x_linear_result.pruned_states) : "null"},
         {"x_bound_all", "null"},
         {"x_bound_left", "null"},
         {"x_bound_right", "null"},
         {"x_bound_incident", "null"},
         {"x_bound_cardinality", "null"},
         {"x_bound_top_k", "null"},
         {"x_bound_vertex_cover", "null"},
         {"y_beam", json_number(s_y)},
         {"y_bound", use_linear_structure_decoder_
                         ? json_number(y_certified_bound) : "null"},
         {"y_bound_beam_certificate", use_linear_structure_decoder_
              ? json_number(y_linear_result.upper_bound) : "null"},
         {"y_bound_pruned_certificate", use_linear_structure_decoder_
              ? json_number(y_linear_result.pruned_upper_bound) : "null"},
         {"y_bound_additive_certificate", use_linear_structure_decoder_
              ? json_number(y_linear_result.additive_upper_bound) : "null"},
         {"y_pruned_states", use_linear_structure_decoder_
              ? std::to_string(y_linear_result.pruned_states) : "null"},
         {"y_bound_all", "null"},
         {"y_bound_left", "null"},
         {"y_bound_right", "null"},
         {"y_bound_incident", "null"},
         {"y_bound_cardinality", "null"},
         {"y_bound_top_k", "null"},
         {"y_bound_vertex_cover", "null"},
         {"z_beam", json_number(s_z)},
         {"z_bound", use_linear_alignment_decoder_
                         ? json_number(z_bound_details.best) : "null"},
         {"z_bound_row", use_linear_alignment_decoder_
                             ? json_number(z_bound_details.row) : "null"},
         {"z_bound_column", use_linear_alignment_decoder_
                                ? json_number(z_bound_details.column) : "null"},
         {"feasible_repair", json_number(feasible_score)},
         {"intersection_repair", json_number(intersection_feasible_score)},
         {"consensus_repair", json_number(consensus_feasible_score)},
         {"best_ub", json_number(best_ub)},
         {"lb", json_number(lb)},
         {"gap", json_number(certified_gap)},
         {"violated", std::to_string(violated)},
         {"polyak_scale", json_number(gm.get_step_size())},
         {"polyak_update", json_number(gm.get_last_update_size())},
         {"q_x_nnz", gm.uses_sparse_structure_storage()
                         ? std::to_string(gm.sparse_q_x().nonzeros()) : "null"},
         {"q_y_nnz", gm.uses_sparse_structure_storage()
                         ? std::to_string(gm.sparse_q_y().nonzeros()) : "null"},
         {"q_z_nnz", gm.uses_sparse_alignment_storage()
                         ? std::to_string(gm.sparse_q_z().nonzeros()) : "null"},
         {"certified_q_x_nnz", use_certified_dual_track
              ? std::to_string(certified_gm->sparse_q_x().nonzeros()) : "null"},
         {"certified_q_y_nnz", use_certified_dual_track
              ? std::to_string(certified_gm->sparse_q_y().nonzeros()) : "null"},
         {"certified_q_z_nnz", use_certified_dual_track
              ? std::to_string(certified_gm->sparse_q_z().nonzeros()) : "null"},
         {"cbp_total", std::to_string(cbp.size())},
         {"cbp_heuristic_added", std::to_string(iteration_added)},
         {"cbp_priced", std::to_string(iteration_priced)},
         {"cbp_removed", std::to_string(iteration_removed)},
         {"selected_pair_matches", std::to_string(w_cbp.size())}});
    if (certified_gap >= 0.0f && certified_gap <= gap_tolerance) {
      stop_reason = "certified_gap";
      break;
    }
    if (violated == 0 &&
        (!use_certified_dual_track || certified_violated == 0)) {
      // With exact subproblem oracles, agreement proves optimality.  A beam
      // oracle can only establish that its own feasible solutions agree; the
      // distinct reason prevents that stationary condition being reported as
      // a certified result when the relaxed UB still has a gap.
      stop_reason = use_linear_structure_decoder_ ||
                    use_linear_alignment_decoder_
                  ? "beam_stationary" : "agreement";
      break;
    }
  }

  if (!best_z.empty()) {
    x = std::move(best_x);
    y = std::move(best_y);
    z = std::move(best_z);
  }
  spdlog::info("Step: {}, BestUB: {}, Violated: {}, LB: {}, Gap: {}",
               t, best_ub, violated, lb, best_ub - lb);

  const double dd_seconds = elapsed_seconds(dd_start);
  metrics_dd_seconds_ += dd_seconds;
  write_metric(
      "dd_summary",
      {{"merge_id", std::to_string(merge_id)},
       {"iterations", std::to_string(iterations_completed)},
       {"seconds", json_number(dd_seconds)},
       {"initial_cbp", std::to_string(initial_cbp)},
       {"final_cbp", std::to_string(cbp.size())},
       {"best_ub", json_number(best_ub)},
       {"lb", json_number(lb)},
       {"gap", json_number(best_ub - lb)},
       {"violated", std::to_string(violated)}},
      {{"stop_reason", stop_reason}});

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
  const uint N1 = aln1.size();
  const uint N2 = aln2.size();
  const RibosumProfile ribosum_x(aln1, fa_);
  const RibosumProfile ribosum_y(aln2, fa_);
  const float x_weight = w_ * 2.0f * N1 / (N1 + N2);
  const float y_weight = w_ * 2.0f * N2 / (N1 + N2);

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
                const float ribosum = w_ribosum_ *
                    ribosum_x.pair_score(i, j, ribosum_y, k, l);
                if (is_valid_cbp(i, j, k, l, p_x, p_y, p_z,
                                 N1, N2, min_th_s, ribosum))
                {
                  cbp.push_back(std::make_pair(std::make_pair(i, j), std::make_pair(k, l)));
                  v_w.push_back(ip.make_variable(ribosum));
                  if (v_x[i][j] < 0)
                    v_x[i][j] = ip.make_variable(
                        x_weight * (p_x.get(i, j) - min_th_s));
                  if (v_y[k][l] < 0)
                    v_y[k][l] = ip.make_variable(
                        y_weight * (p_y.get(k, l) - min_th_s));
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

void DAFS::write_metric(
    const std::string& event,
    const std::vector<std::pair<std::string, std::string>>& raw_fields,
    const std::vector<std::pair<std::string, std::string>>& text_fields) const
{
  if (!metrics_stream_.is_open())
    return;

  metrics_stream_ << "{\"schema_version\":1,\"event\":\""
                  << json_escape(event) << '"';
  for (const auto& [name, value] : raw_fields)
    metrics_stream_ << ",\"" << json_escape(name) << "\":" << value;
  for (const auto& [name, value] : text_fields)
    metrics_stream_ << ",\"" << json_escape(name) << "\":\""
                    << json_escape(value) << '"';
  metrics_stream_ << "}\n";
  // Per-iteration flushing would materially perturb the runtime being
  // measured.  Keep iteration events buffered, but make stage/summary events
  // available promptly so interrupted runs still retain useful diagnostics.
  if (event != "dd_iteration")
    metrics_stream_.flush();
  if (!metrics_stream_)
    throw std::runtime_error("failed to write benchmark metrics to " +
                             metrics_jsonl_path_);
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
    ("ribosum-weight", "Weight of RIBOSUM85-60 pair-pair match scores", cxxopts::value<float>()->default_value("0.075"))
    ("eta", "Initial step width for the subgradient optimization", cxxopts::value<float>()->default_value("0.5"))
    ("m,max-iter", "The maximum number of iteration of the subgradient optimization", cxxopts::value<int>()->default_value("600"), "T")
    ("f,fourway-pct", "Weight of four-way PCT", cxxopts::value<float>()->default_value("0.0"))
    ("v,verbose", "The level of verbose outputs", cxxopts::value<int>()->default_value("0"))
    ("dynamic-cbp", "Use dynamic CBP generation instead of pre-enumeration")
    ("dense-lagrangian", "Force dense Lagrange multiplier storage")
    ("metrics-jsonl", "Write machine-readable benchmark metrics to FILE",
      cxxopts::value<std::string>(), "FILE")
    ;

  options.add_options("Aligning")
    ("a,align-model", "Alignment model for calcualating matching probablities (value=CONTRAlign, ProbCons, LinearAlign)", 
      cxxopts::value<std::string>()->default_value("ProbCons"))
    ("align-beam", "Beam size for LinearAlign probability calculation", cxxopts::value<int>()->default_value("100"))
    ("align-dd-beam", "Beam size for the LinearAlign DD decoder (default: --align-beam)", cxxopts::value<int>())
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
    ("alifold", "Use RNAalifold for profile base-pairing probabilities (disabled by default)")
    ("no-alifold", "Disable RNAalifold for profile base-pairing probabilities (default)")
    ("T,fold-th1", "Threshold for base-pairing probabilities of the conclusive common secondary structures", cxxopts::value<std::vector<float>>())
    ("G,gamma1", "Specify the threshold for base-pairing probabilities of the conclusive common secondary structuresby 1/(gamma+1))", cxxopts::value<std::vector<float>>())
    ("ipknot", "Set optimized parameters for IPknot decoding (--fold-decoder=IPknot -g4,8 -G2,4 --bp-update1)")
    ("bp-update", "Use the iterative update of BPs")
    ("bp-update1", "Use the iterative update of BPs for the final prediction")
    ("fold-aux", "Load base-pairing probability matrices from FILENAME", cxxopts::value<std::string>(), "FILENAME")
    ("linfold-beam", "Beam size for LinearPartition probability calculation", cxxopts::value<int>()->default_value("100"), "SIZE")
    ("fold-dd-beam", "Beam size for the folding DD decoder (default: --linfold-beam)", cxxopts::value<int>(), "SIZE")
    ("fold-final-beam", "Beam size for final consensus folding (default: --linfold-beam)", cxxopts::value<int>(), "SIZE");

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
    w_ribosum_ = res["ribosum-weight"].as<float>();
    if (!std::isfinite(w_ribosum_) || w_ribosum_ < 0.0f)
      throw std::invalid_argument("--ribosum-weight must be finite and non-negative");
    eta0_ = res["eta"].as<float>();
    const int requested_max_iter = res["max-iter"].as<int>();
    if (requested_max_iter < 0)
      throw std::invalid_argument("--max-iter must be non-negative");
    t_max_ = static_cast<uint>(requested_max_iter);
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
    align_model_name_ = align_model;
    const int requested_align_beam = res["align-beam"].as<int>();
    const int requested_align_dd_beam = res.count("align-dd-beam")
        ? res["align-dd-beam"].as<int>() : requested_align_beam;
    const uint align_probability_beam = std::max(1, requested_align_beam);
    const uint align_dd_beam = std::max(1, requested_align_dd_beam);
    align_probability_beam_ = align_probability_beam;
    align_dd_beam_ = align_dd_beam;
    if (requested_align_beam <= 0 && align_model == "LinearAlign")
      spdlog::warn("--align-beam must be positive for linear complexity; using 1");
    if (requested_align_dd_beam <= 0 && align_model == "LinearAlign")
      spdlog::warn("--align-dd-beam must be positive for linear complexity; using 1");

    if (res["align-aux"].count())
      a_model_ = std::make_unique<AUXAlign>(res["align-aux"].as<std::string>(), CUTOFF);
    else if (align_model == "CONTRAlign")
      a_model_ = std::make_unique<CONTRAlign>(th_a_);
    else if (align_model == "ProbCons")
      a_model_ = std::make_unique<ProbCons>(th_a_);
    else if (align_model == "LinearAlign")
      a_model_ = std::make_unique<LinearAlign>(th_a_, align_probability_beam);
    else
      throw "Unknown alignment model: " + align_model;
    assert(a_model_);
    use_linear_alignment_decoder_ = !res["align-aux"].count() &&
                                    align_model == "LinearAlign";
    if (use_linear_alignment_decoder_)
      a_decoder_ = std::make_unique<LinearNeedlemanWunsch>(
          th_a_, align_dd_beam);
    else
      a_decoder_ = std::make_unique<SparseNeedlemanWunsch>(th_a_);

    // options for folding
    w_pct_s_ = res["fold-pct"].as<float>();
    if (res.count("alifold") && res.count("no-alifold"))
      throw std::invalid_argument("--alifold and --no-alifold are mutually exclusive");
    use_alifold_ = res.count("alifold") > 0;
    const std::string fold_model = res["fold-model"].as<std::string>();
    fold_model_name_ = fold_model;
    const int requested_fold_beam = res["linfold-beam"].as<int>();
    const int requested_fold_dd_beam = res.count("fold-dd-beam")
        ? res["fold-dd-beam"].as<int>() : requested_fold_beam;
    const int requested_fold_final_beam = res.count("fold-final-beam")
        ? res["fold-final-beam"].as<int>() : requested_fold_beam;
    const uint fold_probability_beam = std::max(1, requested_fold_beam);
    const uint fold_dd_beam = std::max(1, requested_fold_dd_beam);
    const uint fold_final_beam = std::max(1, requested_fold_final_beam);
    fold_probability_beam_ = fold_probability_beam;
    fold_dd_beam_ = fold_dd_beam;
    fold_final_beam_ = fold_final_beam;
    if (requested_fold_beam <= 0 &&
        (fold_model == "lpv" || fold_model == "lpc" || fold_model == "LinFold"))
      spdlog::warn("--linfold-beam must be positive for linear complexity; using 1");
    if (requested_fold_dd_beam <= 0 &&
        (fold_model == "lpv" || fold_model == "lpc" || fold_model == "LinFold"))
      spdlog::warn("--fold-dd-beam must be positive for linear complexity; using 1");
    if (requested_fold_final_beam <= 0 &&
        (fold_model == "lpv" || fold_model == "lpc" || fold_model == "LinFold"))
      spdlog::warn("--fold-final-beam must be positive for linear complexity; using 1");
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
      s_model_ = std::make_unique<LinFoldWrapper>(CUTOFF, LinFoldWrapper::ModelType::LPV, fold_probability_beam);
    }
    else if (fold_model == "lpc")
    {
      s_model_ = std::make_unique<LinFoldWrapper>(CUTOFF, LinFoldWrapper::ModelType::LPC, fold_probability_beam);
    }
    else
      throw "Unknown folding model: " + fold_model;
    assert(s_model_);
    const bool linear_probability_folding = !res["fold-aux"].count() &&
        (fold_model == "lpv" || fold_model == "lpc" ||
         fold_model == "LinFold");
    // A majority-consensus sequence discards compensatory substitutions and
    // proved both slower and less accurate than the already available mean of
    // the per-sequence BPPs.  Vienna RNAalifold would break linearity, so a
    // linear probability engine deliberately disables both profile surrogates.
    if (linear_probability_folding && use_alifold_) {
      spdlog::info("LinearPartition profile folding disabled; using averaged per-sequence BPPs");
      use_alifold_ = false;
    }
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
        s_decoder_ = std::make_unique<LinearNussinov>(th_s_[0], fold_dd_beam);
        s_decoder1_ = std::make_unique<LinearNussinov>(th_s1[0], fold_final_beam);
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
                 !use_alifold_ ? "disabled" : "Vienna-RNAalifold");
    spdlog::info("RIBOSUM85-60 pair-pair weight: {}", w_ribosum_);

    use_bp_update_ = res["bp-update"].count() > 0;
    use_bp_update1_ = res["bp-update1"].count() > 0 ^ res["ipknot"].count() > 0;

    // Read sequences only after all options have been validated.  Metrics are
    // opened last so a successful parse always produces a self-contained
    // trace, without mixing it into the prediction written to stdout.
    input_path_ = res["input"].as<std::string>();
    Fasta::load(fa_, input_path_.c_str());
    if (res["metrics-jsonl"].count()) {
      metrics_jsonl_path_ = res["metrics-jsonl"].as<std::string>();
      metrics_stream_.open(metrics_jsonl_path_, std::ios::out | std::ios::trunc);
      if (!metrics_stream_)
        throw std::runtime_error("cannot open --metrics-jsonl file: " +
                                 metrics_jsonl_path_);
    }
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
  const auto run_start = Clock::now();
  const uint N = fa_.size();
  size_t total_residues = 0;
  size_t minimum_length = N == 0 ? 0 : std::numeric_limits<size_t>::max();
  size_t maximum_length = 0;
  for (const Fasta& sequence : fa_) {
    total_residues += sequence.size();
    minimum_length = std::min(minimum_length,
                              static_cast<size_t>(sequence.size()));
    maximum_length = std::max(maximum_length,
                              static_cast<size_t>(sequence.size()));
  }
  metrics_merge_id_ = 0;
  metrics_dd_iterations_ = 0;
  metrics_cbp_peak_ = 0;
  metrics_cbp_added_ = 0;
  metrics_cbp_priced_ = 0;
  metrics_cbp_removed_ = 0;
  metrics_dd_seconds_ = 0.0;
  write_metric(
      "run_start",
      {{"sequence_count", std::to_string(N)},
       {"total_residues", std::to_string(total_residues)},
       {"minimum_length", std::to_string(minimum_length)},
       {"maximum_length", std::to_string(maximum_length)},
       {"structure_weight", json_number(w_)},
       {"ribosum_weight", json_number(w_ribosum_)},
       {"alignment_threshold", json_number(th_a_)},
       {"max_iterations", std::to_string(t_max_)},
       {"alignment_beam", std::to_string(align_probability_beam_)},
       {"alignment_dd_beam", std::to_string(align_dd_beam_)},
       {"folding_beam", std::to_string(fold_probability_beam_)},
       {"folding_dd_beam", std::to_string(fold_dd_beam_)},
       {"folding_final_beam", std::to_string(fold_final_beam_)},
       {"alifold", use_alifold_ ? "true" : "false"},
       {"dynamic_cbp", use_dynamic_cbp_ ? "true" : "false"},
       {"sparse_structure_lagrangian",
        use_sparse_structure_lagrangian_ ? "true" : "false"},
       {"sparse_alignment_lagrangian",
        use_sparse_alignment_lagrangian_ ? "true" : "false"}},
      {{"input", input_path_},
       {"alignment_model", align_model_name_},
       {"folding_model", fold_model_name_}});

  // calculate base-pairing probabilities
  auto stage_start = Clock::now();
  s_model_->calculate(fa_, bp_);
  size_t base_pair_probability_nnz = 0;
  for (const BP& probability : bp_)
    base_pair_probability_nnz += sparse_entries(probability);
  write_metric("stage",
               {{"seconds", json_number(elapsed_seconds(stage_start))},
                {"nnz", std::to_string(base_pair_probability_nnz)}},
               {{"name", "base_pair_probabilities"}});
#if 0
  {
    std::ofstream os("bp");
    save_bp(os, bp_);
  }
#endif

  // calculate matching probabilities
  stage_start = Clock::now();
  a_model_->calculate(fa_, mp_);
  size_t matching_probability_nnz = 0;
  for (uint i = 0; i != N; ++i)
    for (uint j = i + 1; j != N; ++j)
      matching_probability_nnz += sparse_entries(mp_[i][j]);
  for (uint i = 0; i != N; ++i)
    for (uint j = i + 1; j != N; ++j)
      transpose_mp(mp_[i][j], mp_[j][i], fa_[i].size(), fa_[j].size());
  write_metric("stage",
               {{"seconds", json_number(elapsed_seconds(stage_start))},
                {"nnz", std::to_string(matching_probability_nnz)}},
               {{"name", "alignment_probabilities"}});
#if 0
  {
    std::ofstream os("mp");
    save_mp(os, mp_);
  }
#endif

  // four-way probabilistic consistency tranformation
  stage_start = Clock::now();
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
  write_metric("stage",
               {{"seconds", json_number(elapsed_seconds(stage_start))}},
               {{"name", "similarity_and_fourway_pct"}});

  // probabilistic consistency tranformation for base-pairing probabilitiy matrix
  stage_start = Clock::now();
  if (w_pct_s_ != 0.0)
    relax_basepairing_probability();

  // probabilistic consistency tranformation for matching probability matrix
  if (w_pct_a_ != 0.0)
    relax_matching_probability();
  write_metric("stage",
               {{"seconds", json_number(elapsed_seconds(stage_start))}},
               {{"name", "probabilistic_consistency"}});
  size_t final_base_pair_probability_nnz = 0;
  for (const BP& probability : bp_)
    final_base_pair_probability_nnz += sparse_entries(probability);
  size_t final_matching_probability_nnz = 0;
  for (uint i = 0; i != N; ++i)
    for (uint j = i + 1; j != N; ++j)
      final_matching_probability_nnz += sparse_entries(mp_[i][j]);
  write_metric(
      "probability_support",
      {{"base_pair_nnz", std::to_string(final_base_pair_probability_nnz)},
       {"alignment_nnz", std::to_string(final_matching_probability_nnz)}});

  // compute the guide tree
  stage_start = Clock::now();
  build_tree();
  write_metric("stage",
               {{"seconds", json_number(elapsed_seconds(stage_start))}},
               {{"name", "guide_tree"}});
  print_tree(std::cout, tree_.size() - 1);
  std::cout << std::endl;

  // compute progressive alignments along with the guide tree
  VU ss;
  ALN aln;
  float s;
  stage_start = Clock::now();
  s = align(ss, aln, tree_.size() - 1);
  write_metric("stage",
               {{"seconds", json_number(elapsed_seconds(stage_start))}},
               {{"name", "progressive_alignment"}});

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
  stage_start = Clock::now();
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
  const double final_folding_seconds = elapsed_seconds(stage_start);
  write_metric("stage",
               {{"seconds", json_number(final_folding_seconds)}},
               {{"name", "final_consensus_folding"}});

  size_t consensus_pairs = 0;
  for (uint i = 0; i < ss.size(); ++i)
    if (ss[i] != -1u && ss[i] > i)
      ++consensus_pairs;
  write_metric(
      "run_summary",
      {{"seconds", json_number(elapsed_seconds(run_start))},
       {"progressive_score", json_number(s)},
       {"alignment_columns", std::to_string(str.size())},
       {"consensus_pairs", std::to_string(consensus_pairs)},
       {"base_pair_probability_initial_nnz",
        std::to_string(base_pair_probability_nnz)},
       {"alignment_probability_initial_nnz",
        std::to_string(matching_probability_nnz)},
       {"base_pair_probability_final_nnz",
        std::to_string(final_base_pair_probability_nnz)},
       {"alignment_probability_final_nnz",
        std::to_string(final_matching_probability_nnz)},
       {"dd_calls", std::to_string(metrics_merge_id_)},
       {"dd_iterations", std::to_string(metrics_dd_iterations_)},
       {"dd_seconds", json_number(metrics_dd_seconds_)},
       {"cbp_peak", std::to_string(metrics_cbp_peak_)},
       {"cbp_added", std::to_string(metrics_cbp_added_)},
       {"cbp_priced", std::to_string(metrics_cbp_priced_)},
       {"cbp_removed", std::to_string(metrics_cbp_removed_)}});

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
  catch (const std::exception& e)
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
                        uint N1, uint N2, float min_th_s,
                        float pair_match_score) const {
    // Use the same logic as the original dafs.cpp:995-1001
    if (p_x.get(i, j) > CUTOFF && p_z.get(i, k) > CUTOFF &&
        p_y.get(k, l) > CUTOFF && p_z.get(j, l) > CUTOFF) {
        
        assert(p_x.get(i, j) <= 1.0);
        assert(p_y.get(k, l) <= 1.0);
        float p = (N1 * p_x.get(i, j) + N2 * p_y.get(k, l)) / (N1 + N2);
        float q = (p_z.get(i, k) + p_z.get(j, l)) / 2;
        // The probability terms occur once for each profile/endpoint in the
        // primal objective, whereas a pair-pair match is attached once to w.
        // Retain a candidate when its total unshifted contribution can be
        // positive.  A zero RIBOSUM weight preserves the historical candidate
        // set for exact backward compatibility.
        const float probability_score =
            w_ * (p - min_th_s) + (q - th_a_);
        if (w_ribosum_ == 0.0f)
          return p - min_th_s > 0.0f && probability_score > 0.0f;
        return 2.0f * probability_score + pair_match_score > 0.0f;
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
                                      const RibosumProfile& ribosum_x,
                                      const RibosumProfile& ribosum_y,
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
                         N1, N2, min_th_s,
                         w_ribosum_ * ribosum_x.pair_score(
                             i, j, ribosum_y, k, l))) {
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
                         N1, N2, min_th_s,
                         w_ribosum_ * ribosum_x.pair_score(
                             i, j, ribosum_y, k, l))) {
            add_cbp_if_new({{i, j}, {k, l}}, cbp, c_x, c_y, c_z);
        }
    }
}

void DAFS::generate_positive_reduced_cost_cbp(
    const GradientManager& gm,
    const GradientManager* alternate_gm,
    const SparseFloatMatrix& p_x,
    const SparseFloatMatrix& p_y,
    const SparseFloatMatrix& p_z,
    const VVU& p_z_forward, const VVU& p_z_reverse,
    const std::vector<std::pair<uint, uint>>& q_x_support,
    const std::vector<std::pair<uint, uint>>& q_y_support,
    uint N1, uint N2, float min_th_s,
    const RibosumProfile& ribosum_x,
    const RibosumProfile& ribosum_y,
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
        const float pair_score = w_ribosum_ *
            ribosum_x.pair_score(i, j, ribosum_y, k, l);
        if (!is_valid_cbp(i, j, k, l, p_x, p_y, p_z,
                          N1, N2, min_th_s, pair_score))
            return;
        const auto positive_reduced_cost = [&](const GradientManager& track) {
            return pair_score
                 + track.q_x(i, j) + track.q_y(k, l)
                 - track.q_z(i, k) - track.q_z(j, l) > 0.0f;
        };
        if (positive_reduced_cost(gm) ||
            (alternate_gm && positive_reduced_cost(*alternate_gm)))
            add_cbp_if_new(candidate, cbp, c_x, c_y, c_z);
    };

    // All positive intrinsic pair-score columns were inserted before the
    // iterations and are never pruned.  For every other missing column the
    // intrinsic score is non-positive; since q_z is non-negative, positive
    // reduced cost therefore implies q_x(i,j)>0 or q_y(k,l)>0.  Searching
    // both supports is exhaustive without a four-dimensional dense scan.
    for (const auto& [i, j] : q_x_support) {
        const bool active = gm.q_x(i, j) > 0.0f ||
            (alternate_gm && alternate_gm->q_x(i, j) > 0.0f);
        if (!active || p_x.get(i, j) <= CUTOFF)
            continue;
        for (const uint k : p_z_forward[i])
            for (const uint l : p_z_forward[j])
                if (k < l && p_y.get(k, l) > CUTOFF)
                    price_candidate(i, j, k, l);
    }

    for (const auto& [k, l] : q_y_support) {
        const bool active = gm.q_y(k, l) > 0.0f ||
            (alternate_gm && alternate_gm->q_y(k, l) > 0.0f);
        if (!active || p_y.get(k, l) <= CUTOFF)
            continue;
        for (const uint i : p_z_reverse[k])
            for (const uint j : p_z_reverse[l])
                if (i < j && p_x.get(i, j) > CUTOFF)
                    price_candidate(i, j, k, l);
    }
}
