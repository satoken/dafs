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
//#include <fstream>
#include "fa.h"
#include "fold.h"
#include "nussinov.h"
#include "ipknot.h"
#include "align.h"
#include "needleman_wunsch.h"
#include "alifold.h"
#include "ip.h"
#include "typedefs.h"

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

#define FOREACH(i, v) for (auto i = std::begin(v); i != std::end(v); ++i)
#define CUTOFF 0.01

#define SPARSE_UPDATE
//#define ADAGRAD
//#define ADAM

class DAFS
{
private:
  // nodes in the guide tree
  typedef std::pair<float, std::pair<uint, uint>> node_t;
  // indices for consensus base pairs
  typedef std::pair<std::pair<uint, uint>, std::pair<uint, uint>> CBP;

public:
  DAFS()
      : use_alifold_(false),
        use_alifold1_(true)
  {
  }

  ~DAFS()
  {
  }

  DAFS &parse_options(int &argc, char **&argv);
  int run();

private:
  void relax_matching_probability();
  void relax_basepairing_probability();
  void relax_fourway_consistency();
  void build_tree();
  void print_tree(std::ostream &os, int i) const;
  void project_alignment(ALN &aln, const ALN &aln1, const ALN &aln2, const VU &z) const;
  void project_secondary_structure(VU &xx, VU &yy, const VU &x, const VU &y, const VU &z) const;
  void average_matching_probability(VVF &posterior, const ALN &aln1, const ALN &aln2) const;
  void average_basepairing_probability(VVF &posterior, const ALN &aln, bool use_alifold) const;
  void update_basepairing_probability(VVF &posterior, const VU &ss, const std::string &str,
                                      const ALN &aln, bool use_alifold) const;
  void align_alignments(ALN &aln, const ALN &aln1, const ALN &aln2) const;
  float align_alignments(VU &ss, ALN &aln, const ALN &aln1, const ALN &aln2) const;
  float solve(VU &x, VU &y, VU &z, const VVF &p_x, const VVF &p_y, const VVF &p_z,
              const ALN &aln1, const ALN &aln2) const
  {
#if defined(WITH_GLPK) || defined(WITH_CPLEX) || defined(WITH_GUROBI)
    return t_max_ != 0 ? solve_by_dd(x, y, z, p_x, p_y, p_z, aln1, aln2) : solve_by_ip(x, y, z, p_x, p_y, p_z, aln1, aln2);
#else
    return solve_by_dd(x, y, z, p_x, p_y, p_z, aln1, aln2);
#endif
  }
  float solve_by_dd(VU &x, VU &y, VU &z, const VVF &p_x, const VVF &p_y, const VVF &p_z,
                    const ALN &aln1, const ALN &aln2) const;
  float solve_by_ip(VU &x, VU &y, VU &z, const VVF &p_x, const VVF &p_y, const VVF &p_z,
                    const ALN &aln1, const ALN &aln2) const;
  void align(ALN &aln, int ch) const;
  float align(VU &ss, ALN &aln, int ch) const;
  float refine(VU &ss, ALN &aln) const;
  void output_verbose(const VU &x, const VU &y, const VU &z, const ALN &aln1, const ALN &aln2) const;
  void output(std::ostream &os, const ALN &aln) const;
  void output(std::ostream &os, ALN::const_iterator b, ALN::const_iterator e) const;

private:
  float w_pct_a_;                   // the weight of PCT for alignment matching probabilities
  float w_pct_s_;                   // the weight of PCT for base-pairing probabilities
  float w_pct_f_;                   // the weight of four-way PCT
  uint n_refinement_;               // the number of the iterative refinement
  uint t_max_;                      // the maximum number of the iteration of the subgradient update
  float th_a_;                      // the threshold for base-pairing probabilities
  VF th_s_;                         // the threshold for alignment matching probabilities
  float w_;                         // the weight for base pairs in the objective function
  float eta0_;                      // the initial step width of the subgradient update
  std::unique_ptr<Align::Model> a_model_;           // alignment model
  std::unique_ptr<Align::Decoder> a_decoder_;       // alignment decoder
  std::unique_ptr<Fold::Model> s_model_;            // folding model
  std::unique_ptr<Fold::Decoder> s_decoder_;        // folding decoder
  std::unique_ptr<Fold::Decoder> s_decoder1_;       // folding decoder for the final folding
  std::vector<Fasta> fa_;           // input sequences
  std::vector<std::vector<MP>> mp_; // alignment matching probability matrices
  std::vector<BP> bp_;              // base-pairing probability matrices
  VVF sim_;                         // simalarity matrix between input sequences
  std::vector<node_t> tree_;        // guide tree
  bool use_alifold_;
  bool use_alifold1_;
  bool use_bp_update_;
  bool use_bp_update1_;
  // bool use_bpscore_;
  uint verbose_;
};

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
    FOREACH(v, mp[i])
    os << " " << v->first << ":" << v->second;
    os << std::endl;
  }
}

static void
print_bp(std::ostream &os, const BP &bp)
{
  for (uint i = 0; i != bp.size(); ++i)
  {
    os << i << ":";
    FOREACH(v, bp[i])
    os << " " << v->first << ":" << v->second;
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
      FOREACH(j, bp[x][i])
      os << " " << j->first << ":" << j->second;
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
        FOREACH(k, mp[x][y][i])
        os << " " << k->first << ":" << k->second;
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
      VVF posterior(L1, VF(L2, 0.0));
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
              posterior[i][j] += p_ik * p_jk * w;
            }
          }
        }
      }

      mp[x][y].resize(L1);
      for (uint i = 0; i != L1; ++i)
      {
        for (uint j = 0; j != L2; ++j)
        {
          float v = posterior[i][j] / sum_w;
          if (v > CUTOFF)
            mp[x][y][i].push_back(std::make_pair(j, v));
        }
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
    VVF p(L1, VF(L1, 0.0));

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
                p[i][j] += p_kl * p_ik * p_jl * w;
            }
          }
        }
      }
    }

    bp[x].resize(L1);
    for (uint i = 0; i != L1 - 1; ++i)
      for (uint j = i + 1; j != L1; ++j)
      {
        float v = p[i][j] / sum_w;
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
      VVF posterior(L1, VF(L2, 0.0));
      assert(L1 == mp_[x][y].size());

      for (uint i = 0; i != L1; ++i)
      {
        for (const auto [k, p_ik] : mp_[x][y][i])
        {
          posterior[i][k] += p_ik * (1.0 - w_pct_f_);
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
                posterior[i][k] += p_ij * p_kl * p_jl * w_pct_f_;
                posterior[j][l] += p_ij * p_kl * p_ik * w_pct_f_;
                ++ll1;
                ++ll2;
              }
            }
          }
        }
      }

      mp[x][y].resize(L1);
      for (uint i = 0; i != L1; ++i)
      {
        for (uint j = 0; j != L2; ++j)
        {
          float v = posterior[i][j];
          if (v > CUTOFF)
            mp[x][y][i].push_back(std::make_pair(j, v));
        }
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
    average_matching_probability(VVF &posterior, const ALN &aln1, const ALN &aln2) const
{
  const uint L1 = aln1.front().second.size();
  const uint L2 = aln2.front().second.size();
  const uint N1 = aln1.size();
  const uint N2 = aln2.size();
  VVF p(L1, VF(L2, 0.0));
  FOREACH(it1, aln1)
  {
    assert(L1 == it1->second.size());
    FOREACH(it2, aln2)
    {
      assert(L2 == it2->second.size());
      const MP &m = mp_[it1->first][it2->first];
      for (uint i = 0, ii = 0; i != L1; ++i)
      {
        if (!it1->second[i])
          continue;
        assert(ii < m.size());
        SV::const_iterator x = m[ii].begin();
        for (uint j = 0, jj = 0; j != L2 && x != m[ii].end(); ++j)
        {
          if (!it2->second[j])
            continue;
          if (jj == x->first)
          {
            p[i][j] += x->second / (N1 * N2);
            ++x;
          }
          ++jj;
        }
        ++ii;
      }
    }
  }
  for (uint i = 0; i != p.size(); ++i)
    for (uint j = 0; j != p[i].size(); ++j)
    {
      if (p[i][j] <= CUTOFF)
        p[i][j] = 0.0;
      if (p[i][j] > 1.0)
        p[i][j] = 1.0;
      //assert(p[i][j]<=1.0);
    }
  std::swap(posterior, p);
}

void DAFS::
    average_basepairing_probability(VVF &posterior, const ALN &aln, bool use_alifold) const
{
  // calculate an averaged base-pairing probabilities
  const uint L = aln.front().second.size();
  const uint N = aln.size();
  VVF p(L, VF(L, 0.0));
  FOREACH(it, aln)
  {
    assert(L == it->second.size());
    uint s = it->first;
    VU idx(fa_[s].size());
    for (uint i = 0, j = 0; i != L; ++i)
      if (it->second[i])
        idx[j++] = i;
    const BP &bp = bp_[s];
    for (uint i = 0; i != bp.size(); ++i)
      for (uint j = 0; j != bp[i].size(); ++j)
        p[idx[i]][idx[bp[i][j].first]] += bp[i][j].second / N;
  }

  // mix base-pairing probabilities by alifold and averaged base-pairing probabilities
  if (use_alifold)
  {
    BP bp;
    Alifold ali(0.0 /*CUTOFF*/);
    ali.fold(aln, fa_, bp);
    assert(L == bp.size());
    for (uint i = 0; i != bp.size(); ++i)
      for (uint j = 0; j != bp[i].size(); ++j)
        p[i][bp[i][j].first] += bp[i][j].second;

    for (uint i = 0; i != L - 1; ++i)
      for (uint j = i + 1; j != L; ++j)
        p[i][j] /= 2;
  }

  // cut off
  for (uint i = 0; i != L - 1; ++i)
    for (uint j = i + 1; j != L; ++j)
    {
      if (p[i][j] <= CUTOFF)
        p[i][j] = 0.0;
      assert(p[i][j] <= 1.0);
    }
  std::swap(p, posterior);
}

void DAFS::
    update_basepairing_probability(VVF &posterior, const VU &ss, const std::string &str,
                                   const ALN &aln, bool use_alifold) const
{
  const uint L = aln.front().second.size();
  const uint N = aln.size();
  const uint plevel = th_s_.size();
  VVF p(L, VF(L, 0.0));

  // calculate an averaged base-pairing probabilities
  //   which are constrained by the previous prediction
  FOREACH(it, aln)
  {
    assert(L == it->second.size());
    // calculate the mapping
    uint s = it->first;
    VU idx(fa_[s].size()); // from the sequence to the alignment
    VU rev(L, -1u);        // from the alignment to the sequence
    for (uint i = 0, j = 0; i != L; ++i)
      if (it->second[i])
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
          p[idx[i]][idx[bp[i][j].first]] += bp[i][j].second / N;
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
      Alifold ali(0.0 /*CUTOFF*/);
      ali.fold(aln, fa_, con, bp);
      assert(L == bp.size());
      for (uint i = 0; i != bp.size(); ++i)
        for (uint j = 0; j != bp[i].size(); ++j)
          p[i][bp[i][j].first] += bp[i][j].second;
    }

    for (uint i = 0; i != L - 1; ++i)
      for (uint j = i + 1; j != L; ++j)
        p[i][j] /= 2;
  }

  // cut off
  for (uint i = 0; i != L - 1; ++i)
    for (uint j = i + 1; j != L; ++j)
    {
      if (p[i][j] <= CUTOFF)
        p[i][j] = 0.0;
      assert(p[i][j] <= 1.0);
    }
  std::swap(p, posterior);
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
    FOREACH(jj, mp[i - 1])
    {
      for (; j - 1 < jj->first; ++j)
      {
        dp[i][j] = dp[i][j - 1];
        tr[i][j] = tr[i][j - 1] + 1;
        if (dp[i][j] < dp[i - 1][j])
        {
          dp[i][j] = dp[i - 1][j];
          tr[i][j] = tr[i - 1][j] + 1;
        }
      }

      dp[i][j] = dp[i - 1][j - 1] + jj->second;
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
  FOREACH(q, aln1)
  {
    p->first = q->first;
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
        p->second[r++] = q->second[i];
        ++k;
      }
      else
        p->second[r++] = q->second[i];
    }
    while (k < L2)
    {
      p->second[r++] = false;
      k++;
    }
    ++p;
  }
  FOREACH(q, aln2)
  {
    p->first = q->first;
    p->second.resize(L, false);
    uint k = 0, r = 0;
    for (uint i = 0; i != z.size(); ++i)
    {
      if (z[i] != -1u)
      {
        while (k < z[i])
          p->second[r++] = q->second[k++];
        p->second[r++] = q->second[k++];
      }
      else
        p->second[r++] = false;
    }
    while (k < L2)
      p->second[r++] = q->second[k++];
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
    align_alignments(ALN &aln, const ALN &aln1, const ALN &aln2) const
{
  // calculate posteriors
  VVF p_x, p_y, p_z;
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
    align_alignments(VU &ss, ALN &aln, const ALN &aln1, const ALN &aln2) const
{
  // calculate posteriors
  VVF p_x, p_y, p_z;
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

#ifdef ADAGRAD
float adagrad_update(float &g2, const float g, const float eta0)
{
  const float eps = 1e-6;
  g2 += g * g;
  return eta0 * g / std::sqrt(g2 + eps);
}
#endif

#ifdef ADAM
float adam_update(int t, float &m, float &v, const float g, float alpha = 0.1)
{
  const float beta1 = 0.9;
  const float beta2 = 0.999;
  const float eps = 1e-8;
  m = beta1 * m + (1 - beta1) * g;
  v = beta2 * v + (1 - beta2) * g * g;
  const float m_hat = m / (1 - std::pow(beta1, t));
  const float v_hat = v / (1 - std::pow(beta2, t));
  return alpha * m_hat / (std::sqrt(v_hat) + eps);
}
#endif

float DAFS::
    solve_by_dd(VU &x, VU &y, VU &z,
                const VVF &p_x, const VVF &p_y, const VVF &p_z,
                const ALN &aln1, const ALN &aln2) const
{
  const uint L1 = p_x.size();
  const uint L2 = p_y.size();
  const uint N1 = aln1.size();
  const uint N2 = aln2.size();

  // enumerate the candidates of consensus base-pairs
  std::vector<CBP> cbp; // consensus base-pairs
  float min_th_s = *std::min_element(th_s_.begin(), th_s_.end());
#ifdef SPARSE_UPDATE
  VVU c_x(L1), c_y(L2), c_z(L1); // project consensus base-pairs into each structure and alignment
#endif
  for (uint i = 0; i != L1 - 1; ++i)
    for (uint j = i + 1; j != L1; ++j)
      if (p_x[i][j] > CUTOFF)
        for (uint k = 0; k != L2 - 1; ++k)
          if (p_z[i][k] > CUTOFF)
            for (uint l = k + 1; l != L2; ++l)
              if (p_y[k][l] > CUTOFF && p_z[j][l] > CUTOFF)
              {
                assert(p_x[i][j] <= 1.0);
                assert(p_y[k][l] <= 1.0);
                float p = (N1 * p_x[i][j] + N2 * p_y[k][l]) / (N1 + N2);
                float q = (p_z[i][k] + p_z[j][l]) / 2;
                if (p - min_th_s > 0.0 && w_ * (p - min_th_s) + (q - th_a_) > 0.0)
                {
                  cbp.push_back(std::make_pair(std::make_pair(i, j), std::make_pair(k, l)));
#ifdef SPARSE_UPDATE
                  c_x[i].push_back(j);
                  c_y[k].push_back(l);
                  c_z[i].push_back(k);
                  c_z[j].push_back(l);
#endif
                }
              }
#ifdef SPARSE_UPDATE
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
#endif

  // precalculate the range for alignment, i.e. alignment envelope
  a_decoder_->initialize(p_z);

  // multipliers
  VVF q_x(L1, VF(L1, 0.0));
  VVF q_y(L2, VF(L2, 0.0));
  VVF q_z(L1, VF(L2, 0.0));

  //uint c=0;
  float c = 0.0;
#if defined ADAGRAD
  VVF g_x(L1, VF(L1, 0.0));
  VVF g_y(L2, VF(L2, 0.0));
  VVF g_z(L1, VF(L2, 0.0));
#elif defined ADAM
  VVF m_x(L1, VF(L1, 0.0)), v_x(L1, VF(L1, 0.0));
  VVF m_y(L2, VF(L2, 0.0)), v_y(L2, VF(L2, 0.0));
  VVF m_z(L1, VF(L2, 0.0)), v_z(L1, VF(L2, 0.0));
#else
  float eta = eta0_;
#endif
  float s_prev = 0.0;
  uint violated = 0;
  uint t;
  for (t = 0; t != t_max_; ++t)
  {
    // solve the subproblems
    float s = 0.0;
    s += s_decoder_->decode(w_ * 2 * N1 / (N1 + N2), p_x, q_x, x);
    s += s_decoder_->decode(w_ * 2 * N2 / (N1 + N2), p_y, q_y, y);
    s += a_decoder_->decode(p_z, q_z, z);

    if (verbose_ >= 2)
      output_verbose(x, y, z, aln1, aln2);

    // update the multipliers
    violated = 0;
    VVI t_x(L1, VI(L1, 0));
    VVI t_y(L2, VI(L2, 0));
    VVI t_z(L1, VI(L2, 0));
    for (uint u = 0; u != cbp.size(); ++u)
    {
      const uint i = cbp[u].first.first, j = cbp[u].first.second;
      const uint k = cbp[u].second.first, l = cbp[u].second.second;
      const float s_w = q_x[i][j] + q_y[k][l] - q_z[i][k] - q_z[j][l];
      const int w_ijkl = s_w > 0.0f ? 1 : 0;
      if (w_ijkl)
      {
        s += s_w;    /* * w_ijkl*/
        t_x[i][j]++; // += w_ijkl;
        t_y[k][l]++; // += w_ijkl;
        t_z[i][k]++; // += w_ijkl;
        t_z[j][l]++; // += w_ijkl;
      }
    }

    // update Lagrangian for x (=q_x)
#ifdef SPARSE_UPDATE // efficient implementation using sparsity
    for (uint i = 0; i != L1; ++i)
    {
      const uint j = x[i];
      if (j != -1u && t_x[i][j] != 1)
      {
        violated++;
#if defined ADAGRAD
        q_x[i][j] -= adagrad_update(g_x[i][j], t_x[i][j] - 1, eta0_);
#elif defined ADAM
        q_x[i][j] -= adam_update(t + 1, m_x[i][j], v_x[i][j], t_x[i][j] - 1, eta0_);
#else
        q_x[i][j] -= eta * (t_x[i][j] - 1);
#endif
      }
      for (uint jj = 0; jj != c_x[i].size(); ++jj)
      {
        const uint j = c_x[i][jj];
        if (x[i] != j && t_x[i][j] != 0)
        {
          violated++;
#if defined ADAGRAD
          q_x[i][j] -= adagrad_update(g_x[i][j], t_x[i][j], eta0_);
#elif defined ADAM
          q_x[i][j] -= adam_update(t + 1, m_x[i][j], v_x[i][j], t_x[i][j], eta0_);
#else
          q_x[i][j] -= eta * t_x[i][j];
#endif
        }
      }
    }
#else // naive implementation
    for (uint i = 0; i != L1 - 1; ++i)
      for (uint j = i + 1; j != L1; ++j)
      {
        const int x_ij = x[i] == j ? 1 : 0;
        if (t_x[i][j] - x_ij != 0)
        {
          violated++;
#if defined ADAGRAD
          q_x[i][j] -= adagrad_update(g_x[i][j], t_x[i][j] - x_ij, eta0_);
#elif defined ADAM
          q_x[i][j] -= adam_update(t + 1, m_x[i][j], v_x[i][j], t_x[i][j] - x_ij, eta0_);
#else
          q_x[i][j] -= eta * (t_x[i][j] - x_ij);
#endif
        }
      }
#endif

    // update Lagrangian for y (=q_y)
#ifdef SPARSE_UPDATE
    for (uint k = 0; k != L2; ++k)
    {
      const uint l = y[k];
      if (l != -1u && t_y[k][l] != 1)
      {
        violated++;
#if defined ADAGRAD
        q_y[k][l] -= adagrad_update(g_y[k][l], t_y[k][l] - 1, eta0_);
#elif defined ADAM
        q_y[k][l] -= adam_update(t + 1, m_y[k][l], v_y[k][l], t_y[k][l] - 1, eta0_);
#else
        q_y[k][l] -= eta * (t_y[k][l] - 1);
#endif
      }
      for (uint ll = 0; ll != c_y[k].size(); ++ll)
      {
        const uint l = c_y[k][ll];
        if (y[k] != l && t_y[k][l] != 0)
        {
          violated++;
#if defined ADAGRAD
          q_y[k][l] -= adagrad_update(g_y[k][l], t_y[k][l], eta0_);
#elif defined ADAM
          q_y[k][l] -= adam_update(t + 1, m_y[k][l], v_y[k][l], t_y[k][l], eta0_);
#else
          q_y[k][l] -= eta * t_y[k][l];
#endif
        }
      }
    }
#else // naive implementation
    for (uint k = 0; k != L2 - 1; ++k)
      for (uint l = k + 1; l != L2; ++l)
      {
        const int y_kl = y[k] == l ? 1 : 0;
        if (t_y[k][l] - y_kl != 0)
        {
          violated++;
#if defined ADAGRAD
          q_y[k][l] -= adagrad_update(g_y[k][l], t_y[k][l] - y_kl, eta0_);
#elif defined ADAM
          q_y[k][l] -= adam_update(t + 1, m_y[k][l], v_y[k][l], t_y[k][l] - y_kl, eta0_);
#else
          q_y[k][l] -= eta * (t_y[k][l] - y_kl);
#endif
        }
      }
#endif

    // update Lagrangian for z (=q_z)
#ifdef SPARSE_UPDATE
    for (uint i = 0; i != L1; ++i)
    {
      const uint k = z[i];
      if (k != -1u) // z_ik==1
      {
        if (t_z[i][k] > 1)
          violated++;
#if defined ADAGRAD
        q_z[i][k] = std::max(0.0f, q_z[i][k] - adagrad_update(g_z[i][k], 1 - t_z[i][k], eta0_));
#elif defined ADAM
        q_z[i][k] = std::max(0.0f, q_z[i][k] - adam_update(t + 1, m_z[i][k], v_z[i][k], 1 - t_z[i][k], eta0_));
#else
        q_z[i][k] = std::max(0.0f, q_z[i][k] - eta * (1 - t_z[i][k]));
#endif
      }
      for (uint kk = 0; kk != c_z[i].size(); ++kk)
      {
        const uint k = c_z[i][kk];
        if (z[i] != k) // z_ik==0
        {
          if (t_z[i][k] > 0)
            violated++;
#if defined ADAGRAD
          q_z[i][k] = std::max(0.0f, q_z[i][k] - adagrad_update(g_z[i][k], -t_z[i][k], eta0_));
#elif defined ADAM
          q_z[i][k] = std::max(0.0f, q_z[i][k] - adam_update(t + 1, m_z[i][k], v_z[i][k], -t_z[i][k], eta0_));
#else
          q_z[i][k] = std::max(0.0f, q_z[i][k] + eta * t_z[i][k]);
#endif
        }
      }
    }
#else // naive implementation
    for (uint i = 0; i != L1; ++i)
      for (uint k = 0; k != L2; ++k)
      {
        const int z_ik = z[i] == k ? 1 : 0;
        if (z_ik - t_z[i][k] < 0)
          violated++;
#if defined ADAGRAD
        q_z[i][k] = std::max(0.0f, q_z[i][k] - adagrad_update(g_z[i][k], z_ik - t_z[i][k], eta0_));
#elif defined ADAM
        q_z[i][k] = std::max(0.0f, q_z[i][k] - adam_update(t + 1, m_z[i][k], v_z[i][k], z_ik - t_z[i][k], eta0_));
#else
        q_z[i][k] = std::max(0.0f, q_z[i][k] - eta * (z_ik - t_z[i][k]));
#endif
      }
#endif

#if !defined ADAGRAD && !defined ADAM
    spdlog::debug("Step: {}, eta: {}, L: {}, Violated: {}", t, eta, s, violated);
#else
    spdlog::debug("Step: {}, L: {}, Violated: {}", t, s, violated);
#endif

    if (violated == 0)
      break; // all constraints were satisfied.

      // update the step width
#if !defined ADAGRAD && !defined ADAM
    if (s > s_prev || t == 0)
    {
      c += std::max(0.0f, 4.0f * cbp.size() - violated) / (4.0 * cbp.size());
      //eta = eta0_/(1.0+sqrt(c));
      eta = eta0_ / (1.0 + c);
    }
#endif
    s_prev = s;
  }
  spdlog::info("Step: {}, Violated: {}", t, violated);

  return s_prev;
}

float DAFS::
    solve_by_ip(VU &x, VU &y, VU &z,
                const VVF &p_x, const VVF &p_y, const VVF &p_z,
                const ALN &aln1, const ALN &aln2) const
{
  const uint L1 = p_x.size();
  const uint L2 = p_y.size();

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
      if (p_z[i][k] > CUTOFF)
        v_z[i][k] = ip.make_variable(p_z[i][k] - th_a_);

  // enumerate the candidates of consensus base-pairs
  std::vector<CBP> cbp;
  for (uint i = 0; i != L1 - 1; ++i)
    for (uint j = i + 1; j != L1; ++j)
      if (p_x[i][j] > CUTOFF)
        for (uint k = 0; k != L2 - 1; ++k)
          if (p_z[i][k] > CUTOFF)
            for (uint l = k + 1; l != L2; ++l)
              if (p_y[k][l] > CUTOFF && p_z[j][l] > CUTOFF)
              {
                assert(p_x[i][j] <= 1.0);
                assert(p_y[k][l] <= 1.0);
                float p = (p_x[i][j] + p_y[k][l]) / 2;
                float q = (p_z[i][k] + p_z[j][l]) / 2;
                if (p - min_th_s > 0.0 && w_ * (p - min_th_s) + (q - th_a_) > 0.0)
                {
                  cbp.push_back(std::make_pair(std::make_pair(i, j), std::make_pair(k, l)));
                  v_w.push_back(ip.make_variable(0.0));
                  if (v_x[i][j] < 0)
                    v_x[i][j] = ip.make_variable(w_ * (p_x[i][j] - min_th_s));
                  if (v_y[k][l] < 0)
                    v_y[k][l] = ip.make_variable(w_ * (p_y[k][l] - min_th_s));
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
    align(ALN &aln, int ch) const
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
    align(VU &ss, ALN &aln, int ch) const
{
  float s = 0.0;
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
    s = align_alignments(ss, aln, aln1, aln2);
  }
  return s;
}

float DAFS::
    refine(VU &ss, ALN &aln) const
{
  VU group[2];
  do
  {
    group[0].clear();
    group[1].clear();
    for (uint i = 0; i != aln.size(); ++i)
      group[rand() % 2].push_back(i);
  } while (group[0].empty() || group[1].empty());

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
    ("w,weight", "Weight of the expected accuracy score for secondary structures", cxxopts::value<float>()->default_value("4.0"))
    ("eta", "Initial step width for the subgradient optimization", cxxopts::value<float>()->default_value("0.5"))
    ("m,max-iter", "The maximum number of iteration of the subgradient optimization", cxxopts::value<int>()->default_value("600"), "T")
    ("f,fourway-pct", "Weight of four-way PCT", cxxopts::value<float>()->default_value("0.0"))
    ("v,verbose", "The level of verbose outputs", cxxopts::value<int>()->default_value("0"));

  options.add_options("Aligning")
    ("a,align-model", "Alignment model for calcualating matching probablities (value=CONTRAlign, ProbCons)", 
      cxxopts::value<std::string>()->default_value("ProbCons"))
    ("p,align-pct", "Weight of PCT for matching probabilities", cxxopts::value<float>()->default_value("0.25"))
    ("u,align-th", "Threshold for matching probabilities", cxxopts::value<float>()->default_value("0.01"))
    ("align-aux", "Load matching probability matrices from FILENAME", cxxopts::value<std::string>(), "FILENAME");

  options.add_options("Folding")
    ("s,fold-model", "Folding model for calculating base-pairing probablities (value=Boltzmann, Vienna, CONTRAfold)",
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
    ("fold-aux", "Load base-pairing probability matrices from FILENAME", cxxopts::value<std::string>(), "FILENAME");

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
    w_ = res["weight"].as<float>();
    eta0_ = res["eta"].as<float>();
    t_max_ = res["max-iter"].as<int>();
    w_pct_f_ = res["fourway-pct"].as<float>();
    verbose_ = res["verbose"].as<int>();
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

    if (res["align-aux"].count())
      a_model_ = std::make_unique<AUXAlign>(res["align-aux"].as<std::string>(), CUTOFF);
    else if (res["align-model"].as<std::string>() == "CONTRAlign")
      a_model_ = std::make_unique<CONTRAlign>(th_a_);
    else if (res["align-model"].as<std::string>() == "ProbCons")
      a_model_ = std::make_unique<ProbCons>(th_a_);
    else
      throw "Unknown alignment model: " + res["align-model"].as<std::string>();
    assert(a_model_);
    a_decoder_ = std::make_unique<SparseNeedlemanWunsch>(th_a_);

    // options for folding
    w_pct_s_ = res["fold-pct"].as<float>();
    use_alifold_ = res["no-alifold"].count() == 0;
    if (res["fold-aux"].count())
      s_model_ = std::make_unique<AUXFold>(res["fold-aux"].as<std::string>(), CUTOFF);
    else if (res["fold-model"].as<std::string>() == "Boltzmann")
      s_model_ = std::make_unique<RNAfold>(true, nullptr, CUTOFF);
    else if (res["fold-model"].as<std::string>() == "Vienna")
      s_model_ = std::make_unique<RNAfold>(false, nullptr, CUTOFF);
    else if (res["fold-model"].as<std::string>() == "CONTRAfold")
      s_model_ = std::make_unique<CONTRAfold>(CUTOFF);
    else
      throw "Unknown folding model: " + res["fold-model"].as<std::string>();
    assert(s_model_);

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

    if (res["fold-decoder"].as<std::string>() == "IPknot" || res["ipknot"].count())
    {
      s_decoder_ = std::make_unique<IPknot>(th_s_);
      s_decoder1_ = std::make_unique<IPknot>(th_s1);
    }
    else if (res["fold-decoder"].as<std::string>() == "Nussinov")
    {
      s_decoder_ = std::make_unique<SparseNussinov>(th_s_[0]);
      s_decoder1_ = std::make_unique<SparseNussinov>(th_s1[0]);
    }
    else
      throw "Unknown folding decoder: " + res["fold-decoder"].as<std::string>();
    assert(s_decoder_);

    use_bp_update_ = res["bp-update"].count() > 0;
    use_bp_update1_ = res["bp-update1"].count() > 0 ^ res["ipknot"].count() > 0;

    // read sequences
    Fasta::load(fa_, res["input"].as<std::string>().c_str());
  }
  catch (cxxopts::option_has_no_value_exception e)
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
  for (uint i = 0; i != N; ++i)
  {
    sim_[i][i] = 1.0;
    for (uint j = i + 1; j != N; ++j)
      sim_[i][j] = sim_[j][i] = calculate_similarity_score(mp_[i][j], fa_[i].size(), fa_[j].size());
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

  std::string str;
  if (s_decoder1_)
  {
    // compute the common secondary structures from the averaged base-pairing matrix
    VVF p;
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
