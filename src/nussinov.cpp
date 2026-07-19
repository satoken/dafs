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
#include "linearfold/fold/linfold.h"
#include "linearfold/param/pair_objective.h"
#include <cassert>
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

template <typename Probability, typename Multiplier>
float
LinearNussinov::
decode_impl(float w, const Probability& probability, uint L,
            Multiplier multiplier, VU& ss)
{
  if (L == 0) {
    ss.clear();
    return 0.0f;
  }

  const std::string dummy_sequence(L, 'a');
  auto params = std::make_unique<PairObjectiveNearestNeighbor>(dummy_sequence);
  LinFold<PairObjectiveNearestNeighbor> decoder(std::move(params));
  LinFold<PairObjectiveNearestNeighbor>::Options options;
  options.beam_size(beam_size_)
         .min_hairpin_loop_length(2)
         .set_allowed_pair('a', 'a')
         .pair_score([&](u_int32_t i, u_int32_t j) {
           return w * (probability(i-1, j-1) - th_) - multiplier(i-1, j-1);
         });

  decoder.compute_viterbi(dummy_sequence, options);
  const std::vector<u_int32_t> pairs = decoder.traceback_viterbi();
  ss.assign(L, -1u);
  float score = 0.0f;
  for (uint i = 1; i <= L; ++i) {
    const uint j = pairs[i];
    if (j > i) {
      ss[i-1] = j-1;
      score += w * (probability(i-1, j-1) - th_) - multiplier(i-1, j-1);
    }
  }
  return score;
}

float
LinearNussinov::
decode(float w, const VVF& p, const VVF& q, VU& ss)
{
  return decode_impl(w, [&](uint i, uint j) { return p[i][j]; }, p.size(),
                     [&](uint i, uint j) { return q[i][j]; }, ss);
}

float
LinearNussinov::
decode(float w, const VVF& p, const SparseFloatMatrix& q, VU& ss)
{
  return decode_impl(w, [&](uint i, uint j) { return p[i][j]; }, p.size(),
                     [&](uint i, uint j) { return q.get(i, j); }, ss);
}

float
LinearNussinov::
decode(float w, const SparseFloatMatrix& p, const VVF& q, VU& ss)
{
  return decode_impl(w, [&](uint i, uint j) { return p.get(i, j); }, p.rows(),
                     [&](uint i, uint j) { return q[i][j]; }, ss);
}

float
LinearNussinov::
decode(float w, const SparseFloatMatrix& p, const SparseFloatMatrix& q, VU& ss)
{
  return decode_impl(w, [&](uint i, uint j) { return p.get(i, j); }, p.rows(),
                     [&](uint i, uint j) { return q.get(i, j); }, ss);
}

float
LinearNussinov::
decode(const VVF& p, VU& ss, std::string& str)
{
  const auto zero = [](uint, uint) { return 0.0f; };
  const float score = decode_impl(1.0f,
                                  [&](uint i, uint j) { return p[i][j]; },
                                  p.size(), zero, ss);
  make_brackets(ss, str);
  return score;
}

float
LinearNussinov::
decode(const SparseFloatMatrix& p, VU& ss, std::string& str)
{
  const auto zero = [](uint, uint) { return 0.0f; };
  const float score = decode_impl(1.0f,
                                  [&](uint i, uint j) { return p.get(i, j); },
                                  p.rows(), zero, ss);
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
