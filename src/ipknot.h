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

#ifndef __INC_IPKNOT_H__
#define __INC_IPKNOT_H__

#include "fold.h"
#include "ip.h"
#include "nussinov.h"

class IPknot : public Fold::Decoder
{
public:
  IPknot(const VF& th, int n_th = 1);
  float decode(float w, const VVF& p, const VVF& q, VU& ss);
  float decode(float w, const SparseFloatMatrix& p, const VVF& q,
               VU& ss) override;
  float decode(float w, const SparseFloatMatrix& p,
               const SparseFloatMatrix& q, VU& ss) override;
  float decode(const VVF& p, VU& ss, std::string& str);
  float decode(const SparseFloatMatrix& p, VU& ss,
               std::string& str) override;
  void make_brackets(const VU& ss, std::string& str) const;

private:
  void make_objective(IP& ip, float w, const VVF& p, const VVF& q);
  void make_objective(IP& ip, float w, const SparseFloatMatrix& p,
                      const VVF& q);
  void make_objective(IP& ip, float w, const SparseFloatMatrix& p,
                      const SparseFloatMatrix& q);
  void make_objective(IP& ip, const VVF& p);
  void make_objective(IP& ip, const SparseFloatMatrix& p);
  void reset_candidates(uint length);
  void add_candidate(IP& ip, uint level, uint i, uint j, float score);
  void make_constraints(IP& ip);
  float solve(IP& ip, VU& ss);
  
private:
  VF th_;
  VF alpha_;
  bool levelwise_;
  bool stacking_constraints_;
  using SparseVariables = std::vector<std::vector<std::pair<uint, int>>>;
  std::vector<SparseVariables> v_left_, v_right_;
  uint n_th_;
  VU plevel_;
};

// Linear-time, fixed-beam IPknot approximation for sparse LinearPartition
// supports. Each threshold level is decoded with a linear Nussinov recurrence;
// later levels receive only pairs crossing an already selected lower-level
// pair. This preserves IPknot's layered output without invoking a MIP solver.
class LinearIPknot : public Fold::Decoder
{
public:
  LinearIPknot(const VF& th, uint beam_size);

  float decode(float w, const VVF& p, const VVF& q, VU& ss) override;
  float decode(float w, const SparseFloatMatrix& p, const VVF& q,
               VU& ss) override;
  float decode(float w, const SparseFloatMatrix& p,
               const SparseFloatMatrix& q, VU& ss) override;
  float decode(const VVF& p, VU& ss, std::string& str) override;
  float decode(const SparseFloatMatrix& p, VU& ss,
               std::string& str) override;
  void make_brackets(const VU& ss, std::string& str) const override;

private:
  float decode_sparse(float w, const SparseFloatMatrix& p,
                      const SparseFloatMatrix& q, VU& ss);
  static void dense_to_sparse(const VVF& dense, SparseFloatMatrix& sparse);
  static void linear_decompose_plevel(const VU& ss, VU& plevel);
  std::string brackets_for(const VU& ss) const;

  VF th_;
  uint beam_size_;
  VU plevel_;
  VU last_structure_;
};

#endif  //  __INC_IPKNOT_H__

// Local Variables:
// mode: C++
// End:
