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

#ifndef __INC_NUSSINOV_H__
#define __INC_NUSSINOV_H__

#include <cstddef>

#include "fold.h"

struct LinearNussinovResult
{
  float score = 0.0f;
  float upper_bound = 0.0f;
  float pruned_upper_bound = 0.0f;
  float additive_upper_bound = 0.0f;
  size_t pruned_states = 0;
};

class Nussinov : public Fold::Decoder
{
public:
  Nussinov(float th) : Fold::Decoder(), th_(th) { }
  float decode(float w, const VVF& p, const VVF& q, VU& ss);
  float score(float w, const VVF& p, const VVF& q, const VU& ss);
  float decode(const VVF& p, VU& ss, std::string& str);
  void make_brackets(const VU& ss, std::string& str) const;

private:
  float th_;
};

class SparseNussinov : public Fold::Decoder
{
public:
  SparseNussinov(float th) : Fold::Decoder(), th_(th) { }
  float decode(float w, const VVF& p, const VVF& q, VU& ss);
  float decode(float w, const VVF& p, const SparseFloatMatrix& q, VU& ss);
  float score(float w, const VVF& p, const VVF& q, const VU& ss);
  float decode(const VVF& p, VU& ss, std::string& str);
  void make_brackets(const VU& ss, std::string& str) const;

private:
  float th_;
};

class LinearNussinov : public Fold::Decoder
{
public:
  LinearNussinov(float th, uint beam_size) : th_(th), beam_size_(beam_size) { }
  float decode(float w, const VVF& p, const VVF& q, VU& ss);
  float decode(float w, const VVF& p, const SparseFloatMatrix& q, VU& ss);
  float decode(float w, const SparseFloatMatrix& p, const VVF& q, VU& ss);
  float decode(float w, const SparseFloatMatrix& p,
               const SparseFloatMatrix& q, VU& ss);
  LinearNussinovResult decode_certified(
      float w, const SparseFloatMatrix& p, const VVF& q, VU& ss);
  LinearNussinovResult decode_certified(
      float w, const SparseFloatMatrix& p,
      const SparseFloatMatrix& q, VU& ss);
  LinearNussinovResult decode_certified(
      float w, const SparseFloatMatrix& p, const VVF& q,
      const std::vector<std::vector<uint>>& pairs_by_right, VU& ss);
  LinearNussinovResult decode_certified(
      float w, const SparseFloatMatrix& p, const SparseFloatMatrix& q,
      const std::vector<std::vector<uint>>& pairs_by_right, VU& ss);
  float decode(const VVF& p, VU& ss, std::string& str);
  float decode(const SparseFloatMatrix& p, VU& ss, std::string& str);
  void make_brackets(const VU& ss, std::string& str) const;

private:
  template <typename PairScore>
  LinearNussinovResult decode_impl(
      uint length,
      const std::vector<std::vector<uint>>& pairs_by_right,
      PairScore pair_score, VU& ss, bool certify);

  float th_;
  uint beam_size_;
};

#endif  //  __INC_NUSSINOV_H__

// Local Variables:
// mode: C++
// End:
