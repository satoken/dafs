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

#ifndef __INC_NEEDLEMAN_WUNSCH_H__
#define __INC_NEEDLEMAN_WUNSCH_H__


#include "align.h"
#include <memory>

class BeamAlign;

class NeedlemanWunsch : public Align::Decoder
{
public:
  NeedlemanWunsch(float th) : Align::Decoder(), th_(th) { }
  void initialize(const VVF& p) { }
  float decode(const VVF& p, const VVF& q, VU& al) const;
  float decode(const VVF& p, VU& al) const;

private:
  float th_;
};

class SparseNeedlemanWunsch : public Align::Decoder
{
public:
  SparseNeedlemanWunsch(float th) : Align::Decoder(), th_(th), env_() { }
  void initialize(const VVF& p);
  float decode(const VVF& p, const VVF& q, VU& al) const;
  float decode(const VVF& p, const SparseFloatMatrix& q, VU& al) const;
  float decode(const VVF& p, VU& al) const;

private:
  float th_;
  std::vector< std::pair<uint,uint> > env_;
};

class LinearNeedlemanWunsch : public Align::Decoder
{
public:
  LinearNeedlemanWunsch(float th, uint beam_size);
  ~LinearNeedlemanWunsch();
  void initialize(const VVF&) { }
  void initialize(const SparseFloatMatrix&) { }
  float decode(const VVF& p, const VVF& q, VU& al) const;
  float decode(const VVF& p, const SparseFloatMatrix& q, VU& al) const;
  float decode(const VVF& p, VU& al) const;
  float decode(const SparseFloatMatrix& p, const VVF& q, VU& al) const;
  float decode(const SparseFloatMatrix& p,
               const SparseFloatMatrix& q, VU& al) const;
  float decode(const SparseFloatMatrix& p, VU& al) const;
  float similarity_score(const MP& p, uint length1, uint length2) const;

private:
  template <typename Probability, typename Multiplier>
  float decode_impl(uint length1, uint length2, Probability probability,
                    Multiplier multiplier, VU& al) const;

  float th_;
  std::unique_ptr<BeamAlign> aligner_;
};

#endif  //  __INC_NEEDLEMAN_WUNSCH_H__

// Local Variables:
// mode: C++
// End:
