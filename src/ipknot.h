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
#include <ipknot/decoder.h>

// Adapt DAFS's weighted posterior/multiplier objective and one-sided structure
// representation to the reusable IPknot decoder's public API.
class IPknot : public Fold::Decoder
{
public:
  IPknot(const VF& th, int n_th = 1, bool allow_isolated = true);
  float decode(float w, const VVF& p, const VVF& q, VU& ss) override;
  float decode(float w, const VVF& p, const SparseFloatMatrix& q,
               VU& ss) override;
  float decode(float w, const SparseFloatMatrix& p, const VVF& q,
               VU& ss) override;
  float decode(float w, const SparseFloatMatrix& p,
               const SparseFloatMatrix& q, VU& ss) override;
  float decode(const VVF& p, VU& ss, std::string& str) override;
  float decode(const SparseFloatMatrix& p, VU& ss,
               std::string& str) override;
  void make_brackets(const VU& ss, std::string& str) const override;

protected:
  explicit IPknot(ipknot::DecoderOptions options);

private:
  struct ScoredPair { uint left, right; double score; };
  float decode_scores(float w, uint length, const std::vector<ScoredPair>& pairs,
                      VU& ss);
  static VU decompose_levels(const VU& ss);

  ipknot::DecoderOptions options_;
  VU levels_, last_structure_;
};

// Bounded beam DD from the same IPknot library, independent of DAFS's outer
// coupling-iteration budget and usable without an integer-programming solver.
class LinearIPknot : public IPknot
{
public:
  LinearIPknot(const VF& th, uint beam_size, bool allow_isolated = true,
               int dp_threads = 1);
};

#endif // __INC_IPKNOT_H__
