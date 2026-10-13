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
#include <memory>

#include "fold.h"

struct LinearNussinovResult
{
  float score = 0.0f;
  float upper_bound = 0.0f;
  float pruned_upper_bound = 0.0f;
  float additive_upper_bound = 0.0f;
  size_t pruned_states = 0;
};

// Independent dense reference implementation; the CLI uses the faster
// IPknot-backed SparseNussinov for exact decoding.
class Nussinov : public Fold::Decoder
{
public:
  Nussinov(float th, bool allow_isolated = true, int dp_threads = 1);
  float decode(float w, const VVF& p, const VVF& q, VU& ss);
  float score(float w, const VVF& p, const VVF& q, const VU& ss);
  float decode(const VVF& p, VU& ss, std::string& str);
  float decode(const SparseFloatMatrix& p,
               const SparseFloatMatrix& pair_bonus,
               VU& ss, std::string& str);
  void make_brackets(const VU& ss, std::string& str) const;
  int used_workers() const { return used_workers_; }
  size_t parallel_steps() const { return parallel_steps_; }

private:
  float th_;
  bool allow_isolated_;
  int dp_threads_, used_workers_ = 1;
  size_t parallel_steps_ = 0;
};

class SparseNussinov : public Fold::Decoder
{
public:
  SparseNussinov(float th, bool allow_isolated = true, int dp_threads = 1);
  ~SparseNussinov() override;
  SparseNussinov(const SparseNussinov&) = delete;
  SparseNussinov& operator=(const SparseNussinov&) = delete;
  float decode(float w, const VVF& p, const VVF& q, VU& ss);
  float decode(float w, const VVF& p, const SparseFloatMatrix& q, VU& ss);
  float decode(float w, const SparseFloatMatrix& p, const VVF& q, VU& ss) override;
  float decode(float w, const SparseFloatMatrix& p,
               const SparseFloatMatrix& q, VU& ss) override;
  float score(float w, const VVF& p, const VVF& q, const VU& ss);
  float decode(const VVF& p, VU& ss, std::string& str);
  float decode(const SparseFloatMatrix& p, VU& ss, std::string& str) override;
  float decode(const SparseFloatMatrix& p,
               const SparseFloatMatrix& pair_bonus,
               VU& ss, std::string& str);
  void make_brackets(const VU& ss, std::string& str) const;

  // Exact closed components omitted by IPknot's decomposition criterion.
  size_t sparsified_candidates() const;
  int used_workers() const;
  size_t parallel_steps() const;

private:
  struct Workspace;
  template <typename PairScore>
  float decode_impl(uint length, const std::vector<std::vector<uint>>* support,
                    PairScore pair_score, VU& ss);
  float th_;
  bool allow_isolated_;
  int dp_threads_;
  std::unique_ptr<Workspace> workspace_;
};

// Chooses the measured faster beam kernel, preserving DAFS certificates.
class LinearNussinov : public Fold::Decoder
{
public:
  static constexpr uint cached_support_minimum_pair_span = 2;

  LinearNussinov(float th, uint beam_size, bool allow_isolated = true,
                 int dp_threads = 1);
  ~LinearNussinov() override;
  LinearNussinov(const LinearNussinov&) = delete;
  LinearNussinov& operator=(const LinearNussinov&) = delete;
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
  float decode(const SparseFloatMatrix& p,
               const SparseFloatMatrix& pair_bonus,
               VU& ss, std::string& str);
  void make_brackets(const VU& ss, std::string& str) const;
  int used_workers() const;
  size_t parallel_steps() const;
  bool used_ipknot_kernel() const;

private:
  struct Workspace;
  struct DafsWorkspace;

  template <typename PairScore>
  LinearNussinovResult decode_impl(
      uint length,
      const std::vector<std::vector<uint>>& pairs_by_right,
      PairScore pair_score, VU& ss, bool certify);

  template <typename PairScore>
  LinearNussinovResult decode_dafs_impl(
      uint length, const std::vector<std::vector<uint>>& pairs_by_right,
      PairScore pair_score, VU& ss, bool certify);

  float th_;
  uint beam_size_;
  bool allow_isolated_;
  int dp_threads_;
  std::unique_ptr<Workspace> workspace_;
};

#endif  //  __INC_NUSSINOV_H__

// Local Variables:
// mode: C++
// End:
