#ifndef DAFS_DD_CERTIFICATE_H
#define DAFS_DD_CERTIFICATE_H
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <vector>
#include "relaxed_bounds.h"

namespace DDCertificate {
inline double up(double x) { return std::nextafter(x, std::numeric_limits<double>::infinity()); }
inline double down(double x) { return std::nextafter(x, -std::numeric_limits<double>::infinity()); }
inline float lower_float(double x) {
  if (!std::isfinite(x)) throw std::runtime_error("non-finite lower bound");
  float f=static_cast<float>(x);
  if (static_cast<double>(f)>x)
    f=std::nextafter(f,-std::numeric_limits<float>::infinity());
  return f;
}
inline double lower_weighted_score(float weight,float probability,float threshold) {
  if (!(weight>=0) || !std::isfinite(weight) || !std::isfinite(probability) || !std::isfinite(threshold))
    throw std::runtime_error("invalid lower-bound coefficient");
  return down(static_cast<double>(weight)*down(static_cast<double>(probability)-threshold));
}

// Enforce monotonicity inside each fixed row block, but drop all constraints
// between blocks. Every full alignment is feasible for this relaxation.
// A 4-pass radix ordering on uint32 column indices avoids log(length).
// Duplicate edges are harmless maxima, and no column-sized dense buffer is
// allocated. Work O(width*M + rows + 1024), storage O(M + rows + width).
template<class Score>
float block_alignment(const std::vector<std::pair<unsigned,unsigned>>& support,
                      unsigned rows, unsigned columns, unsigned width, Score score) {
  if (width==0 || width>64) throw std::invalid_argument("certificate block width must be 1..64");
  std::vector<RelaxedBounds::WeightedEdge> edges;
  edges.reserve(support.size());
  for (const auto& [i,k]:support) {
    if (i>=rows || k>=columns) continue;
    const float v=score(i,k);
    if (!std::isfinite(v)) throw std::runtime_error("non-finite block score");
    if (v>0) edges.push_back({i,k,v});
  }
  std::vector<RelaxedBounds::WeightedEdge> scratch(edges.size());
  for (unsigned shift=0;shift<32;shift+=8) {
    std::array<size_t,256> count{},offset{};
    for (const auto& e:edges) ++count[(e.second>>shift)&255];
    for (size_t b=1;b<256;++b) offset[b]=offset[b-1]+count[b-1];
    for (const auto& e:edges) scratch[offset[(e.second>>shift)&255]++]=e;
    edges.swap(scratch);
  }
  const size_t blocks=(static_cast<size_t>(rows)+width-1)/width;
  std::vector<double> dp(blocks*(width+1),0.0), prefix(dp.size(),0.0);
  std::vector<uint64_t> seen(blocks,std::numeric_limits<uint64_t>::max());
  for (const auto& e:edges) {
    const size_t b=e.first/width, offset=b*(width+1);
    if (seen[b]!=e.second) {
      seen[b]=e.second;
      prefix[offset]=0.0;
      for (unsigned r=1;r<=width;++r)
        prefix[offset+r]=std::max(prefix[offset+r-1],dp[offset+r]);
    }
    const unsigned r=e.first%width;
    // All offers for this column use the frozen prefix, forbidding the
    // reuse of a column even when duplicate/tied entries are interleaved.
    dp[offset+r+1]=std::max(dp[offset+r+1],up(prefix[offset+r]+e.value));
  }
  double total=0.0;
  for (size_t b=0;b<blocks;++b) {
    double best=0.0;
    for (unsigned r=1;r<=width;++r) best=std::max(best,dp[b*(width+1)+r]);
    if (best>0) total=up(total+best);
  }
  return RelaxedBounds::round_up_to_float(total);
}

// Full max-alignment DP performs at most min(rows,columns) nonzero additions
// per path. Free skips mean an optimum contains no negative match. This
// guards double summation of the callback's FLOAT coefficients, not the
// earlier arithmetic that produced those coefficients.
inline float unpruned_alignment(double score, unsigned matches) {
  if (!std::isfinite(score) || score<0)
    throw std::runtime_error("invalid unpruned alignment score");
  const double ne=up(static_cast<double>(matches)*std::numeric_limits<double>::epsilon());
  if (ne>=0.5) throw std::runtime_error("unpruned alignment length overflow");
  const double denominator=std::nextafter(1.0-ne,0.0);
  return RelaxedBounds::round_up_to_float(up(score/denominator));
}
}
#endif
