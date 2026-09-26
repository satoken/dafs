#ifndef DAFS_DD_RECOVERY_H
#define DAFS_DD_RECOVERY_H

#include <array>
#include <algorithm>
#include "sparse_matrix.h"

namespace DDRecovery {
// Only a primal proposal: omitted q-only edges cannot invalidate an upper
// bound, since this average is NEVER used for certification or gradients.
// Fixed support and eight slots give O(8 * support) storage, independent of T.
class AlignmentWindow {
public:
  static constexpr size_t width = 8;
  template<class Lookup>
  void push(const std::vector<std::pair<uint, uint>>& support, Lookup q) {
    auto& slot = slots_[next_];
    slot.resize(support.size());
    for (size_t i = 0; i < support.size(); ++i)
      slot[i] = q(support[i].first, support[i].second);
    next_ = (next_ + 1) % width;
    count_ = std::min(count_ + 1, width);
  }
  void average(const std::vector<std::pair<uint, uint>>& support,
               uint rows, uint columns, SparseFloatMatrix& out) const {
    out.assign(rows, columns);
    if (!count_) return;
    for (size_t i = 0; i < support.size(); ++i) {
      double sum = 0.0;
      for (size_t s = 0; s < count_; ++s) sum += slots_[s][i];
      out.set(support[i].first, support[i].second,
              static_cast<float>(sum / count_));
    }
  }
private:
  std::array<VF, width> slots_;
  size_t next_ = 0, count_ = 0;
};

inline bool monotone(const VU& z, uint columns) {
  uint previous = 0;
  bool have_previous = false;
  for (uint k : z) {
    if (k == -1u) continue;
    if (k >= columns || (have_previous && k <= previous)) return false;
    previous = k;
    have_previous = true;
  }
  return true;
}

// VU may contain either forward-only or symmetric partners. Check endpoint
// uniqueness, noncrossing, and the forward pairs' common monotone mapping.
inline bool structure(const VU& x) {
  VU ends;
  std::vector<bool> used(x.size(), false);
  for (uint i = 0; i < x.size(); ++i) {
    while (!ends.empty() && ends.back() < i) ends.pop_back();
    const uint j = x[i];
    if (j == -1u) continue;
    if (j >= x.size() || j == i) return false;
    if (j < i) { if (x[j] != i) return false; continue; }
    if (used[i] || used[j] || (!ends.empty() && j >= ends.back())) return false;
    used[i] = used[j] = true;
    ends.push_back(j);
  }
  return true;
}

inline bool coupled(const VU& x, const VU& y, const VU& z) {
  if (x.size() != z.size() || !monotone(z, y.size()) ||
      !structure(x) || !structure(y)) return false;
  size_t nx = 0, ny = 0;
  for (uint k = 0; k < y.size(); ++k)
    if (y[k] != -1u && y[k] > k) ++ny;
  for (uint i = 0; i < x.size(); ++i) {
    const uint j = x[i];
    if (j == -1u || j < i) continue;
    ++nx;
    if (z[i] == -1u || z[j] == -1u || y[z[i]] != z[j]) return false;
  }
  return nx == ny;
}
} // namespace DDRecovery
#endif
