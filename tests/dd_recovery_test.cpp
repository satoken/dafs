#include "dd_recovery.h"
#include <stdexcept>
#include <cmath>

void check(bool ok) { if (!ok) throw std::runtime_error("DD recovery test failed"); }
int main() {
  using namespace DDRecovery;
  check(monotone(VU{0, -1u, 2}, 3));
  check(!monotone(VU{0, 0}, 3));
  check(!monotone(VU{2, 1}, 3));
  check(!monotone(VU{3}, 3));
  check(structure(VU{5, 4, -1u, -1u, -1u, -1u}));
  check(structure(VU{5, 4, -1u, -1u, 1, 0}));
  check(!structure(VU{3, 4, -1u, -1u, -1u}));
  check(!structure(VU{4, 4, -1u, -1u, -1u}));
  check(!structure(VU{2, -1u, 3, -1u}));
  check(!structure(VU{-1u, 0}));
  check(coupled(VU{3,-1u,-1u,-1u}, VU{4,-1u,-1u,-1u,-1u}, VU{0,1,3,4}));
  check(!coupled(VU{3,-1u,-1u,-1u}, VU{3,-1u,-1u,-1u,-1u}, VU{0,1,3,4}));
  check(!coupled(VU{3,-1u,-1u,-1u}, VU{4,-1u,-1u,-1u,-1u}, VU{0,1,3,-1u}));
  AlignmentWindow window;
  const std::vector<std::pair<uint,uint>> support{{0,0},{1,2}};
  SparseFloatMatrix average;
  window.average(support, 2, 3, average);
  check(average.nonzeros() == 0);
  for (int t = 1; t <= 40; ++t) {
    window.push(support, [&](uint i, uint) { return float(t * (i+1)); });
    window.average(support, 2, 3, average);
    float expected = float(std::max(1,t-7)+t)/2;
    check(std::abs(average.get(0,0)-expected) < 1e-6);
    check(std::abs(average.get(1,2)-2*expected) < 1e-6);
    check(average.get(0,2) == 0);
  }
}
