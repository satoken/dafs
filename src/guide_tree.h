#ifndef DAFS_GUIDE_TREE_H
#define DAFS_GUIDE_TREE_H

#include <cstddef>
#include <string>
#include <utility>
#include <vector>
#include "fa.h"

namespace guide_tree {
using Node = std::pair<float, std::pair<unsigned int, unsigned int>>;

struct BeamStats {
  size_t band_cells = 0;
  size_t band_width = 0;
  bool certified_exact = false;
  size_t retained_states = 0;
  size_t peak_candidates = 0;
  size_t peak_retained = 0;
};

// MAFFT 7.526 nucleotide --globalpair distance and 0.1 mixed linkage.
// Leaves retain their input indices; parents follow their children.
// beam=0 uses exact affine DP. beam>0 first tries a certified exact band
// (half-width <=16*beam), then caps each state on each antidiagonal if needed.
double mafft_global_distance(const std::string& first, const std::string& second,
                             size_t beam = 0, BeamStats* stats = nullptr);
std::vector<Node> mafft_linkage(std::vector<std::vector<double>> distances);
std::vector<Node> mafft(const std::vector<Fasta>& sequences, size_t beam = 0);

// RAF 1.00 compatibility: input contains only the raw-posterior upper triangle.
// Its lower-triangle reads are intentional and match the benchmarked RAF tree.
std::vector<Node> raf(std::vector<std::vector<float>> distances);
} // namespace guide_tree

#endif
