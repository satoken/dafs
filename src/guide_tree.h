#ifndef DAFS_GUIDE_TREE_H
#define DAFS_GUIDE_TREE_H

#include <string>
#include <utility>
#include <vector>
#include "fa.h"

namespace guide_tree {
using Node = std::pair<float, std::pair<unsigned int, unsigned int>>;

// MAFFT 7.526 nucleotide --globalpair distance and 0.1 mixed linkage.
// Leaves retain their input indices; parents follow their children.
double mafft_global_distance(const std::string& first, const std::string& second);
std::vector<Node> mafft_linkage(std::vector<std::vector<double>> distances);
std::vector<Node> mafft(const std::vector<Fasta>& sequences);

// RAF 1.00 compatibility: input contains only the raw-posterior upper triangle.
// Its lower-triangle reads are intentional and match the benchmarked RAF tree.
std::vector<Node> raf(std::vector<std::vector<float>> distances);
} // namespace guide_tree

#endif
