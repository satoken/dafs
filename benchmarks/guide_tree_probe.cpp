#include "guide_tree.h"
#include <algorithm>
#include <chrono>
#include <functional>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

int main(int argc, char** argv)
{
  if (argc < 2 || argc > 3) {
    std::cerr << "Usage: guide-tree-probe INPUT.fa [BEAM=0]\n";
    return 1;
  }
  try {
    std::vector<Fasta> sequences;
    Fasta::load(sequences, argv[1]);
    size_t beam = 0;
    if (argc > 2) {
      const std::string argument = argv[2];
      if (argument.empty() || argument.find_first_not_of("0123456789") != std::string::npos)
        throw std::invalid_argument("BEAM must be a nonnegative integer");
      beam = std::stoul(argument);
    }
    const size_t n = sequences.size();
    std::vector<std::vector<double>> distances(n, std::vector<double>(n));
    guide_tree::BeamStats total;
    size_t certified_pairs = 0, band_cells = 0, band_width = 0;
    const auto start = std::chrono::steady_clock::now();
    for (size_t i = 0; i < n; ++i)
      for (size_t j = i + 1; j < n; ++j) {
        guide_tree::BeamStats stats;
        distances[i][j] = guide_tree::mafft_global_distance(
            sequences[i].seq(), sequences[j].seq(), beam, &stats);
        total.retained_states += stats.retained_states;
        total.peak_candidates = std::max(total.peak_candidates, stats.peak_candidates);
        total.peak_retained = std::max(total.peak_retained, stats.peak_retained);
        certified_pairs += stats.certified_exact;
        band_cells += stats.band_cells;
        band_width = std::max(band_width, stats.band_width);
      }
    const auto tree = guide_tree::mafft_linkage(distances);
    const double seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - start).count();
    unsigned int next = n;
    std::function<unsigned int(unsigned int)> postorder = [&](unsigned int node) {
      if (node < n) return node;
      const auto left = postorder(tree[node].second.first);
      const auto right = postorder(tree[node].second.second);
      const auto parent = next++;
      std::cout << "m " << parent << ' ' << left << ' ' << right << '\n';
      return parent;
    };
    postorder(tree.size() - 1);
    std::cout << std::setprecision(17);
    for (size_t i = 0; i < n; ++i)
      for (size_t j = i + 1; j < n; ++j)
        std::cout << "d " << i << ' ' << j << ' ' << distances[i][j] << '\n';
    std::cout << "stats " << seconds << ' ' << total.retained_states << ' '
              << total.peak_candidates << ' ' << total.peak_retained << '\n';
    std::cout << "fast " << certified_pairs << ' ' << band_cells << ' ' << band_width << '\n';
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
