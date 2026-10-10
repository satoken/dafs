#include "guide_tree.h"
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>

namespace {
void require(bool condition, const char* message)
{
  if (!condition) throw std::runtime_error(message);
}
void topology(const std::vector<guide_tree::Node>& tree,
              const std::vector<std::pair<unsigned int, unsigned int>>& expected)
{
  const size_t n = expected.size() + 1;
  require(tree.size() == 2 * n - 1, "Wrong guide-tree size");
  for (size_t i = 0; i < expected.size(); ++i)
    require(tree[n + i].second == expected[i], "Wrong ordered guide-tree merge");
}
}

int main()
{
  try {
    using namespace guide_tree;
    require(mafft_global_distance("AGCU", "AGCU") == 0, "Identical distance");
    require(mafft_global_distance("AGCU", "agct") == 0, "U/T and case handling");
    require(mafft_global_distance("NNNX", "AGCU") == 2, "Zero self-score distance");
    require(mafft_global_distance("", "AGCU") == 2, "Empty self-score distance");
    for (const auto& pair : std::vector<std::pair<std::string, std::string>>{
             {"AGCUAGCU", "AGCU"}, {"AGCUAGCU", "UUUUGGGG"},
             {"ARYKMBDHVNX", "AGCUAGCUAGC"}}) {
      const double forward = mafft_global_distance(pair.first, pair.second);
      const double reverse = mafft_global_distance(pair.second, pair.first);
      require(forward >= 0 && std::isfinite(forward), "Non-finite distance");
      require(std::abs(forward - reverse) < 1e-12, "Asymmetric MAFFT distance");
    }
    const std::vector<std::vector<double>> distances = {
        {0, .1, .6, .7}, {0, 0, .5, .65}, {0, 0, 0, .2}, {0, 0, 0, 0}};
    topology(mafft_linkage(distances), {{0, 1}, {2, 3}, {4, 5}});
    topology(mafft_linkage(std::vector<std::vector<double>>(4, std::vector<double>(4))),
             {{0, 1}, {4, 2}, {5, 3}});
    topology(mafft({Fasta("same name", "A"), Fasta("same name", "A"),
                    Fasta("same name", "A"), Fasta("same name", "A")}),
             {{0, 1}, {4, 2}, {5, 3}});
    topology(mafft({Fasta("one", "A")}), {});
    topology(raf({{0, .1f, .6f, .7f}, {0, 0, .5f, .65f},
                  {0, 0, 0, .2f}, {0, 0, 0, 0}}), {{0, 1}, {4, 2}, {5, 3}});
    topology(raf({{0, .2f, .6f, .5f}, {0, 0, .7f, .65f},
                  {0, 0, 0, .1f}, {0, 0, 0, 0}}), {{2, 3}, {0, 1}, {5, 4}});
    topology(raf({{0}}), {});
    bool rejected = false;
    try { mafft({}); } catch (const std::invalid_argument&) { rejected = true; }
    require(rejected, "Empty tree accepted");
    rejected = false;
    try { mafft_global_distance("Z", "A"); }
    catch (const std::invalid_argument&) { rejected = true; }
    require(rejected, "Unsupported nucleotide accepted");
    std::cout << "Guide-tree algorithm tests passed\n";
    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
