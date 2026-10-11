#include "guide_tree.h"
#include <cmath>
#include <iostream>
#include <limits>
#include <random>
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
             {"ARYKMBDHVNX", "AGCUAGCUAGC"},
             {"AAAAAGGGGGCCCCCUUUUU", "GGGGGCCCCC"}}) {
      const double forward = mafft_global_distance(pair.first, pair.second);
      const double reverse = mafft_global_distance(pair.second, pair.first);
      require(forward >= 0 && std::isfinite(forward), "Non-finite distance");
      require(std::abs(forward - reverse) < 1e-12, "Asymmetric MAFFT distance");
      require(forward == mafft_global_distance(pair.first, pair.second, 100),
              "Unpruned beam changed MAFFT gap rules");
    }
    std::mt19937 random(4291);
    const std::string bases = "ACGUNRYX";
    for (size_t trial = 0; trial < 300; ++trial) {
      std::string a(1 + random() % 50, 'A'), b(1 + random() % 50, 'A');
      for (auto& base : a) base = bases[random() % bases.size()];
      for (auto& base : b) base = bases[random() % bases.size()];
      const double exact = mafft_global_distance(a, b);
      require(exact == mafft_global_distance(a, b, 51), "Wide beam differs from exact DP");
      for (const size_t beam : {size_t{1}, size_t{3}, size_t{10}}) {
        BeamStats stats;
        const double approximate = mafft_global_distance(a, b, beam, &stats);
        require(std::isfinite(approximate) && approximate >= exact - 1e-12,
                "Beam score exceeds exact optimum or has no path");
        require(approximate == mafft_global_distance(b, a, beam), "Asymmetric beam distance");
        require(stats.peak_retained <= beam, "Beam cap exceeded on score ties");
        require(stats.peak_candidates <= 3 * beam, "Candidate bound exceeded");
        require(stats.retained_states <= 3 * beam * (a.size() + b.size()),
                "Visited-state bound exceeded");
        require(stats.band_cells <= 66 * beam * (a.size() + b.size()),
                "Certified band exceeded linear work bound");
        if (stats.certified_exact)
          require(approximate == exact, "Incorrect exact-score certificate");
      }
    }
    BeamStats long_stats;
    const std::string long_sequence(5000, 'A');
    require(mafft_global_distance(long_sequence, long_sequence, 3, &long_stats) == 0,
            "Long identical beam distance");
    require(long_stats.retained_states <= 6 * 3 * long_sequence.size(),
            "Long input did not respect linear state bound");
    require(long_stats.certified_exact && long_stats.band_cells > 0 &&
            long_stats.band_cells <= 66 * 3 * 2 * long_sequence.size(),
            "Long identical input did not use bounded certified band");
    BeamStats rejected_certificate;
    mafft_global_distance(std::string(100, 'A'), std::string(100, 'G'), 1, &rejected_certificate);
    require(!rejected_certificate.certified_exact && rejected_certificate.band_cells > 0 &&
            rejected_certificate.retained_states > 0,
            "Unproved band result did not fall back to beam search");
    BeamStats unequal_stats;
    mafft_global_distance("ACGU", std::string(100, 'A') + "ACGU", 3, &unequal_stats);
    require(!unequal_stats.certified_exact && unequal_stats.band_cells == 0 &&
            unequal_stats.retained_states > 0, "Large length difference did not skip band");
    size_t certificates = 0;
    for (size_t trial = 0; trial < 80; ++trial) {
      std::string a(64 + random() % 193, 'A');
      for (auto& base : a) base = "ACGU"[random() % 4];
      std::string b = a;
      for (auto& base : b) if (random() % 10 == 0) base = "ACGU"[random() % 4];
      b.insert(random() % (b.size() + 1), random() % 21, 'G');
      BeamStats stats;
      const double approximate = mafft_global_distance(a, b, 16, &stats);
      const double exact = mafft_global_distance(a, b);
      if (stats.certified_exact) {
        ++certificates;
        require(approximate == exact, "Indel certificate differs from exact DP");
      }
      require(approximate >= exact - 1e-12, "Certified/beam score exceeds exact optimum");
      require(stats.band_cells <= 66 * 16 * (a.size() + b.size()), "Indel band work bound");
    }
    require(certificates > 0, "Indel fixtures never exercised certification");
    for (const size_t length : {size_t{1}, size_t{2}, size_t{50}, size_t{500}})
      require(std::isfinite(mafft_global_distance(std::string(length, 'C'), "AGU", 1)),
              "Narrow beam failed on unequal lengths");
    const std::vector<std::vector<double>> distances = {
        {0, .1, .6, .7}, {0, 0, .5, .65}, {0, 0, 0, .2}, {0, 0, 0, 0}};
    topology(mafft_linkage(distances), {{0, 1}, {2, 3}, {4, 5}});
    topology(mafft_linkage(std::vector<std::vector<double>>(4, std::vector<double>(4))),
             {{0, 1}, {4, 2}, {5, 3}});
    topology(mafft({Fasta("same name", "A"), Fasta("same name", "A"),
                    Fasta("same name", "A"), Fasta("same name", "A")}),
             {{0, 1}, {4, 2}, {5, 3}});
    topology(mafft({Fasta("one", "A")}), {});
    topology(mafft({Fasta("same", "A"), Fasta("same", "A"),
                    Fasta("same", "A"), Fasta("same", "A")}, 1),
             {{0, 1}, {4, 2}, {5, 3}});
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
