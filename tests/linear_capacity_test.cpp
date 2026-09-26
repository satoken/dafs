#include <cmath>
#include <fstream>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "linearfold/fold/linfold.h"
#include "linearfold/param/contrafold.h"
#include "linearfold/param/turner.h"

template <typename Parameter>
bool check_capacity(const std::string& sequence, const char* model_name)
{
  using Fold = LinFold<Parameter>;
  auto model = std::make_unique<Parameter>(sequence);
  Fold fold(std::move(model));
  typename Fold::Options options;
  options.beam_size(100);
  options.probability_cutoff_ = 0.0f;
  options.set_allowed_pair('a', 'u');
  options.set_allowed_pair('u', 'a');
  options.set_allowed_pair('c', 'g');
  options.set_allowed_pair('g', 'c');
  options.set_allowed_pair('g', 'u');
  options.set_allowed_pair('u', 'g');

  fold.compute_inside(sequence, options);
  fold.compute_outside(sequence, options);
  const auto bpp = fold.compute_basepairing_probabilities(sequence, options);
  std::vector<double> incident(sequence.size() + 1, 0.0);
  for (size_t left = 1; left < bpp.size(); ++left)
    for (const auto& [right, probability] : bpp[left]) {
      if (right > sequence.size() || !std::isfinite(probability) ||
          probability < 0.0f || probability > 1.0f + 1e-6f) {
        std::cerr << model_name << " returned invalid BPP\n";
        return false;
      }
      incident[left] += probability;
      incident[right] += probability;
    }
  for (const double mass : incident)
    if (!std::isfinite(mass) || mass > 1.0 + 1e-6) {
      std::cerr << model_name << " exceeded nucleotide capacity: " << mass
                << "\n";
      return false;
    }
  return true;
}

int main(int argc, char** argv)
{
  if (argc != 2) {
    std::cerr << "usage: linear_capacity_test FASTA\n";
    return 2;
  }
  std::ifstream input(argv[1]);
  std::string line, source;
  while (std::getline(input, line)) {
    if (!line.empty() && line[0] != '>')
      source += line;
  }
  if (source.size() != 1400) {
    std::cerr << "unexpected fixture length: " << source.size() << "\n";
    return 2;
  }
  std::string sequence;
  sequence.reserve(2800);
  sequence = source + source;
  if (!check_capacity<CONTRAfoldNearestNeighbor>(sequence, "LPC") ||
      !check_capacity<TurnerNearestNeighbor>(sequence, "LPV"))
    return 1;
  std::cout << "linear partition BPP capacity check passed\n";
}
