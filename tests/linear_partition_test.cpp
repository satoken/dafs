#include <cmath>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include "linearfold/fold/linfold.h"
#include "linearfold/param/pair_objective.h"

using Fold = LinFold<PairObjectiveNearestNeighbor>;

static float pair_score(unsigned i, unsigned j, float lambda = 0.0f) {
  return 0.001f * i + 0.002f * j +
         lambda * (0.013f * i + 0.021f * j);
}

static std::vector<std::vector<int>> enumerate(
    const std::string& sequence, int lo, int hi, const Fold::Options& options) {
  const int n = static_cast<int>(sequence.size());
  if (lo == hi) return {std::vector<int>(n, 0)};

  std::vector<std::vector<int>> result;
  for (auto suffix : enumerate(sequence, lo + 1, hi, options))
    result.push_back(std::move(suffix));
  for (int right = lo + 1; right < hi; ++right) {
    if (!options.allow_paired(sequence, lo + 1, right + 1)) continue;
    for (const auto& inside : enumerate(sequence, lo + 1, right, options))
      for (const auto& suffix : enumerate(sequence, right + 1, hi, options)) {
        auto structure = inside;
        structure[lo] = right;
        structure[right] = lo;
        for (int k = right + 1; k < hi; ++k) structure[k] = suffix[k];
        result.push_back(std::move(structure));
      }
  }
  return result;
}

struct Oracle {
  long double z = 0;
  std::map<std::pair<int, int>, long double> bpp;
};

static Oracle oracle(const std::string& sequence, const Fold::Options& options,
                     float lambda = 0.0f) {
  Oracle answer;
  for (const auto& structure : enumerate(sequence, 0, sequence.size(), options)) {
    long double score = 0;
    for (int i = 0; i < static_cast<int>(sequence.size()); ++i)
      if (structure[i] > i)
        score += pair_score(i + 1, structure[i] + 1, lambda);
    const long double weight = std::exp(score);
    answer.z += weight;
    for (int i = 0; i < static_cast<int>(sequence.size()); ++i)
      if (structure[i] > i)
        answer.bpp[{i + 1, structure[i] + 1}] += weight;
  }
  for (auto& [pair, mass] : answer.bpp) mass /= answer.z;
  return answer;
}

static Fold::Options options(unsigned beam, unsigned max_helix, float lambda = 0.0f) {
  Fold::Options result;
  result.beam_size(beam).max_helix_length(max_helix)
      .min_hairpin_loop_length(3);
  result.set_allowed_pair('a', 'a');
  result.probability_cutoff_ = 0.0f;
  result.pair_score([lambda](unsigned i, unsigned j) {
    return pair_score(i, j, lambda);
  });
  return result;
}

static bool exact_small_cases() {
  for (const unsigned max_helix : {0u, 1u, 2u, 3u, 30u}) {
    auto model = std::make_unique<PairObjectiveNearestNeighbor>("");
    Fold fold(std::move(model));
    for (unsigned n = 0; n <= 12; ++n) {
      const std::string sequence(n, 'a');
      auto opt = options(0, max_helix);
      opt.constraints(std::vector<unsigned>(n + 1, Fold::Options::ANY));
      const auto expected = oracle(sequence, opt);
      const double log_z = fold.compute_inside(sequence, opt);
      fold.compute_outside(sequence, opt);
      const auto actual = fold.compute_basepairing_probabilities(sequence, opt);
      if (std::abs(std::exp(log_z) - expected.z) > 3e-3L) {
        std::cerr << "exact Z failure maxh=" << max_helix << " n=" << n << "\n";
        return false;
      }
      std::map<std::pair<int, int>, float> observed;
      for (unsigned i = 1; i < actual.size(); ++i)
        for (const auto& [j, probability] : actual[i]) observed[{i, j}] = probability;
      for (unsigned i = 1; i <= n; ++i)
        for (unsigned j = i + 1; j <= n; ++j)
          if (std::abs(observed[{i, j}] -
                       (expected.bpp.count({i, j}) ? expected.bpp.at({i, j}) : 0.0L)) > 5e-5L)
            { std::cerr << "exact BPP failure maxh=" << max_helix << " n=" << n << " pair=" << i << "," << j << "\n"; return false; }
    }
  }
  return true;
}

static bool constrained_cases() {
  for (unsigned n = 0; n <= 10; ++n) {
    for (unsigned mode = 0; mode < 4; ++mode) {
      const std::string sequence(n, 'a');
      auto opt = options(0, 30);
      std::vector<unsigned> constraint(n + 1, Fold::Options::ANY);
      if (mode & 1)
        for (unsigned i = 1; i <= n; ++i)
          if (i % 4 == 1) constraint[i] = Fold::Options::UNPAIRED;
      opt.constraints(constraint);
      if (mode & 2)
        opt.position_pair([](unsigned i, unsigned j) {
          return (3 * i + 5 * j) % 7 != 0;
        });
      const auto expected = oracle(sequence, opt);
      auto model = std::make_unique<PairObjectiveNearestNeighbor>(sequence);
      Fold fold(std::move(model));
      const double log_z = fold.compute_inside(sequence, opt);
      fold.compute_outside(sequence, opt);
      const auto actual = fold.compute_basepairing_probabilities(sequence, opt);
      if (std::abs(std::exp(log_z) - expected.z) > 3e-3L) {
        std::cerr << "constrained Z failure mode=" << mode << " n=" << n << "\n";
        return false;
      }
      std::map<std::pair<int, int>, float> observed;
      for (unsigned i = 1; i < actual.size(); ++i)
        for (const auto& [j, probability] : actual[i]) observed[{i, j}] = probability;
      for (unsigned i = 1; i <= n; ++i)
        for (unsigned j = i + 1; j <= n; ++j)
          if (std::abs(observed[{i, j}] -
                       (expected.bpp.count({i, j}) ? expected.bpp.at({i, j}) : 0.0L)) > 7e-5L)
            { std::cerr << "constrained BPP failure mode=" << mode << " n=" << n << " pair=" << i << "," << j << "\n"; return false; }
    }
  }
  return true;
}

static bool derivative_cases() {
  const std::string sequence(14, 'a');
  for (const unsigned beam : {0u, 1u, 2u, 3u, 5u, 100u}) {
    auto run = [&](float lambda) {
      auto opt = options(beam, 30, lambda);
      opt.constraints(std::vector<unsigned>(sequence.size() + 1,
                                           Fold::Options::ANY));
      auto model = std::make_unique<PairObjectiveNearestNeighbor>(sequence);
      Fold fold(std::move(model));
      const float log_z = fold.compute_inside(sequence, opt);
      fold.compute_outside(sequence, opt);
      const auto bpp = fold.compute_basepairing_probabilities(sequence, opt);
      double feature_expectation = 0;
      for (unsigned i = 1; i <= sequence.size(); ++i)
        for (const auto& [j, probability] : bpp[i])
          feature_expectation += probability * (0.013 * i + 0.021 * j);
      return std::pair<double, double>{log_z, feature_expectation};
    };
    const auto minus = run(-1e-3f);
    const auto plus = run(1e-3f);
    const double derivative = (plus.first - minus.first) / 2e-3;
    if (std::abs(derivative - run(0).second) > 1.5e-3) {
      const auto zero = run(0);
      std::cerr << "derivative failure beam=" << beam << " diff=" << (derivative-zero.second)
                << " d=" << derivative << " b=" << zero.second << " zm=" << minus.first
                << " zp=" << plus.first << "\n";
      return false;
    }
  }
  return true;
}

int main() {
  if (!exact_small_cases()) return 1;
  if (!constrained_cases()) return 1;
  if (!derivative_cases()) return 1;
  std::cout << "linear partition exhaustive inside/outside/BPP checks passed\n";
}
