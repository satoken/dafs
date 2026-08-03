#include "nussinov.h"

#include <cmath>
#include <iostream>
#include <random>

const uint Fold::Decoder::n_support_brackets = 4 + 26;
const char* Fold::Decoder::left_brackets = "([{<ABCDEFGHIJKLMNOPQRSTUVWXYZ";
const char* Fold::Decoder::right_brackets = ")]}>abcdefghijklmnopqrstuvwxyz";

namespace {
bool close(float lhs, float rhs)
{
  return std::fabs(lhs-rhs) <= 1e-5f;
}
}

int main()
{
  constexpr uint length = 10;
  constexpr float threshold = 0.2f;
  VVF probability(length, VF(length, 0.0f));
  VVF multiplier(length, VF(length, 0.0f));

  probability[0][9] = 0.90f;
  probability[1][8] = 0.80f;
  probability[3][6] = 0.65f;
  probability[0][5] = 0.75f;
  probability[4][9] = 0.70f;
  // The probability is absent, but a negative multiplier makes this pair
  // profitable.  Sparse support must therefore be the union of p and q.
  multiplier[2][7] = -0.80f;

  Nussinov exact(threshold);
  LinearNussinov linear(threshold, length);
  VU exact_structure, linear_structure;
  const float exact_score =
      exact.decode(1.0f, probability, multiplier, exact_structure);
  const float dense_linear_score =
      linear.decode(1.0f, probability, multiplier, linear_structure);
  if (!close(exact_score, dense_linear_score)) {
    std::cerr << "dense score mismatch: exact=" << exact_score
              << " linear=" << dense_linear_score << '\n';
    for (uint i = 0; i < length; ++i)
      if (exact_structure[i] != -1u || linear_structure[i] != -1u)
        std::cerr << i << ": exact=" << exact_structure[i]
                  << " linear=" << linear_structure[i] << '\n';
    return 1;
  }

  SparseFloatMatrix sparse_probability, sparse_multiplier;
  sparse_probability.assign(length, length);
  sparse_multiplier.assign(length, length);
  for (uint i = 0; i < length; ++i)
    for (uint j = i+1; j < length; ++j) {
      sparse_probability.set(i, j, probability[i][j]);
      sparse_multiplier.set(i, j, multiplier[i][j]);
    }

  const float sparse_linear_score =
      linear.decode(1.0f, sparse_probability, sparse_multiplier,
                    linear_structure);
  if (!close(exact_score, sparse_linear_score)) {
    std::cerr << "sparse score mismatch: exact=" << exact_score
              << " linear=" << sparse_linear_score << '\n';
    return 1;
  }
  if (linear_structure[2] != 7) {
    std::cerr << "q-only profitable pair was omitted\n";
    return 1;
  }

  const LinearNussinovResult full_certificate = linear.decode_certified(
      1.0f, sparse_probability, sparse_multiplier, linear_structure);
  if (!close(full_certificate.score, exact_score) ||
      !close(full_certificate.upper_bound, exact_score) ||
      full_certificate.pruned_states != 0) {
    std::cerr << "unpruned certificate mismatch: exact=" << exact_score
              << " beam=" << full_certificate.score
              << " upper=" << full_certificate.upper_bound
              << " pruned=" << full_certificate.pruned_states << '\n';
    return 1;
  }

  // Compare narrow-beam certificates with the exact cubic decoder over many
  // small deterministic random instances.  Every exact derivation either
  // reaches the root or has a first state removed from the beam; the latter
  // must be covered by the recorded completion bound.
  std::mt19937 generator(20260803u);
  std::uniform_real_distribution<float> score_distribution(0.05f, 1.0f);
  std::bernoulli_distribution edge_distribution(0.45);
  size_t strict_improvements = 0;
  size_t pruned_instances = 0;
  for (uint trial = 0; trial < 2000; ++trial) {
    constexpr uint random_length = 12;
    VVF dense_probability(random_length, VF(random_length, 0.0f));
    VVF dense_multiplier(random_length, VF(random_length, 0.0f));
    SparseFloatMatrix random_probability, random_multiplier;
    random_probability.assign(random_length, random_length);
    random_multiplier.assign(random_length, random_length);
    for (uint i = 0; i < random_length; ++i) {
      for (uint j = i+3; j < random_length; ++j) {
        const bool has_probability = edge_distribution(generator);
        const bool has_q_only_edge = !has_probability &&
            (i + 3*j + trial) % 13 == 0;
        if (!has_probability && !has_q_only_edge)
          continue;
        if (has_probability) {
          const float probability = score_distribution(generator);
          dense_probability[i][j] = probability;
          random_probability.set(i, j, probability);
        }
        if (has_q_only_edge || (i+j+trial) % 7 == 0) {
          const float multiplier = -0.25f * score_distribution(generator);
          dense_multiplier[i][j] = multiplier;
          random_multiplier.set(i, j, multiplier);
        }
      }
    }

    Nussinov random_exact(0.0f);
    VU random_exact_structure;
    const float random_exact_score = random_exact.decode(
        1.0f, dense_probability, dense_multiplier, random_exact_structure);
    VVU cached_support(random_length);
    for (uint i = 0; i < random_length; ++i)
      for (uint j = i+3; j < random_length; ++j)
        if (dense_probability[i][j] != 0.0f ||
            dense_multiplier[i][j] != 0.0f)
          cached_support[j].push_back(i);
    for (const uint beam : {1u, 2u, 3u, 5u}) {
      LinearNussinov random_linear(0.0f, beam);
      VU random_structure;
      const LinearNussinovResult result = random_linear.decode_certified(
          1.0f, random_probability, random_multiplier, random_structure);
      VU cached_structure;
      const LinearNussinovResult cached_result =
          random_linear.decode_certified(
              1.0f, random_probability, random_multiplier,
              cached_support, cached_structure);
      const float tolerance =
          2e-5f * std::max(1.0f, std::fabs(random_exact_score));
      if (result.score > random_exact_score + tolerance ||
          result.upper_bound + tolerance < random_exact_score ||
          result.upper_bound > result.additive_upper_bound + tolerance ||
          result.score > result.upper_bound + tolerance) {
        std::cerr << "invalid beam certificate: trial=" << trial
                  << " beam=" << beam
                  << " lower=" << result.score
                  << " exact=" << random_exact_score
                  << " upper=" << result.upper_bound
                  << " additive=" << result.additive_upper_bound << '\n';
        return 1;
      }
      if (!close(cached_result.score, result.score) ||
          !close(cached_result.upper_bound, result.upper_bound) ||
          cached_result.pruned_states != result.pruned_states ||
          cached_structure != random_structure) {
        std::cerr << "cached support changed certified decoding: trial="
                  << trial << " beam=" << beam << '\n';
        return 1;
      }
      strict_improvements +=
          result.upper_bound + tolerance < result.additive_upper_bound;
      pruned_instances += result.pruned_states != 0;
    }
  }
  if (strict_improvements == 0 || pruned_instances == 0) {
    std::cerr << "random certificate tests did not exercise pruning\n";
    return 1;
  }

  return 0;
}
