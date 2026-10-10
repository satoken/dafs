#include "ipknot.h"

#include <cassert>
#include <cmath>
#include <string>
#include <utility>

// IPknot only needs these Fold symbols to render bracket strings.
const uint Fold::Decoder::n_support_brackets = 30;
const char* Fold::Decoder::left_brackets = "([{<ABCDEFGHIJKLMNOPQRSTUVWXYZ";
const char* Fold::Decoder::right_brackets = ")]}>abcdefghijklmnopqrstuvwxyz";

int main()
{
  if (!ipknot::Decoder::ilp_available())
    return 77;
  constexpr uint length = 10;
  VVF dense_probability(length, VF(length, 0.0f));
  VVF dense_multiplier(length, VF(length, 0.0f));
  SparseFloatMatrix sparse_probability, sparse_multiplier;
  sparse_probability.assign(length, length);
  sparse_multiplier.assign(length, length);

  for (const auto [i, j] : {std::pair<uint, uint>{0, 9},
                            {1, 8}, {2, 7}}) {
    dense_probability[i][j] = 0.9f;
    sparse_probability.set(i, j, 0.9f);
  }
  // DD can make a pair profitable even if LinearPartition did not emit it.
  dense_multiplier[3][6] = -1.5f;
  sparse_multiplier.set(3, 6, -1.5f);

  IPknot dense_decoder(VF{0.2f});
  IPknot sparse_decoder(VF{0.2f});
  VU dense_result, sparse_result;
  const float dense_score = dense_decoder.decode(
      1.0f, dense_probability, dense_multiplier, dense_result);
  const float sparse_score = sparse_decoder.decode(
      1.0f, sparse_probability, sparse_multiplier, sparse_result);
  assert(dense_result == sparse_result);
  assert(sparse_result[3] == 6);
  assert(std::fabs(dense_score - sparse_score) < 1e-5f);
  assert(std::fabs(sparse_score - 3.4f) < 1e-5f);
  const float mixed_score = sparse_decoder.decode(
      1.0f, sparse_probability, dense_multiplier, sparse_result);
  assert(dense_result == sparse_result);
  assert(std::fabs(dense_score - mixed_score) < 1e-5f);
  const float mixed_score2 = sparse_decoder.decode(
      1.0f, dense_probability, sparse_multiplier, sparse_result);
  assert(dense_result == sparse_result);
  assert(std::fabs(dense_score - mixed_score2) < 1e-5f);

  IPknot dense_final(VF{0.2f});
  IPknot sparse_final(VF{0.2f});
  std::string dense_brackets, sparse_brackets;
  dense_final.decode(dense_probability, dense_result, dense_brackets);
  sparse_final.decode(sparse_probability, sparse_result, sparse_brackets);
  assert(dense_result == sparse_result);
  assert(dense_brackets == sparse_brackets);
  assert(dense_brackets == "(((....)))");

  // Cross two stacked stems. Check upstream ILP levels and unaveraged scores.
  ipknot::DecoderOptions options;
  options.backend = ipknot::Backend::ILP;
  options.thresholds = {0.2f, 0.1f};
  options.weights = {1.0f, 1.0f};
  const std::vector<ipknot::PairProbability> pairs{
      {0, 7, 0.9f}, {1, 6, 0.9f}, {3, 10, 0.9f}, {4, 9, 0.9f}};
  const auto reference = ipknot::Decoder(options).decode(15, pairs);
  SparseFloatMatrix crossing;
  crossing.assign(15, 15);
  for (const auto& pair : pairs)
    crossing.set(pair.left, pair.right, pair.probability);
  IPknot layered(options.thresholds);
  const float layered_score = layered.decode(crossing, sparse_result, sparse_brackets);
  assert(std::fabs(layered_score - reference.objective) < 1e-5);
  assert(std::fabs(layered_score - 3.0f) < 1e-5);
  for (uint i = 0; i < sparse_result.size(); ++i) {
    assert(sparse_result[i] == (reference.bpseq[i] > static_cast<int>(i)
                                   ? static_cast<uint>(reference.bpseq[i]) : -1u));
    if (sparse_result[i] != -1u)
      assert(sparse_brackets[i] == Fold::Decoder::left_brackets[reference.levels[i]]);
  }

  // Final profile/refinement scores can exceed a probability of one.
  sparse_probability.scale(2.0f);
  const float profile_score = sparse_final.decode(
      sparse_probability, sparse_result, sparse_brackets);
  assert(std::fabs(profile_score - 4.8f) < 1e-5);

  SparseFloatMatrix lonely, no_multiplier;
  lonely.assign(10, 10);
  no_multiplier.assign(10, 10);
  lonely.set(0, 8, 0.95f);
  assert(sparse_final.decode(lonely, sparse_result, sparse_brackets) == 0.0f);
  IPknot allowing_decoder(VF{0.2f}, 1, true);
  assert(std::fabs(allowing_decoder.decode(lonely, sparse_result, sparse_brackets)
                   - 0.75f) < 1e-5);
  assert(sparse_result[0] == 8 && sparse_brackets == "(.......).");
  assert(std::fabs(allowing_decoder.decode(2.0f, lonely, no_multiplier, sparse_result)
                   - 1.5f) < 1e-5);
  assert(sparse_result[0] == 8);
}
