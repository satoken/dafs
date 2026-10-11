#include "ipknot.h"

#include <cassert>
#include <cmath>
#include <string>
#include <vector>

const uint Fold::Decoder::n_support_brackets = 30;
const char* Fold::Decoder::left_brackets = "([{<ABCDEFGHIJKLMNOPQRSTUVWXYZ";
const char* Fold::Decoder::right_brackets = ")]}>abcdefghijklmnopqrstuvwxyz";

void check_structure_brackets(const VU& structure,
                              const std::string& brackets)
{
  assert(structure.size() == brackets.size());
  std::vector<char> right_endpoint(structure.size(), false);
  const std::string left_brackets(Fold::Decoder::left_brackets);
  const std::string right_brackets(Fold::Decoder::right_brackets);
  for (uint i = 0; i < structure.size(); ++i) {
    const uint j = structure[i];
    if (j == -1u)
      continue;
    assert(i < j && j < structure.size());
    assert(!right_endpoint[j]);
    right_endpoint[j] = true;
    const std::string::size_type level = left_brackets.find(brackets[i]);
    assert(level != std::string::npos);
    assert(brackets[j] == right_brackets[level]);
  }
  for (uint i = 0; i < structure.size(); ++i) {
    if (structure[i] == -1u && !right_endpoint[i])
      assert(brackets[i] == '.');
    if (structure[i] != -1u)
      assert(!right_endpoint[i]);
  }
}

int main()
{
  constexpr uint length = 15; // Include trailing unpaired positions.
  SparseFloatMatrix posterior;
  posterior.assign(length, length);
  const std::vector<ipknot::PairProbability> pairs{
      {0, 7, 0.9f}, {1, 6, 0.9f}, {3, 10, 0.9f}, {4, 9, 0.9f}};
  for (const auto& pair : pairs)
    posterior.set(pair.left, pair.right, pair.probability);

  // Compare against the upstream decoder, including its level assignments.
  ipknot::DecoderOptions options;
  options.thresholds = {0.2f, 0.1f};
  options.weights = {1.0f, 1.0f};
  options.dd.beam = 32;
  const auto reference = ipknot::Decoder(options).decode(length, pairs);
  LinearIPknot decoder(options.thresholds, 32);
  VU structure;
  std::string brackets;
  const float score = decoder.decode(posterior, structure, brackets);
  assert(structure.size() == length && brackets.size() == length);
  assert(structure[0] == 7 && structure[1] == 6);
  assert(structure[3] == 10 && structure[4] == 9);
  for (uint i = 0; i < length; ++i) {
    assert(structure[i] == (reference.bpseq[i] > static_cast<int>(i)
                               ? static_cast<uint>(reference.bpseq[i]) : -1u));
    if (structure[i] != -1u)
      assert(brackets[i] == Fold::Decoder::left_brackets[reference.levels[i]]);
  }
  check_structure_brackets(structure, brackets);
  assert(std::fabs(score - reference.objective) < 1e-5);
  assert(std::fabs(score - 3.0f) < 1e-5);

  // Projected structures have no cached levels and still need valid brackets.
  VU projected(length, -1u);
  projected[0] = 5;
  projected[2] = 8;
  decoder.make_brackets(projected, brackets);
  check_structure_brackets(projected, brackets);

  // Signed multipliers can introduce stacked pairs outside the BPP support,
  // with effective scores above one. Test all dense/sparse combinations.
  SparseFloatMatrix p, q;
  p.assign(12, 12);
  q.assign(12, 12);
  p.set(0, 11, 0.9f);
  p.set(1, 10, 0.85f);
  q.set(0, 11, 0.15f);
  q.set(2, 9, -1.5f);
  q.set(3, 8, -1.25f);
  LinearIPknot weighted_decoder(VF{0.2f}, 32);
  VU sparse_result, dense_result, mixed_result;
  const float sparse_score = weighted_decoder.decode(2.0f, p, q, sparse_result);
  const float dense_score = weighted_decoder.decode(2.0f, p.dense(), q.dense(), dense_result);
  assert(sparse_result == dense_result);
  assert(sparse_result[2] == 9 && sparse_result[3] == 8);
  assert(std::fabs(sparse_score - 4.5f) < 1e-5);
  assert(std::fabs(sparse_score - dense_score) < 1e-5);
  assert(std::fabs(weighted_decoder.decode(2.0f, p, q.dense(), mixed_result)
                   - sparse_score) < 1e-5);
  assert(mixed_result == sparse_result);
  assert(std::fabs(weighted_decoder.decode(2.0f, p.dense(), q, mixed_result)
                   - sparse_score) < 1e-5);
  assert(mixed_result == sparse_result);

  // Zero structure weight still leaves the multiplier objective active.
  assert(std::fabs(weighted_decoder.decode(0.0f, p, q, structure) - 2.75f) < 1e-5);
  assert(structure[2] == 9 && structure[3] == 8);

  // The new DD decoder enforces IPknot's no-lonely-pairs rule.
  SparseFloatMatrix lonely;
  lonely.assign(length, length);
  lonely.set(0, 8, 0.95f);
  assert(decoder.decode(lonely, structure, brackets) == 0.0f);
  assert(brackets == std::string(length, '.'));
  LinearIPknot allowing_decoder(options.thresholds, 32, true);
  assert(std::fabs(allowing_decoder.decode(lonely, structure, brackets) - 0.75f) < 1e-5);
  assert(structure[0] == 8);
  check_structure_brackets(structure, brackets);
  SparseFloatMatrix no_multiplier;
  no_multiplier.assign(length, length);
  assert(std::fabs(allowing_decoder.decode(2.0f, lonely, no_multiplier, structure)
                   - 1.5f) < 1e-5);
  assert(structure[0] == 8);

  // Allow two isolated crossing pairs on separate levels.
  lonely.set(3, 11, 0.9f);
  assert(decoder.decode(lonely, structure, brackets) == 0.0f);
  assert(std::fabs(allowing_decoder.decode(lonely, structure, brackets) - 1.55f) < 1e-5);
  assert(structure[0] == 8 && structure[3] == 11);
  assert(brackets[0] != brackets[3]);
  check_structure_brackets(structure, brackets);

  // Empty inputs and decoder reuse must not retain old structures or levels.
  SparseFloatMatrix empty;
  empty.assign(0, 0);
  assert(decoder.decode(empty, structure, brackets) == 0.0f);
  assert(structure.empty() && brackets.empty());
  if (!ipknot::Decoder::ilp_available()) {
    bool rejected = false;
    try {
      IPknot exact(VF{0.2f});
      exact.decode(empty, structure, brackets);
    } catch (const std::runtime_error&) {
      rejected = true;
    }
    assert(rejected);
  }
}
