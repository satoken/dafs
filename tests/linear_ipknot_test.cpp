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
  constexpr uint length = 12;
  SparseFloatMatrix posterior;
  posterior.assign(length, length);
  // The first layer is non-crossing.  (1,8) crosses both selected arcs and
  // is therefore available to the second IPknot layer.
  posterior.set(0, 5, 0.95f);
  posterior.set(6, 11, 0.90f);
  posterior.set(1, 8, 0.92f);

  LinearIPknot decoder(VF{0.2f, 0.1f}, 32);
  VU structure;
  std::string brackets;
  const float score = decoder.decode(posterior, structure, brackets);
  assert(structure.size() == length);
  assert(structure[0] == 5 && structure[5] == -1u);
  assert(structure[6] == 11 && structure[11] == -1u);
  assert(structure[1] == 8 && structure[8] == -1u);
  check_structure_brackets(structure, brackets);
  assert(std::isfinite(score) && score > 0.0f);

  // make_brackets() also receives projected structures during verbose
  // decoding.  This structure has a crossing pair and therefore exercises
  // the linear fallback coloring instead of the cached decoder levels.
  VU projected(length, -1u);
  projected[0] = 5;
  projected[2] = 8;
  std::string projected_brackets;
  decoder.make_brackets(projected, projected_brackets);
  check_structure_brackets(projected, projected_brackets);

  // A larger right endpoint must not hide a smaller crossing witness.  The
  // pair (6,10) crosses (1,8), even though the already selected (2,100)
  // interval extends beyond its right endpoint.
  constexpr uint long_length = 101;
  SparseFloatMatrix long_posterior;
  long_posterior.assign(long_length, long_length);
  long_posterior.set(2, 100, 0.95f);
  long_posterior.set(1, 8, 0.94f);
  long_posterior.set(6, 10, 0.85f);

  LinearIPknot long_decoder(VF{0.9f, 0.8f, 0.8f}, 64);
  VU long_structure;
  std::string long_brackets;
  const float long_score =
      long_decoder.decode(long_posterior, long_structure, long_brackets);
  assert(long_structure[2] == 100 && long_structure[100] == -1u);
  assert(long_structure[1] == 8 && long_structure[8] == -1u);
  assert(long_structure[6] == 10 && long_structure[10] == -1u);
  check_structure_brackets(long_structure, long_brackets);
  const float expected_long_score = (0.95f - 0.9f) +
                                    (0.94f - 0.8f) +
                                    (0.85f - 0.8f);
  assert(std::isfinite(long_score));
  assert(std::fabs(long_score - expected_long_score) < 1e-4f);

  // A raw sparse row can contain many alternatives for one right endpoint.
  // The support union and deduplication must remain bounded and return a
  // valid matching even when all alternatives share that endpoint.
  constexpr uint crowded_length = 128;
  SparseFloatMatrix crowded_posterior;
  crowded_posterior.assign(crowded_length, crowded_length);
  for (uint left = 0; left + 2 < crowded_length; ++left)
    crowded_posterior.set(left, crowded_length - 1,
                          0.2f + 0.005f * left);
  LinearIPknot crowded_decoder(VF{0.1f}, 4);
  VU crowded_structure;
  const float crowded_score =
      crowded_decoder.decode(crowded_posterior, crowded_structure,
                             brackets);
  assert(crowded_structure[124] == crowded_length - 1);
  assert(std::isfinite(crowded_score));
}
