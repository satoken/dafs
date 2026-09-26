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

  IPknot dense_final(VF{0.2f});
  IPknot sparse_final(VF{0.2f});
  std::string dense_brackets, sparse_brackets;
  dense_final.decode(dense_probability, dense_result, dense_brackets);
  sparse_final.decode(sparse_probability, sparse_result, sparse_brackets);
  assert(dense_result == sparse_result);
  assert(dense_brackets == sparse_brackets);
}
