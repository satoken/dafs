#include "nussinov.h"

#include <iostream>
#include <stdexcept>

const uint Fold::Decoder::n_support_brackets = 30;
const char* Fold::Decoder::left_brackets = "([{<ABCDEFGHIJKLMNOPQRSTUVWXYZ";
const char* Fold::Decoder::right_brackets = ")]}>abcdefghijklmnopqrstuvwxyz";

namespace {
void require(bool value, const char* message)
{
  if (!value) throw std::runtime_error(message);
}

SparseFloatMatrix sparse(const VVF& matrix)
{
  SparseFloatMatrix out;
  out.assign(matrix.size(), matrix.size());
  for (uint i = 0; i < matrix.size(); ++i)
    for (uint j = 0; j < matrix.size(); ++j)
      if (matrix[i][j] != 0) out.set(i, j, matrix[i][j]);
  return out;
}

void tied_crossing_pairs()
{
  // Old DAFS prefers skipping the left endpoint on an exact tie.
  VVF p(6, VF(6)), q = p;
  p[0][4] = p[1][5] = .9f;
  const auto ps = sparse(p), qs = sparse(q);
  DAFSNussinov decoder(.2f);
  Fold::Decoder& polymorphic = decoder;
  VU ss;
  const auto check = [&](float score) {
    require(score == .9f - .2f, "legacy tied score changed");
    require(ss == VU({-1u, 5, -1u, -1u, -1u, -1u}),
            "legacy skip-left tie choice changed");
  };
  check(decoder.decode(1, p, q, ss));
  check(decoder.decode(1, p, qs, ss));
  check(polymorphic.decode(1, ps, q, ss));
  check(polymorphic.decode(1, ps, qs, ss));
  std::string brackets;
  check(decoder.decode(p, ss, brackets));
  require(brackets == ".(...)", "legacy final traceback changed");
  check(polymorphic.decode(ps, ss, brackets));
  check(decoder.decode(ps, qs, ss, brackets));
}

void float_accumulation()
{
  // Both additions round away a unit under the old nested float DP.
  // Widening the chart to double changes the returned score and tie choice.
  VVF p(10, VF(10)), q = p;
  p[0][9] = p[1][8] = p[2][7] = 1;
  q[1][8] = -16777215.0f;
  const auto ps = sparse(p), qs = sparse(q);
  DAFSNussinov decoder(0);
  Fold::Decoder& polymorphic = decoder;
  VU ss;
  const auto check = [&](float score) {
    require(score == 16777216.0f, "legacy float accumulation widened");
    require(ss[0] == -1u && ss[1] == 8 && ss[2] == 7,
            "legacy float near-tie choice changed");
  };
  check(decoder.decode(1, p, q, ss));
  check(decoder.decode(1, p, qs, ss));
  check(polymorphic.decode(1, ps, q, ss));
  check(polymorphic.decode(1, ps, qs, ss));
  SparseFloatMatrix bonus = qs;
  bonus.set(1, 8, 16777215.0f);
  std::string brackets;
  check(decoder.decode(ps, bonus, ss, brackets));
}

void empty_and_short()
{
  for (bool allow : {true, false}) {
    DAFSNussinov decoder(.2f, allow);
    Fold::Decoder& polymorphic = decoder;
    for (uint n : {0u, 1u, 2u, 3u}) {
      VVF p(n, VF(n)), q = p;
      const auto ps = sparse(p), qs = sparse(q);
      VU ss(20, 3);
      std::string brackets(20, '(');
      require(decoder.decode(1, p, q, ss) == 0 && ss == VU(n, -1u),
              "empty/short DAFS weighted decode failed");
      require(polymorphic.decode(1, ps, qs, ss) == 0 && ss == VU(n, -1u),
              "empty/short DAFS sparse decode failed");
      require(polymorphic.decode(ps, ss, brackets) == 0 && brackets == std::string(n, '.'),
              "empty/short DAFS final decode failed");
      require(decoder.decode(ps, qs, ss, brackets) == 0 && brackets == std::string(n, '.'),
              "empty/short DAFS bonus decode failed");
    }
    require(decoder.used_workers() == 1 && decoder.parallel_steps() == 0,
            "legacy DAFS kernel is no longer serial");
  }
}
}

int main()
{
  try {
    tied_crossing_pairs();
    float_accumulation();
    empty_and_short();
    std::cout << "DAFS default float scores, tie choices and adapters passed\n";
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
