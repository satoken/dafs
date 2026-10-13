#include "nussinov.h"
#include "dd_recovery.h"

#include <cmath>
#include <iostream>
#include <random>
#include <stdexcept>

const uint Fold::Decoder::n_support_brackets = 30;
const char* Fold::Decoder::left_brackets = "([{<ABCDEFGHIJKLMNOPQRSTUVWXYZ";
const char* Fold::Decoder::right_brackets = ")]}>abcdefghijklmnopqrstuvwxyz";

namespace {
void require(bool value, const char* message)
{
  if (!value) throw std::runtime_error(message);
}

bool close(double a, double b)
{
  return std::fabs(a-b) <= 3e-5 * std::max(1.0, std::fabs(b));
}

SparseFloatMatrix sparse(const VVF& matrix)
{
  SparseFloatMatrix result;
  result.assign(matrix.size(), matrix.size());
  for (uint i = 0; i < matrix.size(); ++i)
    for (uint j = 0; j < matrix.size(); ++j)
      if (matrix[i][j] != 0.0f) result.set(i, j, matrix[i][j]);
  return result;
}

void validate(double actual, double expected, const VU& ss,
              const VVF& p, const VVF& q, float w, float th, bool allow)
{
  require(close(actual, expected), "sparsified score differs from independent dense DP");
  require(ss.size() == p.size() && DDRecovery::structure(ss, allow),
          "invalid sparsified traceback");
  double score = 0;
  for (uint i = 0; i < ss.size(); ++i) if (ss[i] != -1u) {
    const uint j = ss[i];
    require(j-i >= 3, "IPknot minimum span replaced DAFS minimum span");
    if (!allow) require(p[i][j] != 0 || q[i][j] != 0, "artificial stack partner");
    score += w * (p[i][j]-th) - q[i][j];
  }
  require(close(actual, score), "sparsified traceback has wrong score");
}

void randomized()
{
  std::mt19937 rng(38433200);
  std::uniform_real_distribution<float> probability(0.01f, 1.0f), multiplier(-1, 1);
  for (bool allow : {true, false}) {
    SparseNussinov decoder(0.3f, allow);
    Nussinov reference(0.3f, allow);
    Fold::Decoder& polymorphic = decoder;
    for (int trial = 0; trial < 100; ++trial) {
      const bool very_sparse = trial%5 == 0;
      const uint n = very_sparse ? 64 : 6 + rng()%25;
      VVF p(n, VF(n)), q = p, bonus = p;
      for (uint i = 0; i < n; ++i) for (uint j = i+1; j < n; ++j) {
        if (very_sparse ? rng()%200 == 0 : rng()%3 != 0) p[i][j] = probability(rng);
        if (very_sparse ? rng()%200 == 0 : rng()%3 != 0) q[i][j] = multiplier(rng);
        if (very_sparse ? rng()%200 == 0 : rng()%3 != 0) bonus[i][j] = multiplier(rng);
      }
      const float w = trial%3 == 0 ? 0.0f : 1.7f;
      const auto ps = sparse(p), qs = sparse(q), bs = sparse(bonus);
      VU ss, expected_ss;
      const float expected = reference.decode(w, p, q, expected_ss);
      validate(decoder.decode(w, p, q, ss), expected, ss, p, q, w, 0.3f, allow);
      validate(decoder.decode(w, p, qs, ss), expected, ss, p, q, w, 0.3f, allow);
      validate(polymorphic.decode(w, ps, q, ss), expected, ss, p, q, w, 0.3f, allow);
      validate(polymorphic.decode(w, ps, qs, ss), expected, ss, p, q, w, 0.3f, allow);
      std::string brackets, expected_brackets;
      const float final_expected = reference.decode(p, expected_ss, expected_brackets);
      const VVF zero(n, VF(n));
      validate(decoder.decode(p, ss, brackets), final_expected, ss, p, zero, 1, 0.3f, allow);
      validate(polymorphic.decode(ps, ss, brackets), final_expected, ss, p, zero, 1, 0.3f, allow);
      const float bonus_expected = reference.decode(ps, bs, expected_ss, expected_brackets);
      for (auto& row : bonus) for (auto& value : row) value = -value;
      validate(decoder.decode(ps, bs, ss, brackets), bonus_expected, ss, p, bonus, 1, 0.3f, allow);
      require(brackets.size() == n, "wrong bracket length");
    }
  }
}

void candidates_and_reuse()
{
  constexpr uint length = 60;
  VVF p(length, VF(length)), q = p;
  for (uint i = 0; i < length; ++i)
    for (uint j = i+3; j < length; ++j) p[i][j] = 1;
  for (bool allow : {true, false}) {
    SparseNussinov decoder(0, allow);
    VU ss;
    require(decoder.sparsified_candidates() == 0, "uninitialized pruning count");
    for (int reuse = 0; reuse < 2; ++reuse) {
      const double score = decoder.decode(1, p, q, ss);
      validate(score, 29, ss, p, q, 1, 0, allow);
      require(decoder.sparsified_candidates() > 100, "closed candidates were not sparsified");
    }
    // Same graph, changed weights must not retain old paired states or traces.
    for (uint i = 0; i < length; ++i)
      for (uint j = i+3; j < length; ++j) q[i][j] = 2;
    require(decoder.decode(1, p, q, ss) == 0, "reused weights retained old optimum");
    require(decoder.sparsified_candidates() == 0, "reused pruning count leaked");
    q.assign(length, VF(length));
    for (uint n : {0u, 1u, 2u, 3u, 5u, length}) {
      SparseFloatMatrix empty;
      empty.assign(n, n);
      std::string brackets;
      require(decoder.decode(empty, ss, brackets) == 0 && brackets == std::string(n, '.'),
              "empty/short graph or length change retained old states");
      require(decoder.sparsified_candidates() == 0, "empty graph pruning count leaked");
    }
  }
  // Absent entries can be profitable with a negative threshold. Preserve
  // that historical search space only when isolated pairs are allowed.
  SparseFloatMatrix empty;
  empty.assign(8, 8);
  SparseNussinov allowed(-0.1f), stacked(-0.1f, false);
  VU ss;
  require(close(allowed.decode(1, empty, empty, ss), 0.3), "negative threshold lost absent pairs");
  require(stacked.decode(1, empty, empty, ss) == 0, "negative threshold fabricated stack support");
}

void long_sparse()
{
  constexpr uint length = 20000;
  SparseFloatMatrix p, q;
  p.assign(length, length); q.assign(length, length);
  p.set(0, length-1, 0.95f);
  p.set(1, length-2, 0.1f); // negative child must support its profitable parent
  for (bool allow : {true, false}) {
    SparseNussinov decoder(0.2f, allow);
    Fold::Decoder& polymorphic = decoder;
    VU ss;
    require(close(polymorphic.decode(1, p, q, ss), allow ? 0.75 : 0.65),
            "long sparse graph lost signed helix");
    require(DDRecovery::structure(ss, allow), "long sparse traceback is invalid");
    require(ss[0] == length-1 && ss[1] == (allow ? -1u : length-2),
            "long sparse trace omitted distant pairs");
  }
}
} // namespace

int main()
{
  try {
    randomized(); candidates_and_reuse(); long_sparse();
    std::cout << "Sparsified exact DP, dense/sparse adapters, signed weights, workspace reuse and long sparse input passed\n";
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
