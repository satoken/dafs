#include "nussinov.h"
#include "dd_recovery.h"

#include <cmath>
#include <functional>
#include <iostream>
#include <map>
#include <random>
#include <stdexcept>

const uint Fold::Decoder::n_support_brackets = 30;
const char* Fold::Decoder::left_brackets = "([{<ABCDEFGHIJKLMNOPQRSTUVWXYZ";
const char* Fold::Decoder::right_brackets = ")]}>abcdefghijklmnopqrstuvwxyz";

namespace {
void require(bool value, const std::string& message)
{
  if (!value) throw std::runtime_error(message);
}

bool close(double a, double b)
{
  return std::fabs(a-b) <= 3e-5 * std::max(1.0, std::fabs(b));
}

SparseFloatMatrix to_sparse(const VVF& dense)
{
  SparseFloatMatrix sparse;
  sparse.assign(dense.size(), dense.size());
  for (uint i = 0; i < dense.size(); ++i)
    for (uint j = 0; j < dense.size(); ++j)
      if (dense[i][j] != 0.0f) sparse.set(i, j, dense[i][j]);
  return sparse;
}

using Pairs = std::vector<std::pair<uint, uint>>;
using Configurations = std::vector<Pairs>;

// Enumerate all noncrossing matchings before independently checking stack
// support. This oracle does not use the production two-state recurrence.
Configurations enumerate(int length, int span)
{
  std::map<std::pair<int, int>, Configurations> cache;
  std::function<const Configurations&(int, int)> visit = [&](int i, int j)
      -> const Configurations& {
    const auto key = std::make_pair(i, j);
    const auto found = cache.find(key);
    if (found != cache.end()) return found->second;
    Configurations out;
    if (i > j) out.emplace_back();
    else {
      out = visit(i+1, j);
      for (int k = i+span; k <= j; ++k)
        for (const auto& inside : visit(i+1, k-1))
          for (const auto& suffix : visit(k+1, j)) {
            Pairs pairs{{static_cast<uint>(i), static_cast<uint>(k)}};
            pairs.insert(pairs.end(), inside.begin(), inside.end());
            pairs.insert(pairs.end(), suffix.begin(), suffix.end());
            out.push_back(std::move(pairs));
          }
    }
    return cache.emplace(key, std::move(out)).first->second;
  };
  return visit(0, length-1);
}

double oracle(const Configurations& configurations, float w,
              const VVF& p, const VVF& q, float threshold)
{
  double best = 0.0;
  for (const auto& pairs : configurations) {
    VU structure(p.size(), -1u);
    double score = 0.0;
    bool supported = true;
    for (const auto [i, j] : pairs) {
      structure[i] = j;
      supported &= p[i][j] != 0.0f || q[i][j] != 0.0f;
      score += static_cast<double>(w * (p[i][j] - threshold) - q[i][j]);
    }
    if (supported && DDRecovery::structure(structure, false))
      best = std::max(best, score);
  }
  return best;
}

void validate(float score, const VU& structure, float w,
              const VVF& p, const VVF& q, float threshold)
{
  require(structure.size() == p.size(), "wrong structure length");
  require(DDRecovery::structure(structure, false), "isolated or crossing output pair");
  double rescored = 0.0;
  for (uint i = 0; i < structure.size(); ++i) {
    const uint j = structure[i];
    if (j == -1u) continue;
    require(p[i][j] != 0.0f || q[i][j] != 0.0f, "selected missing candidate");
    rescored += static_cast<double>(w * (p[i][j] - threshold) - q[i][j]);
  }
  require(close(score, rescored), "traceback score differs from returned score");
}

void fixtures()
{
  constexpr uint length = 10;
  VVF p(length, VF(length)), q = p;
  p[0][9] = 0.95f;
  Nussinov exact(0.2f), allowed_exact(0.2f, true);
  SparseNussinov sparse(0.2f), allowed_sparse(0.2f, true);
  LinearNussinov linear(0.2f, length), allowed_linear(0.2f, length, true);
  VU ss;
  std::string brackets;
  for (Fold::Decoder* decoder : std::initializer_list<Fold::Decoder*>{&exact, &sparse, &linear}) {
    require(decoder->decode(p, ss, brackets) == 0.0f, "default allows lone pair");
    require(brackets == "..........", "lone pair rendered by default");
  }
  for (Fold::Decoder* decoder : std::initializer_list<Fold::Decoder*>{&allowed_exact, &allowed_sparse, &allowed_linear}) {
    require(close(decoder->decode(p, ss, brackets), 0.75), "allow switch lost lone pair");
    require(brackets == "(........)", "allowed lone pair traceback");
  }
  // A below-threshold neighbor must be retained if the combined helix wins.
  p[1][8] = 0.1f;
  SparseFloatMatrix ps, qs, bonus;
  ps = to_sparse(p); qs = to_sparse(q); bonus.assign(length, length);
  for (Fold::Decoder* decoder : std::initializer_list<Fold::Decoder*>{&exact, &sparse, &linear}) {
    require(close(decoder->decode(p, ss, brackets), 0.65), "negative-scoring stack lost");
    require(brackets == "((......))", "stacked helix traceback");
    require(close(decoder->decode(ps, bonus, ss, brackets), 0.65), "bonus overload lost stack");
  }
  // A multiplier-only neighboring pair is also a valid stack candidate.
  p[1][8] = 0.0f; q[1][8] = -0.1f;
  ps = to_sparse(p); qs = to_sparse(q);
  require(close(exact.decode(1.0f, p, q, ss), 0.65), "dense q-only pair lost");
  require(close(sparse.decode(1.0f, p, qs, ss), 0.65), "sparse q-only pair lost");
  require(close(linear.decode(1.0f, ps, qs, ss), 0.65), "linear q-only pair lost");
  bonus.set(1, 8, 0.1f);
  for (Fold::Decoder* decoder : std::initializer_list<Fold::Decoder*>{&exact, &sparse, &linear})
    require(close(decoder->decode(ps, bonus, ss, brackets), 0.65), "bonus-only stack lost");

  // Reuse the paired-state workspace across empty, short, and longer inputs.
  for (uint size : {0u, 1u, 4u, 10u, 18u, 10u}) {
    SparseFloatMatrix probability;
    probability.assign(size, size);
    if (size >= 10) {
      probability.set(0, size-1, 0.95f);
      probability.set(1, size-2, 0.1f);
    }
    const float value = linear.decode(probability, ss, brackets);
    require(close(value, size >= 10 ? 0.65 : 0.0), "workspace reuse changed score");
    require(DDRecovery::structure(ss, false), "workspace reuse invalid traceback");
  }
}

void randomized()
{
  constexpr uint length = 12;
  constexpr float threshold = 0.3f, w = 1.7f;
  const auto normal = enumerate(length, 3), cached = enumerate(length, 2);
  std::mt19937 generator(89213);
  std::uniform_real_distribution<float> probability(0.01f, 1.0f), multiplier(-0.7f, 0.7f);
  std::bernoulli_distribution present(0.7);
  size_t pruned = 0;
  for (uint trial = 0; trial < 120; ++trial) {
    VVF p(length, VF(length)), q = p;
    VVU support(length);
    const bool span_two = trial % 2;
    for (uint i = 0; i < length; ++i)
      for (uint j = i + (span_two ? 2 : 3); j < length; ++j) {
        if (present(generator)) p[i][j] = probability(generator);
        if (present(generator)) q[i][j] = multiplier(generator);
        if (p[i][j] != 0 || q[i][j] != 0) support[j].push_back(i);
      }
    SparseFloatMatrix ps, qs;
    ps = to_sparse(p); qs = to_sparse(q);
    const double best = oracle(span_two ? cached : normal, w, p, q, threshold);
    VU ss;
    if (!span_two) {
      Nussinov exact(threshold);
      SparseNussinov sparse(threshold);
      require(close(exact.decode(w, p, q, ss), best), "dense DP differs from enumeration");
      validate(best, ss, w, p, q, threshold);
      require(close(sparse.decode(w, p, q, ss), best), "sparse DP differs from enumeration");
      validate(best, ss, w, p, q, threshold);
      require(close(sparse.decode(w, p, qs, ss), best), "mixed DP differs from enumeration");
      LinearNussinov wide(threshold, length);
      require(close(wide.decode(w, p, q, ss), best), "wide dense beam differs from enumeration");
      require(close(wide.decode(w, p, qs, ss), best), "wide mixed beam differs from enumeration");
      require(close(wide.decode(w, ps, q, ss), best), "wide sparse/dense beam differs");
      require(close(wide.decode(w, ps, qs, ss), best), "wide sparse beam differs from enumeration");
      validate(best, ss, w, p, q, threshold);
    }
    // Malformed cached entries must be rejected before accessing matrices.
    for (uint right = 0; right < length; ++right) {
      support[right].push_back(right);
      support[right].push_back(length+right);
      if (right > 0) support[right].push_back(right-1);
    }
    for (uint beam : {1u, 2u, 5u, length}) {
      LinearNussinov decoder(threshold, beam);
      const auto result = decoder.decode_certified(w, ps, qs, support, ss);
      validate(result.score, ss, w, p, q, threshold);
      require(result.score <= best + 3e-5, "beam exceeds exact optimum");
      require(result.upper_bound + 3e-5 >= best, "certificate below exact optimum");
      require(result.upper_bound <= result.additive_upper_bound + 3e-5, "invalid additive bound");
      require(result.score <= result.upper_bound + 3e-5, "beam exceeds certificate");
      if (beam == length) require(close(result.score, best), "wide cached beam differs from enumeration");
      pruned += result.pruned_states;
    }
  }
  require(pruned > 0, "random tests did not exercise beam pruning");
}

void repair()
{
  VU x(8, -1u), y(9, -1u), z{0, 2, 3, 4, 5, 6, 7, 8};
  x[0] = 7; x[1] = 6;
  y[0] = 8; y[2] = 7;
  require(DDRecovery::coupled(x, y, z), "invalid repair fixture");
  require(!DDRecovery::coupled(x, y, z, false), "gap-broken stack accepted");
  DDRecovery::remove_isolated(x, y, z);
  require(DDRecovery::coupled(x, y, z, false), "repair left an isolated pair");
  require(x[0] == -1u && x[1] == -1u, "repair missed a cascade");
  x[1] = 6; x[2] = 5; y[2] = 7; y[3] = 6;
  DDRecovery::remove_isolated(x, y, z);
  require(x[1] == 6 && x[2] == 5 && DDRecovery::coupled(x, y, z, false),
          "repair removed a valid shared helix");
}

void long_sparse()
{
  constexpr uint length = 20000;
  SparseFloatMatrix p, q;
  p.assign(length, length); q.assign(length, length);
  for (uint i = 0; i < length; i += 8) {
    p.set(i, i+7, 0.9f);
    p.set(i+1, i+6, 0.9f);
  }
  LinearNussinov decoder(0.2f, 4);
  VU ss;
  const auto result = decoder.decode_certified(1.0f, p, q, ss);
  require(close(result.score, 0.7 * (length/4)), "long sparse helices lost");
  require(result.upper_bound + 0.01 >= result.score &&
          DDRecovery::structure(ss, false), "invalid long sparse result");
}
} // namespace

int main()
{
  try {
    fixtures(); randomized(); repair(); long_sparse();
    std::cout << "Isolated-pair switches, exhaustive oracles, beam certificates, and repair passed\n";
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
