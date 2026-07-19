#include "linearalign.h"
#include "linfold_wrapper.h"

#include <exception>
#include <cmath>
#include <iostream>
#include <string>

// The production definitions live in fold.cpp/align.cpp together with the
// non-linear model implementations.  Keep this focused test independent of
// those optional libraries while satisfying the two base-class vtables.
void Fold::Model::calculate(const std::vector<Fasta>& fa, std::vector<BP>& bp)
{
  bp.resize(fa.size());
  for (size_t i = 0; i < fa.size(); ++i)
    calculate(fa[i].seq(), bp[i]);
}

void Align::Model::calculate(
    const std::vector<Fasta>& fa, std::vector<std::vector<MP>>& mp)
{
  mp.resize(fa.size(), std::vector<MP>(fa.size()));
  for (size_t i = 0; i < fa.size(); ++i)
    for (size_t j = i + 1; j < fa.size(); ++j)
      calculate(fa[i].seq(), fa[j].seq(), mp[i][j]);
}

namespace
{
template <typename Function>
bool throws_with(Function&& function, const std::string& expected)
{
  try {
    function();
  } catch (const std::exception& e) {
    if (std::string(e.what()).find(expected) != std::string::npos)
      return true;
    std::cerr << "unexpected exception: " << e.what() << '\n';
    return false;
  }
  std::cerr << "expected an exception containing: " << expected << '\n';
  return false;
}
}

int main()
{
  LinFoldWrapper fold(0.01f, LinFoldWrapper::ModelType::LPC, 10);
  BP bp;
  if (!throws_with(
          [&] { fold.calculate("ACGU", "???", bp); },
          "LinearPartition LPC failed"))
    return 1;

  LinearAlign align(0.01f, 10);
  align.setHMMParameters(nullptr, nullptr);
  MP mp;
  if (!throws_with(
          [&] { align.calculate("ACGU", "ACGU", mp); },
          "LinearAlign failed"))
    return 1;

  const std::string sequence = "GGGAAACCC";
  for (const auto model : {LinFoldWrapper::ModelType::LPV,
                           LinFoldWrapper::ModelType::LPC}) {
    LinFoldWrapper working_fold(0.01f, model, 100);
    working_fold.calculate(sequence, bp);
    size_t pair_count = 0;
    std::vector<float> marginal(sequence.size(), 0.0f);
    if (bp.size() != sequence.size())
      return 1;
    for (size_t i = 0; i < bp.size(); ++i) {
      for (const auto& [j, probability] : bp[i]) {
        if (j <= i || j >= sequence.size() || !std::isfinite(probability) ||
            probability <= 0.0f || probability > 1.0f)
          return 1;
        marginal[i] += probability;
        marginal[j] += probability;
        ++pair_count;
      }
    }
    if (pair_count == 0)
      return 1;
    for (const float probability : marginal)
      if (probability > 1.0001f)
        return 1;
  }

  LinearAlign working_align(0.01f, 100);
  working_align.calculate(sequence, "GGGAAUCCC", mp);
  size_t match_count = 0;
  for (size_t i = 0; i < mp.size(); ++i) {
    for (const auto& [j, probability] : mp[i]) {
      if (j >= sequence.size() || !std::isfinite(probability) ||
          probability <= 0.0f || probability > 1.0f)
        return 1;
      ++match_count;
    }
  }
  if (match_count == 0)
    return 1;

  return 0;
}
