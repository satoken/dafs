// Direct LinearAlifold profile probe for round4.
//
// The program prints one JSON object per case.  LogZ is recorded as both a
// decimal value and its IEEE-754 bits; every returned BPP is recorded the same
// way.  The wrapper is deliberately reused across the length-boundary and
// changed-profile cases so stale dense caches are observable.

#include "linfold_wrapper.h"

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// The standalone probe links the same wrapper object as the CTest oracle.
// Supply the non-linear base overload here so no unrelated DAFS folding
// objects (RNAfold/CONTRAfold) enter the probe link.
void Fold::Model::calculate(const std::vector<Fasta>& fa, std::vector<BP>& bp)
{
  bp.resize(fa.size());
  for (size_t index = 0; index < fa.size(); ++index)
    calculate(fa[index].seq(), bp[index]);
}

namespace {

struct Case {
  std::string name;
  std::vector<std::string> aligned;
  std::string constraint;
};

uint64_t bits(double value)
{
  uint64_t result = 0;
  std::memcpy(&result, &value, sizeof(result));
  return result;
}

uint32_t bits(float value)
{
  uint32_t result = 0;
  std::memcpy(&result, &value, sizeof(result));
  return result;
}

std::string hex64(uint64_t value)
{
  std::ostringstream output;
  output << std::hex << std::setw(16) << std::setfill('0') << value;
  return output.str();
}

std::string hex32(uint32_t value)
{
  std::ostringstream output;
  output << std::hex << std::setw(8) << std::setfill('0') << value;
  return output.str();
}

char complement(char base)
{
  switch (base) {
  case 'A': return 'U';
  case 'C': return 'G';
  case 'G': return 'C';
  case 'U': return 'A';
  default: return 'A';
  }
}

std::vector<std::string> profile_rows(size_t length, unsigned variant)
{
  std::vector<std::string> rows(3, std::string(length, 'A'));
  static constexpr char alphabet[] = "ACGU";
  const size_t stem = std::min<size_t>(8, length / 2);
  for (size_t row = 0; row < rows.size(); ++row) {
    for (size_t column = 0; column < length; ++column)
      rows[row][column] = alphabet[(column * 7 + row * 3 + variant) % 4];
    for (size_t offset = 0; offset < stem; ++offset) {
      const char left = alphabet[(offset + row + variant) % 4];
      rows[row][offset] = left;
      rows[row][length - 1 - offset] = complement(left);
    }
  }
  // Keep the profile genuinely covarying while preserving a supported pair
  // in both directions for the first and last stem columns.
  if (length >= 12) {
    rows[1][0] = 'G';
    rows[1][length - 1] = 'C';
    rows[2][0] = 'C';
    rows[2][length - 1] = 'G';
  }
  return rows;
}

std::vector<std::string> gapped_rows(size_t length, unsigned variant)
{
  auto rows = profile_rows(length, variant);
  for (size_t row = 0; row < rows.size(); ++row) {
    const size_t first_gap = 11 + row;
    const size_t second_gap = 37 + 2 * row;
    if (first_gap < length) rows[row][first_gap] = '-';
    if (second_gap < length) rows[row][second_gap] = '-';
  }
  return rows;
}

Case make_case(std::string name, std::vector<std::string> rows,
               std::string constraint = {})
{
  const size_t length = rows.front().size();
  for (const auto& row : rows)
    if (row.size() != length)
      throw std::runtime_error("profile rows have different lengths");
  return {std::move(name), std::move(rows), std::move(constraint)};
}

ALN alignment_for(const std::vector<std::string>& aligned,
                  std::vector<Fasta>& sequences)
{
  ALN alignment;
  sequences.clear();
  sequences.reserve(aligned.size());
  alignment.reserve(aligned.size());
  for (size_t row = 0; row < aligned.size(); ++row) {
    std::string sequence;
    std::vector<bool> present;
    sequence.reserve(aligned[row].size());
    present.reserve(aligned[row].size());
    for (char base : aligned[row]) {
      const bool present_here = base != '-';
      present.push_back(present_here);
      if (present_here) sequence.push_back(base);
    }
    sequences.emplace_back("row" + std::to_string(row), sequence);
    alignment.emplace_back(static_cast<unsigned>(row), std::move(present));
  }
  return alignment;
}

void print_case(const std::string& run_name, const Case& test,
                LinFoldWrapper& wrapper)
{
  std::vector<Fasta> sequences;
  const ALN alignment = alignment_for(test.aligned, sequences);
  BP bpp;
  const double logz = wrapper.calculate_profile(
      alignment, sequences, bpp, test.constraint);

  struct PairValue { unsigned i, j; float probability; };
  std::vector<PairValue> pairs;
  for (size_t i = 0; i < bpp.size(); ++i)
    for (const auto& [j, probability] : bpp[i])
      pairs.push_back({static_cast<unsigned>(i), j, probability});
  std::sort(pairs.begin(), pairs.end(), [](const PairValue& left,
                                           const PairValue& right) {
    return left.i != right.i ? left.i < right.i : left.j < right.j;
  });

  std::cout << "{\"run\":\"" << run_name
            << "\",\"case\":\"" << test.name
            << "\",\"length\":" << test.aligned.front().size()
            << ",\"rows\":" << test.aligned.size()
            << ",\"constraint\":\"" << (test.constraint.empty() ? "none" : test.constraint)
            << "\",\"logz\":" << std::setprecision(17) << logz
            << ",\"logz_bits\":\"" << hex64(bits(logz)) << "\""
            << ",\"pair_count\":" << pairs.size() << ",\"pairs\":[";
  for (size_t index = 0; index < pairs.size(); ++index) {
    if (index) std::cout << ',';
    const auto& pair = pairs[index];
    std::cout << "{\"i\":" << pair.i << ",\"j\":" << pair.j
              << ",\"p\":" << std::setprecision(9) << pair.probability
              << ",\"p_bits\":\"" << hex32(bits(pair.probability)) << "\"}";
  }
  std::cout << "]}\n";
}

template <typename Options>
void allow_canonical_pairs_probe(Options& options)
{
  options.set_allowed_pair('a', 'u');
  options.set_allowed_pair('c', 'g');
  options.set_allowed_pair('g', 'u');
}

std::string sequence_for(size_t length, unsigned variant)
{
  static constexpr char alphabet[] = "ACGU";
  std::string sequence(length, 'A');
  for (size_t index = 0; index < length; ++index)
    sequence[index] = alphabet[(index * 7 + variant * 3) % 4];
  const size_t stem = std::min<size_t>(8, length / 2);
  for (size_t offset = 0; offset < stem; ++offset) {
    const char left = alphabet[(offset + variant) % 4];
    sequence[offset] = left;
    sequence[length - 1 - offset] = complement(left);
  }
  return sequence;
}

template <typename NearestNeighbor>
void print_sequence_case(const std::string& run_name,
                         const std::string& case_name,
                         const std::string& sequence)
{
  using Engine = LinFold<NearestNeighbor>;
  auto params = std::make_unique<NearestNeighbor>(sequence);
  Engine engine(std::move(params));
  typename Engine::Options options;
  options.beam_size(100);
  allow_canonical_pairs_probe(options);
  const auto logz_score = engine.compute_inside(sequence, options);
  engine.compute_outside(sequence, options);
  const auto bpp = engine.compute_basepairing_probabilities(sequence, options);

  struct PairValue { unsigned i, j; float probability; };
  std::vector<PairValue> pairs;
  for (size_t i = 0; i < bpp.size(); ++i)
    for (const auto& [j, probability] : bpp[i])
      pairs.push_back({static_cast<unsigned>(i), j, probability});
  std::sort(pairs.begin(), pairs.end(), [](const PairValue& left,
                                           const PairValue& right) {
    return left.i != right.i ? left.i < right.i : left.j < right.j;
  });

  const double logz = static_cast<double>(logz_score);
  std::cout << "{\"run\":\"" << run_name
            << "\",\"case\":\"" << case_name
            << "\",\"length\":" << sequence.size()
            << ",\"rows\":1,\"constraint\":\"none\""
            << ",\"logz\":" << std::setprecision(17) << logz
            << ",\"logz_bits\":\"" << hex64(bits(logz)) << "\""
            << ",\"pair_count\":" << pairs.size() << ",\"pairs\":[";
  for (size_t index = 0; index < pairs.size(); ++index) {
    if (index) std::cout << ',';
    const auto& pair = pairs[index];
    std::cout << "{\"i\":" << pair.i << ",\"j\":" << pair.j
              << ",\"p\":" << std::setprecision(9) << pair.probability
              << ",\"p_bits\":\"" << hex32(bits(pair.probability)) << "\"}";
  }
  std::cout << "]}\n";
}

std::vector<Case> boundary_cases()
{
  std::vector<Case> cases;
  cases.push_back(make_case("stem24", profile_rows(24, 0)));
  cases.push_back(make_case("boundary63", profile_rows(63, 0)));
  cases.push_back(make_case("boundary64", profile_rows(64, 0)));
  cases.push_back(make_case("boundary1024", profile_rows(1024, 0)));
  cases.push_back(make_case("boundary1025", profile_rows(1025, 0)));
  // Same length, new sequence/profile: exercises cache invalidation on a
  // reused LinFoldWrapper and a different pair eligibility pattern.
  cases.push_back(make_case("reused64_variant", profile_rows(64, 1)));
  cases.push_back(make_case("reused64_gapped", gapped_rows(64, 2)));
  cases.push_back(make_case("reused64_all_unpaired", profile_rows(64, 3),
                            std::string(64, '.')));
  return cases;
}

}  // namespace

int main()
{
  // Profile folding is independent of LPV/LPC; LPV selects the same
  // ProfileNearestNeighbor path used by LinearAlifold.
  const auto cases = boundary_cases();
  LinFoldWrapper cutoff0(0.0f, LinFoldWrapper::ModelType::LPV, 100);
  for (const auto& test : cases) print_case("lpv_cutoff0", test, cutoff0);

  // Run selected cases at all production probability cutoffs.  The Python
  // comparator reports any pair crossing a cutoff between builds.
  for (const float cutoff : {0.001f, 0.01f, 0.05f, 0.1f}) {
    LinFoldWrapper wrapper(cutoff, LinFoldWrapper::ModelType::LPV, 100);
    for (const auto& test : cases)
      if (test.name == "stem24" || test.name == "boundary64" ||
          test.name == "reused64_gapped") {
        print_case("lpv_cutoff" + std::to_string(cutoff), test, wrapper);
    }
  }

  // Direct one-sequence LPC/LPV checks exercise the non-profile path and
  // record the complete returned BPP map.  Length 63 stays on the bounded
  // fallback path; 64 and 80 exercise dense storage after the threshold.
  for (const size_t length : {size_t(63), size_t(64), size_t(80)}) {
    print_sequence_case<TurnerNearestNeighbor>(
        "direct_lpv", "seq" + std::to_string(length),
        sequence_for(length, 0));
    print_sequence_case<CONTRAfoldNearestNeighbor>(
        "direct_lpc", "seq" + std::to_string(length),
        sequence_for(length, 1));
  }
}
