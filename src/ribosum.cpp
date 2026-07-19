#include "ribosum.h"

#include <cctype>
#include <stdexcept>
#include <utility>

namespace
{
// RIBOSUM85-60 base-pair substitution scores from the official Infernal
// matrix distribution (Eddy/Rivas Lab).  Only the lower triangle is used;
// the substitution matrix is symmetric.
constexpr std::array<std::array<float, 16>, 16> PAIR_SCORES = {{
  {{-2.488349f}},
  {{-7.042094f, -2.108879f}},
  {{-8.238017f, -8.895066f, -0.803423f}},
  {{-4.317463f, -2.038769f, -5.133726f, 4.492700f}},
  {{-8.842033f, -9.372576f, -10.407814f, -5.564446f, -5.125312f}},
  {{-14.373559f, -9.081234f, -14.496253f, -6.705747f, -10.448895f, -3.593042f}},
  {{-4.678511f, -5.856881f, -4.567757f, 1.673203f, -3.567043f, -5.704751f, 5.360799f}},
  {{-12.640723f, -10.446549f, -10.140939f, -5.173809f, -8.485494f, -5.771068f, -4.963137f, -2.275363f}},
  {{-6.858807f, -9.728310f, -8.609094f, -5.328293f, -7.981155f, -12.429263f, -5.996523f, -7.708326f, -1.046177f}},
  {{-5.030739f, -3.812753f, -5.770632f, 2.704820f, -5.949981f, -3.701762f, 2.112560f, -5.842817f, -4.876907f, 5.616325f}},
  {{-8.393424f, -11.052794f, -5.383659f, -5.607690f, -11.357831f, -12.578234f, -4.664596f, -13.694318f, -8.670446f, -4.130013f, -1.975120f}},
  {{-5.835310f, -4.720463f, -6.596424f, 0.593694f, -7.929741f, -7.873675f, -0.270488f, -5.612507f, -6.094950f, 1.205393f, -5.767859f, 3.468458f}},
  {{-4.006673f, -5.324754f, -5.430794f, 1.608648f, -2.415577f, -6.876950f, 2.748442f, -4.717212f, -5.847454f, 1.596571f, -5.746373f, -0.565905f, 4.967781f}},
  {{-11.323969f, -8.665760f, -8.871719f, -4.812943f, -7.084955f, -7.402874f, -4.909116f, -3.834964f, -6.628842f, -4.485037f, -12.010797f, -5.302508f, -2.981881f, -3.208242f}},
  {{-6.161755f, -6.925209f, -5.941616f, -0.505944f, -5.629131f, -8.412817f, 1.319530f, -7.352963f, -7.551171f, -0.077872f, -4.273159f, -2.088271f, 1.136473f, -4.762134f, 3.364318f}},
  {{-9.048210f, -7.827377f, -11.073954f, -2.979395f, -8.393641f, -5.406597f, -3.671990f, -5.212174f, -11.540895f, -3.899205f, -10.786610f, -4.444826f, -3.387513f, -5.975262f, -4.277820f, -0.018299f}}
}};
}

float Ribosum85_60::pair_score(unsigned first_pair, unsigned second_pair)
{
  if (first_pair >= 16 || second_pair >= 16)
    throw std::out_of_range("RIBOSUM pair index is out of range");
  if (first_pair < second_pair)
    std::swap(first_pair, second_pair);
  return PAIR_SCORES[first_pair][second_pair];
}

int Ribosum85_60::nucleotide_code(char nucleotide)
{
  switch (std::toupper(static_cast<unsigned char>(nucleotide))) {
  case 'A': return 0;
  case 'C': return 1;
  case 'G': return 2;
  case 'T':
  case 'U': return 3;
  default: return -1;
  }
}

RibosumProfile::RibosumProfile(
    const ALN& alignment, const std::vector<Fasta>& sequences)
  : profile_size_(static_cast<unsigned>(alignment.size()))
{
  if (alignment.empty())
    throw std::invalid_argument("RIBOSUM profile alignment is empty");
  const size_t length = alignment.front().second.size();
  columns_.assign(length, std::vector<std::int8_t>(alignment.size(), -1));

  for (size_t profile_row = 0; profile_row < alignment.size(); ++profile_row) {
    const auto& [sequence_index, present] = alignment[profile_row];
    if (sequence_index >= sequences.size() || present.size() != length)
      throw std::invalid_argument("invalid alignment in RIBOSUM profile");
    const std::string& sequence = sequences[sequence_index].seq();
    size_t sequence_position = 0;
    for (size_t column = 0; column < length; ++column) {
      if (!present[column])
        continue;
      if (sequence_position >= sequence.size())
        throw std::invalid_argument("alignment exceeds sequence in RIBOSUM profile");
      columns_[column][profile_row] = static_cast<std::int8_t>(
          Ribosum85_60::nucleotide_code(sequence[sequence_position++]));
    }
    if (sequence_position != sequence.size())
      throw std::invalid_argument("alignment does not consume sequence in RIBOSUM profile");
  }
}

std::uint64_t RibosumProfile::key(unsigned i, unsigned j)
{
  return (static_cast<std::uint64_t>(i) << 32) | j;
}

RibosumProfile::PairDistribution
RibosumProfile::pair_distribution(unsigned i, unsigned j) const
{
  if (i >= columns_.size() || j >= columns_.size())
    throw std::out_of_range("RIBOSUM profile column is out of range");
  const std::uint64_t cache_key = key(i, j);
  const auto found = cache_.find(cache_key);
  if (found != cache_.end())
    return found->second;

  PairDistribution distribution{};
  for (unsigned row = 0; row < profile_size_; ++row) {
    const int first = columns_[i][row];
    const int second = columns_[j][row];
    if (first >= 0 && second >= 0)
      distribution[4 * first + second] += 1.0f / profile_size_;
  }
  cache_.emplace(cache_key, distribution);
  return distribution;
}

float RibosumProfile::pair_score(
    unsigned i, unsigned j, const RibosumProfile& other,
    unsigned k, unsigned l) const
{
  const PairDistribution first = pair_distribution(i, j);
  const PairDistribution second = other.pair_distribution(k, l);
  double score = 0.0;
  for (unsigned x = 0; x < 16; ++x)
    for (unsigned y = 0; y < 16; ++y)
      score += static_cast<double>(first[x]) * second[y] *
               Ribosum85_60::pair_score(x, y);
  return static_cast<float>(score);
}
