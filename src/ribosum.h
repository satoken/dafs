#ifndef DAFS_RIBOSUM_H
#define DAFS_RIBOSUM_H

#include <array>
#include <cstdint>
#include <unordered_map>
#include <vector>

#include "fa.h"
#include "typedefs.h"

namespace Ribosum85_60
{
// Ordered pair alphabet: AA AC AG AU CA CC CG CU GA GC GG GU UA UC UG UU.
float pair_score(unsigned first_pair, unsigned second_pair);
int nucleotide_code(char nucleotide);
}

// A progressive-alignment profile represented as distributions of ordered
// nucleotide pairs at two alignment columns.  Frequencies are divided by the
// total profile size, so sequences with a gap/ambiguous residue contribute
// zero rather than changing the scale between guide-tree levels.
class RibosumProfile
{
public:
  RibosumProfile(const ALN& alignment, const std::vector<Fasta>& sequences);

  float pair_score(unsigned i, unsigned j,
                   const RibosumProfile& other,
                   unsigned k, unsigned l) const;

private:
  using PairDistribution = std::array<float, 16>;

  PairDistribution pair_distribution(unsigned i, unsigned j) const;
  static std::uint64_t key(unsigned i, unsigned j);

  std::vector<std::vector<std::int8_t>> columns_;
  unsigned profile_size_;
  mutable std::unordered_map<std::uint64_t, PairDistribution> cache_;
};

#endif
