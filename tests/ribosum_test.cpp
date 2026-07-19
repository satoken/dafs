#include "ribosum.h"

#include <cmath>
#include <iostream>

namespace
{
bool close(float observed, float expected)
{
  if (std::abs(observed - expected) <= 1e-5f)
    return true;
  std::cerr << "observed " << observed << ", expected " << expected << '\n';
  return false;
}
}

int main()
{
  // Ordered-pair indices: AU=3, CG=6, GC=9, GU=11, UA=12, UG=14.
  if (!close(Ribosum85_60::pair_score(3, 3), 4.492700f) ||
      !close(Ribosum85_60::pair_score(3, 6), 1.673203f) ||
      !close(Ribosum85_60::pair_score(6, 3), 1.673203f) ||
      !close(Ribosum85_60::pair_score(9, 9), 5.616325f) ||
      !close(Ribosum85_60::pair_score(11, 14), -2.088271f))
    return 1;

  std::vector<Fasta> sequences = {
      Fasta("au", "AU"), Fasta("cg", "CG"),
      Fasta("gc", "GC"), Fasta("gap", "AU")};
  const ALN au = {{0, {true, true}}};
  const ALN cg = {{1, {true, true}}};
  const ALN mixed = {{2, {true, true, false}}, {3, {true, false, true}}};
  const RibosumProfile au_profile(au, sequences);
  const RibosumProfile cg_profile(cg, sequences);
  const RibosumProfile mixed_profile(mixed, sequences);

  if (!close(au_profile.pair_score(0, 1, cg_profile, 0, 1), 1.673203f))
    return 1;
  // The second member has a gap at column 1, so only half of the profile
  // contributes; GC versus AU is 2.704820 in RIBOSUM85-60.
  if (!close(mixed_profile.pair_score(0, 1, au_profile, 0, 1),
             0.5f * 2.704820f))
    return 1;
  return 0;
}
