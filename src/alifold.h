// $Id:$

#ifndef __INC_ALIFOLD_H__
#define __INC_ALIFOLD_H__

#include <vector>
#include <string>
#include "fa.h"
#include "typedefs.h"

class Alifold
{
public:
  Alifold(float th) : th_(th) { }
  void fold(const ALN& aln, const std::vector<Fasta>& fa, BP& bp) const;
  void fold(const ALN& aln, const std::vector<Fasta>& fa,
            const std::string& str, BP& bp) const;

private:
  static char** alloc_aln(const ALN& aln, const std::vector<Fasta>& fa);
  static void free_aln(char** seqs);
  
private:
  float th_;
};

#endif  //  __INC_ALIFOLD_H__

// Local Variables:
// mode: C++
// End:
