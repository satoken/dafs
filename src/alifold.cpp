// $Id:$

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include "alifold.h"

namespace Vienna {
extern "C" {
#include <ViennaRNA/fold.h>
#include <ViennaRNA/fold_vars.h>
#include <ViennaRNA/part_func.h>
#include <ViennaRNA/alifold.h>
#include <ViennaRNA/aln_util.h>
#include <ViennaRNA/utils.h>
#include <ViennaRNA/PS_dot.h>
};
};

#define FOREACH(itr, i, v) for (itr i=(v).begin(); i!=(v).end(); ++i)

#ifdef HAVE_VIENNA18
  typedef Vienna::plist pair_info;
#else
  typedef Vienna::pair_info pair_info;
#endif

void
Alifold::
fold(const ALN& aln, const std::vector<Fasta>& fa, BP& bp) const
{
  const uint L=aln.front().second.size();
  //const uint N=aln.size();

  std::string res(L+1, ' ');
  char** seqs = alloc_aln(aln, fa);
    
  // scaling parameters to avoid overflow
  double min_en = Vienna::alifold(seqs, &res[0]);
  double kT = (Vienna::temperature+273.15)*1.98717/1000.; /* in Kcal */
  Vienna::pf_scale = exp(-(1.07*min_en)/kT/L);
  Vienna::free_alifold_arrays();

  pair_info* pi;
  Vienna::alipf_fold(seqs, NULL, &pi);

  bp.resize(L);
  for (uint k=0; pi[k].i!=0; ++k)
    if (pi[k].p>th_)
      bp[pi[k].i-1].push_back(std::make_pair(pi[k].j-1, pi[k].p));

  free(pi);
  Vienna::free_alipf_arrays();
  free_aln(seqs);
}

void
Alifold::
fold(const ALN& aln, const std::vector<Fasta>& fa, const std::string& str, BP& bp) const
  {
    const uint L=aln.front().second.size();
    std::string p(str);
    std::replace(p.begin(), p.end(), '.', 'x');
    std::replace(p.begin(), p.end(), '?', '.');

    int bk = Vienna::fold_constrained;
    Vienna::fold_constrained = 1;

    char** seqs = alloc_aln(aln, fa);

    // scaling parameters to avoid overflow
    std::string res(p);
    double min_en = Vienna::alifold(seqs, &res[0]);
    double kT = (Vienna::temperature+273.15)*1.98717/1000.; /* in Kcal */
    Vienna::pf_scale = exp(-(1.07*min_en)/kT/L);
    Vienna::free_alifold_arrays();

    pair_info* pi;
    Vienna::alipf_fold(seqs, &p[0], &pi);

    bp.resize(L);
    for (uint k=0; pi[k].i!=0; ++k)
      if (pi[k].p>th_)
        bp[pi[k].i-1].push_back(std::make_pair(pi[k].j-1, pi[k].p));

    free(pi);
    Vienna::free_alipf_arrays();
    free_aln(seqs);
    Vienna::fold_constrained = bk;
  }

//static
char**
Alifold::
alloc_aln(const ALN& aln, const std::vector<Fasta>& fa)
{
  const uint L=aln.front().second.size();
  const uint N=aln.size();

  char** seqs = new char*[N+1];
  seqs[N] = NULL;
  char** s = seqs;
  FOREACH (ALN::const_iterator, it, aln)
  {
    *s = new char[L+1];
    for (uint i=0, j=0; i!=L; ++i)
      (*s)[i] = it->second[i] ? fa[it->first].seq()[j++] : '-';
    (*s)[L] = 0;
    ++s;
  }
  return seqs;
}

//static
void
Alifold::free_aln(char** seqs)
{
  for (char** s=seqs; *s!=NULL; ++s)
    delete[] *s;
  delete[] seqs;
}
