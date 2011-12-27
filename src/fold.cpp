// $Id:$

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "fold.h"
#include <algorithm>

namespace Vienna {
extern "C" {
#include <ViennaRNA/fold.h>
#include <ViennaRNA/fold_vars.h>
#include <ViennaRNA/part_func.h>
#include <ViennaRNA/alifold.h>
#include <ViennaRNA/aln_util.h>
#include <ViennaRNA/utils.h>
#include <ViennaRNA/PS_dot.h>
  extern void read_parameter_file(const char fname[]);
};
};

extern "C" {
#include "boltzmann_param.h"
};

// some constants
const uint Fold::Decoder::n_support_brackets=4+26;
const char* Fold::Decoder::left_brackets ="([{<ABCDEFGHIJKLMNOPQRSTUVWXYZ";
const char* Fold::Decoder::right_brackets=")]}>abcdefghijklmnopqrstuvwxyz";


RNAfold::
RNAfold(bool bl, const char* param, float th)
  : Fold::Model(th)
{
  if (bl) copy_boltzmann_parameters();
  if (param) Vienna::read_parameter_file(param);
}

void
RNAfold::
calculate(const std::string& seq, BP& bp)
{
  uint L=seq.size();
  bp.resize(L);
#if 0
  std::string str(seq.size()+1, '.');
  float min_en = Vienna::fold(const_cast<char*>(seq.c_str()), &str[0]);
  float sfact = 1.07;
  float kT = (Vienna::temperature+273.15)*1.98717/1000.; /* in Kcal */
  Vienna::pf_scale = exp(-(sfact*min_en)/kT/seq.size());
#else
  Vienna::pf_scale = -1;
#endif
  Vienna::init_pf_fold(L);
  Vienna::pf_fold(const_cast<char*>(seq.c_str()), NULL);
  for (uint i=0; i!=L-1; ++i)
    for (uint j=i+1; j!=L; ++j)
    {
      const float& p = Vienna::pr[Vienna::iindx[i+1]-(j+1)];
      if (p>threshold())
        bp[i].push_back(std::make_pair(j, p));
    }
  Vienna::free_pf_arrays();
}

void
RNAfold::
calculate(const std::string& seq, const std::string& str, BP& bp)
{
  std::string p(str);
  std::replace(p.begin(), p.end(), '.', 'x');
  std::replace(p.begin(), p.end(), '?', '.');

  int bk = Vienna::fold_constrained;
  Vienna::fold_constrained = 1;

  uint L=seq.size();
  bp.resize(L);
#if 0
  std::string str(seq.size()+1, '.');
  float min_en = Vienna::fold(const_cast<char*>(seq.c_str()), &str[0]);
  float sfact = 1.07;
  float kT = (Vienna::temperature+273.15)*1.98717/1000.; /* in Kcal */
  Vienna::pf_scale = exp(-(sfact*min_en)/kT/seq.size());
#else
  Vienna::pf_scale = -1;
#endif
  Vienna::init_pf_fold(L);
  Vienna::pf_fold(const_cast<char*>(seq.c_str()), &p[0]);
  for (uint i=0; i!=L-1; ++i)
    for (uint j=i+1; j!=L; ++j)
    {
      const float& p = Vienna::pr[Vienna::iindx[i+1]-(j+1)];
      if (p>threshold())
        bp[i].push_back(std::make_pair(j, p));
    }
  Vienna::free_pf_arrays();
  Vienna::fold_constrained = bk;
}

CONTRAfold::
CONTRAfold(float th)
  : Fold::Model(th), CONTRAFOLD::CONTRAfold<float>()
{
}

void
CONTRAfold::
calculate(const std::string& seq, BP& bp)
{
  bp.resize(seq.size());
  std::vector<float> posterior;
  ComputePosterior(seq, posterior);
  for (uint i=0, k=0; i!=seq.size()+1; ++i)
  {
    for (uint j=i; j!=seq.size()+1; ++j, ++k)
    {
      if (i!=0 && posterior[k]>threshold())
        bp[i-1].push_back(std::make_pair(j-1, posterior[k]));
    }
  }
}

void
CONTRAfold::
calculate(const std::string& seq, const std::string& str, BP& bp)
{
  bp.resize(seq.size());
  std::vector<float> posterior;
  SetConstraint(str);
  ComputePosterior(seq, posterior);
  for (uint i=0, k=0; i!=seq.size()+1; ++i)
  {
    for (uint j=i; j!=seq.size()+1; ++j, ++k)
    {
      if (i!=0 && posterior[k]>threshold())
        bp[i-1].push_back(std::make_pair(j-1, posterior[k]));
    }
  }
}
