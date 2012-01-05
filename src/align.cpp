// $Id:$

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "align.h"
#include <sys/errno.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iterator>
#include <cstdio>

void
Align::Model::
calculate(const std::vector<Fasta>& fa, std::vector<std::vector<MP> >& mp)
{
  const uint N=fa.size();
  mp.resize(N, std::vector<MP>(N));
  for (uint i=0; i!=N; ++i)
  {
    mp[i][i].resize(fa[i].size());
    for (uint x=0; x!=fa[i].size(); ++x)
      mp[i][i][x].push_back(std::make_pair(x, 1.0f));
    for (uint j=i+1; j!=N; ++j)
    {
      this->calculate(fa[i].seq(), fa[j].seq(), mp[i][j]);
      //transpose_mp(mp_[i][j], mp_[j][i], fa_[i].size(), fa_[j].size());
    }
  }
}

ProbCons::
ProbCons(float th)
  : Align::Model(th), PROBCONS::Probcons()
{
}

void
ProbCons::
calculate(const std::string& seq1, const std::string& seq2, MP& mp)
{
  std::vector<float> posterior;
  ComputePosterior(seq1, seq2, posterior, threshold());
  const uint L1 = seq1.size();
  const uint L2 = seq2.size();

  mp.resize(L1);
  for (uint i=0; i!=L1; ++i)
  {
    for (uint j=0; j!=L2; ++j)
    {
      const float& p=posterior[(L2+1)*(i+1)+(j+1)];
      if (p>threshold())
        mp[i].push_back(std::make_pair(j,p));
    }
  }
}

CONTRAlign::
CONTRAlign(float th)
  : Align::Model(th), CONTRALIGN::CONTRAlign<float>()
{
}

void
CONTRAlign::
calculate(const std::string& seq1, const std::string& seq2, MP& mp)
{
  std::vector<float> posterior;
  ComputePosterior(seq1, seq2, posterior, threshold());
  const uint L1 = seq1.size();
  const uint L2 = seq2.size();

  mp.resize(L1);
  for (uint i=0; i!=L1; ++i)
  {
    for (uint j=0; j!=L2; ++j)
    {
      const float& p=posterior[(L2+1)*(i+1)+(j+1)];
      if (p>threshold())
        mp[i].push_back(std::make_pair(j,p));
    }
  }
}

#if 0
PartAlign::
PartAlign(float th, const std::string& arg)
  : Align::Model(th), PARTALIGN::PartAlign<LogValue<float> >()
{
  float sm[10] = {                // ribosum85_60
    2.22, -1.86, -1.46, -1.39,
    /*  */ 1.16, -2.48, -1.05,
    /*         */ 1.03, -1.74,
    /*                */ 1.65
  };
  float alpha=1.0, beta=1.0, gap=-10.0, ext=-5.0;
  if (!arg.empty())
  {
    bool suc=true;
    const char* p=arg.c_str();
    suc &= sscanf(p, "%f", &alpha)==1; while (*p && *p!=',') ++p;
    suc &= sscanf(++p, "%f", &beta)==1; while (*p && *p!=',') ++p;
    suc &= sscanf(++p, "%f", &gap)==1; while (*p && *p!=',') ++p;
    suc &= sscanf(++p, "%f", &ext)==1; while (*p && *p!=',') ++p;
    for (uint i=0; i!=10; ++i)
    {
      suc &= sscanf(++p, "%f", &sm[i])==1; while (*p && *p!=',') ++p;
    }
    //std::cout << suc << std::endl;
  }
  set_parameters(alpha, beta, gap, ext);
  set_scoring_matrix(sm);
  //std::copy(sm, sm+10, std::ostream_iterator<float>(std::cout, " "));
  //std::cout << std::endl;
}

void
PartAlign::
calculate(const std::string& seq1, const std::string& seq2, MP& mp)
{
  const uint L1 = seq1.size();
  const uint L2 = seq2.size();
  VVF posterior;
  load_sequences(seq1, seq2);
  compute_forward();
  compute_backward();
  compute_posterior(posterior);

  mp.resize(L1);
  for (uint i=0; i!=L1; ++i)
  {
    for (uint j=0; j!=L2; ++j)
    {
      const float& p=posterior[i][j];
      if (p>threshold())
        mp[i].push_back(std::make_pair(j,p));
    }
  }
}

void
PartAlign::
calculate(const std::string& seq1, const std::string& seq2,
          const BP& bp1, const BP& bp2, MP& mp)
{
  const uint L1 = seq1.size();
  const uint L2 = seq2.size();
  VVF posterior;
  load_sequences(seq1, seq2, bp1, bp2);
  compute_forward();
  compute_backward();
  compute_posterior(posterior);

  mp.resize(L1);
  for (uint i=0; i!=L1; ++i)
  {
    for (uint j=0; j!=L2; ++j)
    {
      const float& p=posterior[i][j];
      if (p>threshold())
        mp[i].push_back(std::make_pair(j,p));
    }
  }
}
#endif

AUXAlign::
AUXAlign(const std::string& file, float th)
  : Align::Model(th),
    file_(file)
{
}

void
AUXAlign::
calculate(const std::string& seq1, const std::string& seq2, MP& mp)
{
  throw "not supported";
}

static
void
load_mp(std::istream& is, std::vector<std::vector<MP> >& mp)
{
  std::string s, t;
  uint x, y, i, k;
  float p;
  while (std::getline(is, s))
  {
    std::istringstream ss(s);
    if (s[0]=='>')
    {
      ss >> t >> x >> y;
      assert(x<y && x<mp.size() && y<mp[x].size());
    }
    else
    {
      ss >> i;
      if (i>=mp[x][y].size()) mp[x][y].resize(i+1);
      while (ss >> t)
      {
        if (sscanf(t.c_str(), "%d:%f", &k, &p)==2)
        {
          mp[x][y][i].push_back(std::make_pair(k, p));
        }
      }
    }
  }
}

void
AUXAlign::
calculate(const std::vector<Fasta>& fa, std::vector<std::vector<MP> >& mp)
{
  const uint N=fa.size();
  mp.resize(N, std::vector<MP>(N));
  std::ifstream is(file_.c_str());
  if (is.is_open())
    load_mp(is, mp);
  else
    throw strerror(errno);
}
