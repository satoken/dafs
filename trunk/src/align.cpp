// $Id:$

#include "config.h"
#include "align.h"

ProbCons::
ProbCons(float th)
  : Align(th), PROBCONS::Probcons()
{
}

float
ProbCons::
align(const std::string& seq1, const std::string& seq2, MP& mp)
{
  std::vector<float> posterior;
  ComputePosterior(seq1, seq2, posterior, threshold());
  const uint L1 = seq1.size();
  const uint L2 = seq2.size();

  std::vector<float> dp((L1+1)*(L2+1), 0.0);
  mp.resize(L1);
  for (uint i=0; i!=L1; ++i)
  {
    for (uint j=0; j!=L2; ++j)
    {
      const float& p=posterior[(L2+1)*(i+1)+(j+1)];
      if (p>threshold())
        mp[i].push_back(std::make_pair(j,p));
      float v = dp[(L2+1)*i+j] + p;
      v = std::max(v, dp[(L2+1)*(i+1)+j]);
      v = std::max(v, dp[(L2+1)*i+(j+1)]);
      dp[(L2+1)*(i+1)+(j+1)] = v;
    }
  }
  return dp.back();
}

CONTRAlign::
CONTRAlign(float th)
  : Align(th), CONTRALIGN::CONTRAlign<float>()
{
}

float
CONTRAlign::
align(const std::string& seq1, const std::string& seq2, MP& mp)
{
  std::vector<float> posterior;
  ComputePosterior(seq1, seq2, posterior, threshold());
  const uint L1 = seq1.size();
  const uint L2 = seq2.size();

  std::vector<float> dp((L1+1)*(L2+1), 0.0);
  mp.resize(L1);
  for (uint i=0; i!=L1; ++i)
  {
    for (uint j=0; j!=L2; ++j)
    {
      const float& p=posterior[(L2+1)*(i+1)+(j+1)];
      if (p>threshold())
        mp[i].push_back(std::make_pair(j,p));
      float v = dp[(L2+1)*i+j] + p;
      v = std::max(v, dp[(L2+1)*(i+1)+j]);
      v = std::max(v, dp[(L2+1)*i+(j+1)]);
      dp[(L2+1)*(i+1)+(j+1)] = v;
    }
  }
  return dp.back();
}

PartAlign::
PartAlign(float th)
  : Align(th), PARTALIGN::PartAlign<LogValue<float> >()
{
  float ribosum85_60[10] = {
    2.22, -1.86, -1.46, -1.39,
    /*  */ 1.16, -2.48, -1.05,
    /*         */ 1.03, -1.74,
    /*                */ 1.65
  };
  set_parameters(1.0, 1.0, -10.0, -5.0);
  set_scoring_matrix(ribosum85_60);
}

float
PartAlign::
align(const std::string& seq1, const std::string& seq2, MP& mp)
{
  const uint L1 = seq1.size();
  const uint L2 = seq2.size();
  VVF posterior;
  load_sequences(seq1, seq2);
  compute_forward();
  compute_backward();
  compute_posterior(posterior);

  VF dp((L1+1)*(L2+1), 0.0);
  mp.resize(L1);
  for (uint i=0; i!=L1; ++i)
  {
    for (uint j=0; j!=L2; ++j)
    {
      const float& p=posterior[i][j];
      if (p>threshold())
        mp[i].push_back(std::make_pair(j,p));
      float v = dp[(L2+1)*i+j] + p;
      v = std::max(v, dp[(L2+1)*(i+1)+j]);
      v = std::max(v, dp[(L2+1)*i+(j+1)]);
      dp[(L2+1)*(i+1)+(j+1)] = v;
    }
  }
  return dp.back();
}

float
PartAlign::
align(const std::string& seq1, const std::string& seq2,
      const BP& bp1, const BP& bp2, MP& mp)
{
  const uint L1 = seq1.size();
  const uint L2 = seq2.size();
  VVF posterior;
  load_sequences(seq1, seq2, bp1, bp2);
  compute_forward();
  compute_backward();
  compute_posterior(posterior);

  VF dp((L1+1)*(L2+1), 0.0);
  mp.resize(L1);
  for (uint i=0; i!=L1; ++i)
  {
    for (uint j=0; j!=L2; ++j)
    {
      const float& p=posterior[i][j];
      if (p>threshold())
        mp[i].push_back(std::make_pair(j,p));
      float v = dp[(L2+1)*i+j] + p;
      v = std::max(v, dp[(L2+1)*(i+1)+j]);
      v = std::max(v, dp[(L2+1)*i+(j+1)]);
      dp[(L2+1)*(i+1)+(j+1)] = v;
    }
  }
  return dp.back();
}
