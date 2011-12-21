// $Id:$

#ifndef __INC_IPKNOT_H__
#define __INC_IPKNOT_H__

#include "fold.h"
#include "ip.h"

class IPknot : public Fold::Decoder
{
public:
  IPknot(float w, const VF& th, int n_th = 1);
  float decode(const VVF& p, const VVF& q, VU& ss);
  float decode(const VVF& p, VU& ss, std::string& str);
  void make_parenthsis(const VU& ss, std::string& str) const;

private:
  void make_objective(IP& ip, const VVF& p, const VVF& q);
  void make_objective(IP& ip, const VVF& p);
  void make_constraints(IP& ip);
  float solve(IP& ip, VU& ss);
  
private:
  float weight_;
  VF th_;
  VF alpha_;
  bool levelwise_;
  bool stacking_constraints_;
  std::vector<VVI> v_;
  std::vector<VVI> w_;
  uint n_th_;
  VU plevel_;
};

#endif  //  __INC_IPKNOT_H__
