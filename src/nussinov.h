// $Id:$

#ifndef __INC_NUSSINOV_H__
#define __INC_NUSSINOV_H__

#include "fold.h"

class Nussinov : public Fold::Decoder
{
public:
  Nussinov(float w, float th) : Fold::Decoder(), w_(w), th_(th) { }
  float decode(const VVF& p, const VVF& q, VU& ss) const;
  float decode(const VVF& p, VU& ss) const;

private:
  float w_;
  float th_;
};

class SparseNussinov : public Fold::Decoder
{
public:
  SparseNussinov(float w, float th) : Fold::Decoder(), w_(w), th_(th) { }
  float decode(const VVF& p, const VVF& q, VU& ss) const;
  float decode(const VVF& p, VU& ss) const;

private:
  float w_;
  float th_;
};

#endif  //  __INC_NUSSINOV_H__

// Local Variables:
// mode: C++
// End:
