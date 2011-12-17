// $Id:$

#ifndef __INC_NEEDLEMAN_WUNSCH_H__
#define __INC_NEEDLEMAN_WUNSCH_H__


#include "align.h"

class NeedlemanWunsch : public Align::Decoder
{
public:
  NeedlemanWunsch(float th) : Align::Decoder(), th_(th) { }
  void initialize(const VVF& p) { }
  float decode(const VVF& p, const VVF& q, VU& al) const;
  float decode(const VVF& p, VU& al) const;

private:
  float th_;
};

class SparseNeedlemanWunsch : public Align::Decoder
{
public:
  SparseNeedlemanWunsch(float th) : Align::Decoder(), th_(th), env_() { }
  void initialize(const VVF& p);
  float decode(const VVF& p, const VVF& q, VU& al) const;
  float decode(const VVF& p, VU& al) const;

private:
  float th_;
  std::vector< std::pair<uint,uint> > env_;
};

#endif  //  __INC_NEEDLEMAN_WUNSCH_H__

// Local Variables:
// mode: C++
// End:
