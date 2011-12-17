// $Id:$

#ifndef __INC_FOLD_H__
#define __INC_FOLD_H__

#include <string>
#include <vector>

#include "typedefs.h"
#include "contrafold/contrafold.h"

namespace Fold 
{
  // base class for modeling probability distribution of secondary structures
  class Model
  {
  public:
    Model(float th) : th_(th) { }
    virtual ~Model() { }
    virtual void calculate(const std::string& seq, BP& bp) = 0;
    float threshold() const { return th_; }

  private:
    float th_;
  };

  class Decoder
  {
  public:
    Decoder() { }
    virtual ~Decoder() { }
    virtual float decode(const VVF& p, const VVF& q, VU& ss) const = 0;
    virtual float decode(const VVF& p, VU& ss) const = 0;
  };
}

class RNAfold : public Fold::Model
{
public:
  RNAfold(bool bl, const char* param, float th);
  void calculate(const std::string& seq, BP& bp);
};

class CONTRAfold : public Fold::Model, CONTRAFOLD::CONTRAfold<float>
{
public:
  CONTRAfold(float th);
  void calculate(const std::string& seq, BP& bp);
};

#endif  //  __INC_FOLD_H__

// Local Variables:
// mode: C++
// End:
