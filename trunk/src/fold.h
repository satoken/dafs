// $Id:$

#ifndef __INC_FOLD_H__
#define __INC_FOLD_H__

#include <string>
#include <vector>

#include "typedefs.h"
#include "contrafold/contrafold.h"

// base class
class Fold
{
public:
  Fold(float th) : th_(th) { }
  virtual void fold(const std::string& seq, BP& bp) = 0;
  float threshold() const { return th_; }

private:
  float th_;
};

class RNAfold : public Fold
{
public:
  RNAfold(bool bl, const char* param, float th);
  void fold(const std::string& seq, BP& bp);
};

class CONTRAfold : public Fold, CONTRAFOLD::CONTRAfold<float>
{
public:
  CONTRAfold(float th);
  void fold(const std::string& seq, BP& bp);
};

#endif  //  __INC_FOLD_H__
