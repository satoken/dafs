// $Id:$

#ifndef __INC_ALIGN_H__
#define __INC_ALIGN_H__

#include <string>
#include <vector>

#include "typedefs.h"
#include "probconsRNA/probcons.h"
#include "contralign/contralign.h"
#include "partalign/partalign.h"

// base class
class Align
{
public:
  Align(float th) : th_(th) { }
  virtual float align(const std::string& seq1, const std::string& seq2, MP& mp) = 0;
  virtual float align(const std::string& seq1, const std::string& seq2,
                      const BP& bp1, const BP& bp2, MP& mp)
  {
    throw "not implemented";
  }
  float threshold() const { return th_; };

private:
  float th_;
};

class ProbCons : public Align, PROBCONS::Probcons
{
public:
  ProbCons(float th);
  ~ProbCons() { }
  float align(const std::string& seq1, const std::string& seq2, MP& mp);
};

class CONTRAlign : public Align, CONTRALIGN::CONTRAlign<float>
{
public:
  CONTRAlign(float th);
  ~CONTRAlign() { }
  float align(const std::string& seq1, const std::string& seq2, MP& mp);
};

class PartAlign : public Align, PARTALIGN::PartAlign<LogValue<float> >
{
public:
  PartAlign(float th);
  ~PartAlign() { }
  float align(const std::string& seq1, const std::string& seq2, MP& mp);
  float align(const std::string& seq1, const std::string& seq2,
              const BP& bp1, const BP& bp2, MP& mp);
};

#endif  //  __INC_ALIGN_H__
