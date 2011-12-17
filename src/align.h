// $Id:$

#ifndef __INC_ALIGN_H__
#define __INC_ALIGN_H__

#include <string>
#include <vector>

#include "typedefs.h"
#include "probconsRNA/probcons.h"
#include "contralign/contralign.h"
#if 0
#include "partalign/partalign.h"
#endif

namespace Align
{
  // base class for modeling probability distirubtion of alignments
  class Model
  {
  public:
    Model(float th) : th_(th) { }
    virtual ~Model() { }
    virtual void calculate(const std::string& seq1, const std::string& seq2, MP& mp) = 0;
    virtual void calculate(const std::string& seq1, const std::string& seq2,
                           const BP& bp1, const BP& bp2, MP& mp)
    {
      return calculate(seq1, seq2, mp);
    }
    float threshold() const { return th_; };

  private:
    float th_;
  };

  class Decoder
  {
  public:
    Decoder() { }
    virtual ~Decoder() { }
    virtual void initialize(const VVF& p) { }
    virtual float decode(const VVF& p, const VVF& q, VU& al) const = 0;
    virtual float decode(const VVF& p, VU& al) const = 0;
  };
}

class ProbCons : public Align::Model, PROBCONS::Probcons
{
public:
  ProbCons(float th);
  ~ProbCons() { }
  void calculate(const std::string& seq1, const std::string& seq2, MP& mp);
};

class CONTRAlign : public Align::Model, CONTRALIGN::CONTRAlign<float>
{
public:
  CONTRAlign(float th);
  ~CONTRAlign() { }
  void calculate(const std::string& seq1, const std::string& seq2, MP& mp);
};

#if 0
class PartAlign : public Align::Model, PARTALIGN::PartAlign<LogValue<float> >
{
public:
  PartAlign(float th, const std::string& arg);
  ~PartAlign() { }
  void calculate(const std::string& seq1, const std::string& seq2, MP& mp);
  void calculate(const std::string& seq1, const std::string& seq2,
                 const BP& bp1, const BP& bp2, MP& mp);
};
#endif

#endif  //  __INC_ALIGN_H__

// Local Variables:
// mode: C++
// End:
