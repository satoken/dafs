#pragma once

#if RNA
namespace CONTRALIGN_RNA
#else
namespace CONTRALIGN
#endif
{
  template <class T>
  class CONTRAlign {
  private:
    struct Impl;

  public:
    CONTRAlign ();
    ~CONTRAlign ();

    const T* ComputePosterior (const std::string& seq1, const std::string& seq2, float th=0.0);
    const T* ComputePosterior (const std::string& seq1, const std::string& seq2, 
			       std::vector<T>& p, float th=0.0);
  private:
    Impl* impl_;
  };
}