#include "nussinov.h"
#include "dd_recovery.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <random>
#include <stdexcept>
#include <type_traits>
#ifdef _OPENMP
#include <omp.h>
#endif

const uint Fold::Decoder::n_support_brackets = 30;
const char* Fold::Decoder::left_brackets = "([{<ABCDEFGHIJKLMNOPQRSTUVWXYZ";
const char* Fold::Decoder::right_brackets = ")]}>abcdefghijklmnopqrstuvwxyz";

namespace {
bool sanitizer_smoke = false;

void require(bool value, const char* message)
{
  if (!value) throw std::runtime_error(message);
}
void same(const LinearNussinovResult& a, const LinearNussinovResult& b)
{
  require(a.score == b.score && a.upper_bound == b.upper_bound &&
          a.pruned_upper_bound == b.pruned_upper_bound &&
          a.additive_upper_bound == b.additive_upper_bound &&
          a.pruned_states == b.pruned_states,
          "serial/parallel beam score or certificate differs");
}
SparseFloatMatrix sparse(const VVF& dense)
{
  SparseFloatMatrix result;
  result.assign(dense.size(), dense.size());
  for (uint i=0; i<dense.size(); ++i)
    for (uint j=i+2; j<dense.size(); ++j)
      if (dense[i][j] != 0) result.set(i,j,dense[i][j]);
  return result;
}
template <typename Serial, typename Parallel>
void adapters(Serial& serial, Parallel& parallel, const VVF& p, const VVF& q,
              bool allow)
{
  const auto ps=sparse(p), qs=sparse(q);
  VU a,b;
  const auto check=[&](float x, float y) {
    require(x == y && a == b, "serial/parallel adapter result differs");
    require(DDRecovery::structure(b,allow), "invalid parallel structure");
  };
  check(serial.decode(1.7f,p,q,a),parallel.decode(1.7f,p,q,b));
  Fold::Decoder& sa=serial;
  Fold::Decoder& pa=parallel;
  if constexpr (!std::is_same_v<Serial,Nussinov>)
    check(sa.decode(1.7f,ps,qs,a),pa.decode(1.7f,ps,qs,b));
  check(sa.decode(1.7f,ps,q,a),pa.decode(1.7f,ps,q,b));
  if constexpr (!std::is_same_v<Serial,Nussinov>)
    check(sa.decode(1.7f,p,qs,a),pa.decode(1.7f,p,qs,b));
  std::string astr,bstr;
  check(sa.decode(p,a,astr),pa.decode(p,b,bstr));
  require(astr==bstr,"parallel final brackets differ");
  check(sa.decode(ps,a,astr),pa.decode(ps,b,bstr));
  check(sa.decode(ps,qs,a,astr),pa.decode(ps,qs,b,bstr));
  require(astr==bstr,"parallel bonus brackets differ");
}

void validation()
{
  const auto rejects=[](const auto& make) {
    bool rejected=false;
    try { make(); } catch (const std::invalid_argument&) { rejected=true; }
    require(rejected,"invalid thread count accepted");
  };
  rejects([] { Nussinov d(.2f,true,-1); });
  rejects([] { SparseNussinov d(.2f,true,-1); });
  rejects([] { LinearNussinov d(.2f,100,true,-1); });
#ifndef _OPENMP
  rejects([] { Nussinov d(.2f,true,2); });
  rejects([] { SparseNussinov d(.2f,true,2); });
  rejects([] { LinearNussinov d(.2f,100,true,2); });
#endif
}

void exact_and_adapters()
{
  std::mt19937 rng(9137);
#ifdef _OPENMP
  const std::vector<int> requests = sanitizer_smoke
      ? std::vector<int>{4} : std::vector<int>{0,2,4};
#else
  const std::vector<int> requests{0,1};
#endif
  for (bool allow : {true,false}) for (int threads : requests) {
    Nussinov dense(.2f,allow),dense_parallel(.2f,allow,threads);
    SparseNussinov exact(.2f,allow),exact_parallel(.2f,allow,threads);
    LinearNussinov beam(.2f,13,allow),beam_parallel(.2f,13,allow,threads);
    for (uint n : {64,7,0,33}) {
      VVF p(n,VF(n)),q=p;
      for (uint i=0;i<n;++i) for (uint j=i+3;j<n;++j) {
        p[i][j]=.25f*(1+rng()%4);
        q[i][j]=.25f*int(rng()%9)-1;
      }
      adapters(dense,dense_parallel,p,q,allow);
      adapters(exact,exact_parallel,p,q,allow);
      adapters(beam,beam_parallel,p,q,allow);
      require(beam_parallel.used_workers()==1 && beam_parallel.parallel_steps()==0,
              "small beam should fall back to serial");
      if (n>1) {
#ifdef _OPENMP
        require(dense_parallel.used_workers()>1 && dense_parallel.parallel_steps()==n-1,
                "dense exact parallel path was not exercised");
        // Last call used sparse bonus support; this dense call forces wavefront.
        VU a,b;
        require(exact.decode(1,p,q,a)==exact_parallel.decode(1,p,q,b) && a==b,
                "parallel exact wavefront differs");
        require(exact_parallel.used_workers()>1 && exact_parallel.parallel_steps()==n-1,
                "IPknot wavefront parallel path was not exercised");
#endif
      }
    }
  }
}

void beam_parallel()
{
  std::mt19937 rng(38433200);
#ifdef _OPENMP
  const std::vector<int> requests = sanitizer_smoke
      ? std::vector<int>{4} : std::vector<int>{0,2,4};
#else
  const std::vector<int> requests{0};
#endif
  const uint n=sanitizer_smoke ? 320 : 384;
  for (bool allow : {true,false}) {
    LinearNussinov serial(.2f,256,allow);
    std::vector<std::unique_ptr<LinearNussinov>> decoders;
    for (int threads : requests)
      decoders.push_back(std::make_unique<LinearNussinov>(.2f,256,allow,threads));
    for (int trial=sanitizer_smoke ? 2 : 0;trial<3;++trial) {
      std::cout << "Beam case: allow=" << allow << ", trial=" << trial << std::endl;
      VVF p(n,VF(n)),q=p;
      std::vector<std::vector<uint>> support(n);
      for (uint j=0;j<n;++j) {
        for (uint i=0;i+2<=j;++i) {
          p[i][j]=1;
          q[i][j]=trial==0 ? 0 : .125f*int(rng()%9)-.5f;
          if (trial==2 && rng()%5==0) q[i][j]=1.125f; // signed stack partners
          support[j].push_back(i);
        }
        if (j>5) support[j].push_back(1); // duplicate endpoints
        support[j].push_back(n+7); // invalid cached entries
        support[j].push_back(j);
        std::shuffle(support[j].begin(),support[j].end(),rng);
      }
      const auto ps=sparse(p),qs=sparse(q);
      VU expected_ss,actual_ss;
      const auto expected=serial.decode_certified(1.0f,ps,qs,support,expected_ss);
      require(DDRecovery::structure(expected_ss,allow),"invalid serial structure");
      for (const auto& decoder : decoders) {
        auto actual=decoder->decode_certified(1.0f,ps,qs,support,actual_ss);
        same(expected,actual);
        require(expected_ss==actual_ss,"parallel cached traceback differs");
#ifdef _OPENMP
        require(decoder->used_workers()>1 && decoder->parallel_steps()==n,
                "beam parallel path was not exercised");
#endif
        actual=decoder->decode_certified(1.0f,ps,q,support,actual_ss);
        same(expected,actual);
        require(expected_ss==actual_ss,"parallel dense-multiplier traceback differs");
      }
      // Standard overloads use span >=3, whereas cached supports also allow 2.
      if (!sanitizer_smoke) adapters(serial,*decoders.back(),p,q,allow);
      else {
        // All adapter formats were covered above on smaller signed inputs;
        // exercise the uncertified parallel kernel once on this wide chart.
        const float expected_score=serial.decode(1,p,q,expected_ss);
        const float actual_score=decoders.back()->decode(1,p,q,actual_ss);
        require(expected_score==actual_score && expected_ss==actual_ss,
                "uncertified parallel traceback differs");
      }
    }
  }
}

void pruned_parallel_certificate()
{
  const uint n=sanitizer_smoke ? 600 : 768;
  const uint beam=sanitizer_smoke ? 64 : 128;
  SparseFloatMatrix p,q;
  p.assign(n,n);q.assign(n,n);
  for (uint i=0;i<n;++i)
    for (uint j=i+3;j<n;++j) p.set(i,j,1);
  for (bool allow : {true,false}) {
    LinearNussinov serial(.2f,beam,allow),parallel(.2f,beam,allow,0);
    VU a,b;
    const auto expected=serial.decode_certified(1,p,q,a);
    const auto actual=parallel.decode_certified(1,p,q,b);
    same(expected,actual);
    require(a==b && actual.pruned_states>0,"parallel pruning/certificate not exercised");
#ifdef _OPENMP
    require(parallel.used_workers()==4,"four-worker beam path not exercised");
#endif
    require(DDRecovery::structure(b,allow),"invalid pruned parallel traceback");
  }
}

void kernel_selection()
{
  // Reuse one object across IPknot -> DAFS -> changed IPknot -> empty ->
  // IPknot. This exercises both caches and same-length graph rebuilding.
  const uint n=sanitizer_smoke ? 256 : 320;
  const uint beam=sanitizer_smoke ? 8 : 64;
  for (bool allow : {true,false}) {
    LinearNussinov reused(.2f,beam,allow,0);
    for (int phase=0;phase<5;++phase) {
      const uint length=phase==3 ? 0 : n;
      SparseFloatMatrix p,q;
      p.assign(length,length);q.assign(length,length);
      if (length && phase==1) {
        p.set(0,n-1,.9f);p.set(1,n-2,.1f);
      } else {
        for (uint j=3;j<length;++j) for (uint i=0;i+3<=j;++i) {
          if (phase==2 && (i+j)%5==0) continue;
          p.set(i,j,.125f*(1+(i+3*j)%7));
          q.set(i,j,.125f*int((i+j)%5)-.25f);
        }
      }
      LinearNussinov fresh(.2f,beam,allow,0);
      VU a,b;
      const auto expected=fresh.decode_certified(1,p,q,a);
      const auto actual=reused.decode_certified(1,p,q,b);
      same(expected,actual);
      require(a==b,"switching beam kernels changed traceback");
      require(reused.used_ipknot_kernel()==(phase==0 || phase==2 || phase==4),
              "automatic beam kernel selection differs");
      require(DDRecovery::structure(b,allow),"invalid switched-kernel structure");
    }
  }
}

void sparse_long()
{
  const uint n=sanitizer_smoke ? 2000 : 20000;
  SparseFloatMatrix p,q;
  p.assign(n,n);q.assign(n,n);
  p.set(0,n-1,.95f);p.set(1,n-2,.1f);
  for (bool allow : {true,false}) {
    SparseNussinov serial(.2f,allow),parallel(.2f,allow,0);
    VU a,b;
    require(serial.decode(1,p,q,a)==parallel.decode(1,p,q,b) && a==b,
            "long sparse exact differs");
    require(parallel.used_workers()==1 && parallel.parallel_steps()==0,
            "long sparse exact allocated a parallel dense chart");
    LinearNussinov beam(.2f,4,allow),beam_parallel(.2f,4,allow,0);
    same(beam.decode_certified(1,p,q,a),beam_parallel.decode_certified(1,p,q,b));
    require(a==b && beam_parallel.used_workers()==1,"long sparse beam differs");
  }
}
} // namespace
int main(int argc, char** argv)
{
  try {
    if (argc==2 && std::string(argv[1])=="--sanitizer-smoke") sanitizer_smoke=true;
    else require(argc==1,"unknown test argument");
#ifdef _OPENMP
    omp_set_dynamic(0);
    omp_set_num_threads(4); // bounds the explicit auto-runtime test
#endif
    validation();
    std::cout << "Checking exact and adapter parity" << std::endl;
    exact_and_adapters();
    std::cout << "Checking parallel beam parity" << std::endl;
    beam_parallel();
    std::cout << "Checking parallel pruning certificates" << std::endl;
    pruned_parallel_certificate();
    std::cout << "Checking automatic beam kernel selection" << std::endl;
    kernel_selection();
    std::cout << "Checking long sparse fallback" << std::endl;
    sparse_long();
    std::cout << "Nussinov serial/parallel parity passed\n";
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';return 1;
  }
}
