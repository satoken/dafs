#include "dd_certificate.h"
#include <random>
#include <iostream>

int main() {
  for(double v:{0.0,0.1,-0.1,1e20,-1e20,1e-30,-1e-30})
    if(static_cast<double>(DDCertificate::lower_float(v))>v) return 5;
  for(float w:{0.0f,0.3f,4.0f}) for(float p:{0.0f,0.2f,0.99f})
    if(DDCertificate::lower_weighted_score(w,p,0.1f)>
        static_cast<long double>(w)*(static_cast<long double>(p)-0.1f)) return 6;
  std::mt19937 rng(987);
  for (unsigned trial=0;trial<2000;++trial) {
    const unsigned n=1+rng()%97,m=1+rng()%29;
    std::vector<float> a(n*m,-1.0f);
    std::vector<std::pair<unsigned,unsigned>> support;
    for(unsigned i=0;i<n;++i) for(unsigned j=0;j<m;++j) {
      if(rng()%3==0) continue;
      a[i*m+j]=(static_cast<int>(rng()%33)-8)/8.0f;
      support.emplace_back(i,j);
      if(rng()%4==0) support.emplace_back(i,j);
    }
    std::shuffle(support.begin(),support.end(),rng);
    auto exact=[&](unsigned start,unsigned end) {
      std::vector<double> dp((end-start+1)*(m+1),0);
      for(unsigned i=1;i<=end-start;++i) for(unsigned j=1;j<=m;++j)
        dp[i*(m+1)+j]=std::max({dp[(i-1)*(m+1)+j],dp[i*(m+1)+j-1],
            dp[(i-1)*(m+1)+j-1]+a[(start+i-1)*m+j-1]});
      return dp.back();
    };
    for(unsigned width:{1u,2u,4u,8u,16u,32u,64u}) {
      double expected=0;
      for(unsigned i=0;i<n;i+=width) expected+=exact(i,std::min(n,i+width));
      const float got=DDCertificate::block_alignment(support,n,m,width,
          [&](unsigned i,unsigned j){return a[i*m+j];});
      if(got<expected || got<exact(0,n) || got>expected+1e-4) {
        std::cerr<<"block oracle failure "<<trial<<' '<<width<<' '<<expected<<' '<<got<<'\n';return 1;
      }
    }
    if(DDCertificate::unpruned_alignment(exact(0,n),std::min(n,m))<exact(0,n)) return 2;
  }
  try { DDCertificate::block_alignment({},0,0,65,[](unsigned,unsigned){return 0.0f;});return 3; }
  catch(const std::invalid_argument&) {}
  if(DDCertificate::block_alignment({},0,0,8,[](unsigned,unsigned){return 0.0f;})!=0) return 4;
  return 0;
}
