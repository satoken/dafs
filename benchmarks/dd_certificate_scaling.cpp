#include "dd_certificate.h"
#include <chrono>
#include <iostream>
int main() {
  std::cout<<"length,edges,width,median_seconds\n";
  volatile float sink=0;
  for (unsigned n:{1000,2000,4000,8000,16000}) {
    std::vector<std::pair<unsigned,unsigned>> support;
    for(unsigned i=0;i<n;++i) for(unsigned d=0;d<8;++d)
      support.emplace_back(i,(i+17*d)%n);
    for(unsigned w:{8,32}) {
      std::vector<double> times;
      for(unsigned rep=0;rep<11;++rep) {
        const auto start=std::chrono::steady_clock::now();
        sink=DDCertificate::block_alignment(support,n,n,w,[](unsigned i,unsigned j){return 0.1f+(i+j)%13*0.05f;});
        times.push_back(std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count());
      }
      std::sort(times.begin(),times.end());
      std::cout<<n<<','<<support.size()<<','<<w<<','<<times[5]<<'\n';
    }
  }
  return sink<0;
}
