// Compare the complete sparse model with the bundled nonlinear CONTRAlign.
#include <iostream>
#include <random>
#include <stdexcept>
#include "linearcontralign5.h"
#include "contralign/contralign.h"

void require(bool value,const char* message) {if(!value)throw std::runtime_error(message);}

int main() {
  LinearCONTRAlign5::Model model;
  CONTRALIGN::CONTRAlign<float> reference;
  std::mt19937 rng(42);
  const std::string alphabet="ACGUNacgu";
  std::vector<std::pair<std::string,std::string>> examples={
    {"A","A"},{"A","C"},{"A","AAAAAAAAAAAA"},{"AAAAAAAAAAAA","A"},
    {"ACGUACGU","ACGUUUUUUUUUACGU"},{"NNNAU","AUNNN"},
    {"GCGCAAAAGCGC","GCGCGCGC"},{"GGGGCCCC","GGAGCCCC"}};
  for(int case_id=0;case_id<300;++case_id) {
    const int n=1+rng()%45,m=1+rng()%45;
    std::string x,y;
    for(int i=0;i<n;++i)x+=alphabet[rng()%alphabet.size()];
    for(int j=0;j<m;++j)y+=alphabet[rng()%alphabet.size()];
    examples.emplace_back(x,y);
  }
  double largest_error=0,total_error=0;size_t values=0;
  double second_gap_visits=0;
  for(const auto& example:examples) {
    const auto& x=example.first;const auto& y=example.second;
    std::vector<float> dense;
    reference.ComputePosterior(x,y,dense,0);
    MP exact;LinearCONTRAlign5::Diagnostics full;
    model.calculate(x,y,0,0,exact,&full);
    require(full.pruned_cells==0,"Unpruned diagnostic unexpectedly pruned");
    require(full.retained_cells==(x.size()+1)*(y.size()+1),"Missing unpruned lattice cells");
    std::vector<std::vector<float>> actual(x.size(),std::vector<float>(y.size()));
    for(size_t i=0;i<exact.size();++i)for(const auto& p:exact[i])actual[i][p.first]=p.second;
    for(size_t i=0;i<x.size();++i)for(size_t j=0;j<y.size();++j) {
      const double error=std::abs(actual[i][j]-dense[(i+1)*(y.size()+1)+j+1]);
      largest_error=std::max(largest_error,error);total_error+=error;++values;
    }
    second_gap_visits+=full.expected_state_visits[LinearCONTRAlign5::INS2_X]
                      +full.expected_state_visits[LinearCONTRAlign5::INS2_Y];
    for(unsigned beam:{1u,2u,7u,50u}) {
      MP sparse;LinearCONTRAlign5::Diagnostics stats;
      model.calculate(x,y,beam,0,sparse,&stats);
      require(stats.forward_log_z<=full.forward_log_z+1e-8,"Beam partition exceeds complete model");
      require(stats.peak_layer_cells<=beam,"Beam memory limit exceeded");
      std::vector<double> columns(y.size(),0);
      for(const auto& row:sparse) {
        double sum=0;
        for(const auto& p:row) {require(p.second>=0 && p.second<=1,"Invalid posterior");sum+=p.second;columns[p.first]+=p.second;}
        require(sum<=1+1e-6,"Posterior row exceeds one");
      }
      for(double sum:columns)require(sum<=1+1e-6,"Posterior column exceeds one");
      if(beam>=std::min(x.size(),y.size())+1)
        require(sparse==exact,"A sufficient beam must reproduce the unpruned calculation");
    }
  }
  // The reference uses polynomial Fast_LogPlus/Fast_Exp and float arithmetic;
  // the new implementation uses double log1p/exp. Compare their posteriors,
  // not bit patterns. Omitting the second affine family fails this test.
  require(largest_error<0.002,"Full-model posterior differs from bundled CONTRAlign");
  require(second_gap_visits>10,"Second affine-gap states are inactive");
  for(size_t length:{1000u,2000u,4000u}) {
    std::string x(length,'A'),y(length,'U');MP p;LinearCONTRAlign5::Diagnostics d;
    model.calculate(x,y,20,.01,p,&d);
    require(d.retained_cells<=20*(2*length+1),"Retained cells are not linear at fixed beam");
    std::cout<<"scaling length="<<length<<" cells="<<d.retained_cells<<" peak_layer="<<d.peak_layer_cells<<'\n';
  }
  std::cout<<"PASS pairs="<<examples.size()<<" probabilities="<<values
           <<" max_absolute_error="<<largest_error<<" mean_absolute_error="<<total_error/values
           <<" second_gap_expected_visits="<<second_gap_visits<<'\n';
}
