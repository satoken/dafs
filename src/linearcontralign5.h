// Sparse forward/backward for the complete bundled CONTRAlign RNA model.
// The five states, all 24 parameters, permitted transitions, and boundary
// conventions follow contralign/InferenceEngine.ipp. Only lattice cells are
// beam-pruned; no affine-gap state or feature is removed from the model.
#ifndef DAFS_LINEARCONTRALIGN5_H
#define DAFS_LINEARCONTRALIGN5_H

#include <algorithm>
#include <array>
#include <cmath>
#include <cctype>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>
#include "typedefs.h"
#include "contralign/Defaults.ipp"

namespace LinearCONTRAlign5 {
enum State { MATCH, INS_X, INS_Y, INS2_X, INS2_Y, STATES };
struct Diagnostics {
  double forward_log_z=0, backward_log_z=0;
  size_t retained_cells=0, peak_layer_cells=0, pruned_cells=0;
  std::array<double,STATES> expected_state_visits{};
};

inline double log_add(double a,double b) {
  if (a==-std::numeric_limits<double>::infinity()) return b;
  if (b==-std::numeric_limits<double>::infinity()) return a;
  if (a<b) std::swap(a,b);
  return a+std::log1p(std::exp(b-a));
}

inline int nucleotide(char c) {
  switch (std::toupper(static_cast<unsigned char>(c))) {
    case 'A':return 0; case 'C':return 1; case 'G':return 2; case 'U':return 3;
    default:return 4; // Same zero-scored unknown alphabet entry as CONTRAlign.
  }
}

class Model {
  // Read the exact float defaults used by the bundled nonlinear wrapper.
  const std::vector<float> weights_=CONTRALIGN::GetDefaultRNAValues<float>();
  std::array<std::array<float,5>,5> match_{};
  std::array<float,5> insert_{};
  std::array<float,STATES> single_{};
  std::array<std::array<float,STATES>,STATES> transition_{};

public:
  Model() {
    if (weights_.size()!=24) throw std::logic_error("Unexpected CONTRAlign RNA parameter count");
    size_t p=0;
    for (int a=0;a<4;++a) for(int b=a;b<4;++b)
      match_[a][b]=match_[b][a]=weights_[p++];
    for(int a=0;a<4;++a) insert_[a]=weights_[p++];
    single_={weights_[14],weights_[15],weights_[15],weights_[16],weights_[16]};
    transition_[MATCH][MATCH]=weights_[17];
    for (int group=0;group<2;++group) {
      const int x=group?INS2_X:INS_X,y=group?INS2_Y:INS_Y,k=group?21:18;
      transition_[MATCH][x]=transition_[MATCH][y]=transition_[x][MATCH]=transition_[y][MATCH]=weights_[k];
      transition_[x][x]=transition_[y][y]=weights_[k+1];
      transition_[x][y]=transition_[y][x]=weights_[k+2];
    }
  }

  // Parameters/state-space are unchanged at finite beam. A zero beam is an
  // unpruned diagnostic mode; production CLI requires a positive beam.
  void calculate(const std::string& x,const std::string& y,unsigned beam,
                 float cutoff,MP& posterior,Diagnostics* diagnostics=nullptr) const {
    const size_t n=x.size(),m=y.size();
    const double negative=-std::numeric_limits<double>::infinity();
    struct Cell {
      std::array<double,STATES> forward,backward;
      Cell() {forward.fill(-std::numeric_limits<double>::infinity());backward=forward;}
    };
    using Layer=std::map<size_t,Cell>; // i is the key; j=layer-i.
    std::vector<Layer> layers(n+m+1);
    layers[0][0].forward.fill(0.0); // CONTRAlign initializes all five start states.
    Diagnostics stats;
    const auto emit=[&](size_t i,size_t j,int state) {
      const float base=state==MATCH ? match_[nucleotide(x[i-1])][nucleotide(y[j-1])]
          : insert_[nucleotide((state==INS_X||state==INS2_X)?x[i-1]:y[j-1])];
      // Float addition order agrees with ScoreMatch/ScoreInsert in the wrapper.
      return base+single_[state];
    };
    const auto edges=[&](size_t i,size_t j,int from,const auto& visit) {
      if (i<n && j<m && (from==MATCH || i || j)) {
        const float score=emit(i+1,j+1,MATCH)+((i||j)?transition_[from][MATCH]:0.0f);
        visit(i+1,j+1,MATCH,score);
      }
      for(int group=0;group<2;++group) {
        const int sx=group?INS2_X:INS_X,sy=group?INS2_Y:INS_Y;
        if (i<n && (from==sx || (j && (from==MATCH||from==sy)))) {
          const float score=emit(i+1,j,sx)+((i||j)?transition_[from][sx]:0.0f);
          visit(i+1,j,sx,score);
        }
        if (j<m && (from==sy || (i && (from==MATCH||from==sx)))) {
          const float score=emit(i,j+1,sy)+((i||j)?transition_[from][sy]:0.0f);
          visit(i,j+1,sy,score);
        }
      }
    };
    for(size_t step=0;step<=n+m;++step) {
      auto& layer=layers[step];
      if (beam && layer.size()>beam) {
        std::vector<std::pair<double,size_t>> rank;rank.reserve(layer.size());
        for (const auto& entry:layer) {
          double value=negative;
          for(double f:entry.second.forward)value=log_add(value,f);
          rank.emplace_back(value,entry.first);
        }
        // Stable tie policy; every retained cell keeps all five states.
        std::sort(rank.begin(),rank.end(),[](const auto& a,const auto& b) {
          return a.first!=b.first?a.first>b.first:a.second<b.second;
        });
        stats.pruned_cells+=rank.size()-beam;
        for(size_t k=beam;k<rank.size();++k)layer.erase(rank[k].second);
      }
      stats.retained_cells+=layer.size();
      stats.peak_layer_cells=std::max(stats.peak_layer_cells,layer.size());
      for(const auto& entry:layer) {
        const size_t i=entry.first,j=step-i;
        for(int state=0;state<STATES;++state) {
          const double value=entry.second.forward[state];
          if (!std::isfinite(value))continue;
          edges(i,j,state,[&](size_t ni,size_t nj,int target,float score) {
            auto& f=layers[ni+nj][ni].forward[target];f=log_add(f,value+score);
          });
        }
      }
    }
    auto terminal=layers[n+m].find(n);
    if (terminal==layers[n+m].end())throw std::runtime_error("No retained CONTRAlign terminal");
    double z=negative;
    for(double f:terminal->second.forward)z=log_add(z,f);
    if (!std::isfinite(z))throw std::runtime_error("Non-finite CONTRAlign partition function");
    terminal->second.backward.fill(0.0); // No terminal transition feature.
    for(size_t reverse=n+m;reverse>0;--reverse) {
      const size_t step=reverse-1;
      for(auto& entry:layers[step]) {
        const size_t i=entry.first,j=step-i;
        for(int state=0;state<STATES;++state) {
          double value=negative;
          edges(i,j,state,[&](size_t ni,size_t nj,int target,float score) {
            const auto target_cell=layers[ni+nj].find(ni);
            if(target_cell!=layers[ni+nj].end())
              value=log_add(value,score+target_cell->second.backward[target]);
          });
          entry.second.backward[state]=value;
        }
      }
    }
    double backwards=negative;
    for(double b:layers[0].at(0).backward)backwards=log_add(backwards,b);
    if(std::abs(backwards-z)>1e-8*std::max(1.0,std::abs(z)))
      throw std::runtime_error("CONTRAlign forward/backward partition mismatch");
    posterior.assign(n,{});
    for(size_t step=1;step<=n+m;++step) for(const auto& entry:layers[step]) {
      const size_t i=entry.first,j=step-i;
      for(int state=0;state<STATES;++state) {
        const double p=std::exp(entry.second.forward[state]+entry.second.backward[state]-z);
        stats.expected_state_visits[state]+=p;
        if(state==MATCH && i && j) {
          if(!std::isfinite(p)||p>1.0+1e-7)throw std::runtime_error("Invalid CONTRAlign posterior");
          const float clipped=static_cast<float>(std::min(1.0,p));
          if(clipped>cutoff)posterior[i-1].emplace_back(j-1,clipped);
        }
      }
    }
    for(auto& row:posterior)std::sort(row.begin(),row.end());
    stats.forward_log_z=z;stats.backward_log_z=backwards;
    if(diagnostics)*diagnostics=stats;
  }
};
} // namespace LinearCONTRAlign5
#endif
