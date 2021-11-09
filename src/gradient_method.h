/*
 * Copyright (C) 2016 Takuro Anamizu
 *
 * This file is part of DMSA.
 *
 * DMSA is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * DMSA is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with DMSA.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef __INC_TOPOLOGICAL_SORT_H__
#define __INC_TOPOLOGICAL_SORT_H__

#include "typedefs.h"
#include <cmath>

//Align_Graph
#include "alignment_graph.h"


class GradientMethod{
private:
  float eta_;
  VVN clips_;
  VVE cycles_;
  // [Clip]
  VF lagmul_clip_; // Lagrangian Multiplier
  VF state1_clip_, state2_clip_;
  // [Cycle]
  VF lagmul_cycle_;
  VF state1_cycle_, state2_cycle_;
public:
  void add_clip(VVN& clips){
    for(uint i=0; i<clips.size(); i++){
      VN clip = clips[i];
      if( !util::include(clips_, clip) ){
	clips_.push_back(clip);
	lagmul_clip_.push_back(0.0);
	state1_clip_.push_back(0);
	state2_clip_.push_back(0);
      }
    }
  }
  void add_cycle(VVE& cycles){
    for(uint i=0; i<cycles.size(); i++){
      VE cycle = cycles[i];
      if( !util::include(cycles_, cycle) ){
	cycles_.push_back(cycle);
	lagmul_cycle_.push_back(0.0);
	state1_cycle_.push_back(0);
	state2_cycle_.push_back(0);
      }
    }
  }

  float adagrad(float grad, float& state_1){
    const float eps = 1e-8;
    state_1 += grad*grad;
    const float eta = eta_ / (sqrt(state_1) +eps);
    return eta * grad;
  }
  float adam(float grad, float& state_1, float& state_2, float t){
    const float alpha = (eta_==0.5) ? 0.001 : eta_;
    const float beta = 0.9;
    const float gamma = 0.999;
    const float eps = 1e-8;
    state_1 += (1-beta) * (grad -state_1);
    state_2 += (1-gamma) * (grad*grad -state_2);
    const float m = state_1 /(1 -std::pow(beta, t));
    const float v = state_2 /(1 -std::pow(gamma, t));
    return alpha * m/(std::sqrt(v) +eps);
  }

  void update(edge e, float lagmul, VVVVF& lambda, bool sgn){
    node n1 = e.first, n2 = e.second;
    if(n1.first > n2.first) swap(n1,n2);
    const uint k1 = n1.first, i1 = n1.second;
    const uint k2 = n2.first, i2 = n2.second;
    lambda[k1][k2][i1][i2] += lagmul * (sgn ? 1.0 : -1.0);
    // lambda[k1][k2][i1][i2] -= 0.5;
  }
  
  void update(VVVVF& lambda, Align_Graph& g, float t){
    //DEBUG: 
    /*
    std::cout << " clip(total): " << clips_.size() << std::endl
	      << " cycle(total): " << cycles_.size() << std::endl;
    */
    for(uint i=0; i<clips_.size(); i++){
      VN clip = clips_[i];
      const node inner = clip[0];
      const node outer_1 = clip[1];
      const node outer_2 = clip[2];
      // calculate gradient
      const bool z01 = g.isAdjacent(inner, outer_1);
      const bool z02 = g.isAdjacent(inner, outer_2);
      const bool z12 = g.isAdjacent(outer_1, outer_2);
      const int grad = (1 -z01 -z02 +z12);
      // update lagmul
      lagmul_clip_[i] = std::max(lagmul_clip_[i]
				 //-adagrad(grad, state1_clip_[i]) ,0.0f);
				 -adam(grad, state1_clip_[i], state2_clip_[i], t) ,0.0f);
      // update lambda
      float lagmul = lagmul_clip_[i];
      update( edge(inner, outer_1), lagmul, lambda, false);
      update( edge(inner, outer_2), lagmul, lambda, false);
      update( edge(outer_1, outer_2), lagmul, lambda, true);
      //
    }
    for(uint i=0; i<cycles_.size(); i++){
      VE cycle = cycles_[i];
      // calculate gradient
      int grad = cycle.size() -1;
      for(uint j=0; j<cycle.size(); j++){
	const node n1 = cycle[j].first, n2 = cycle[j].second;
	grad -= g.isAdjacent(n1,n2);
      }
      // update lagmul
      lagmul_cycle_[i] = std::max(lagmul_cycle_[i]
				  //-adagrad(grad, state1_cycle_[i]) ,0.0f);
				  -adam(grad, state1_cycle_[i], state2_cycle_[i], t) ,0.0f);
      // update lambda
      for(uint j=0; j<cycle.size(); j++){
	update(cycle[j], lagmul_cycle_[i], lambda, false);
      }
      //
    }
  }
  
  GradientMethod(float eta){
    eta_ = eta;
  }
};


#endif  //  __INC_TOPOLOGICAL_SORT_H__

// Local Variables:
// mode: C++
// End: