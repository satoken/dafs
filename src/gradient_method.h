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

#pragma once

#include "typedefs.h"
#include <cmath>
#include <map>

//Align_Graph
#include "alignment_graph.h"

class GradientMethod
{
private:
  float eta_;
  std::map<VN, std::tuple<float, float, float>> clips_;
  std::map<VE, std::tuple<float, float, float>> cycles_;

public:
  void add_clip(const VVN &clips)
  {
    for (const auto &clip : clips)
    {
      clips_.insert({clip, {0.0f, 0.0f, 0.0f}});
    }
  }

  void add_cycle(const VVE &cycles)
  {
    for (const auto &cycle : cycles)
    {
      cycles_.insert({cycle, {0.0f, 0.0f, 0.0f}});
    }
  }

  auto num_clips() const 
  {
    return clips_.size();
  }

  auto num_cycles() const
  {
    return cycles_.size();
  }

  void update(VVVVF &lambda, Align_Graph &g, uint t)
  {
    for (auto &[clip, vals] : clips_)
    {
      const node &inner = clip[0];
      const node &outer_1 = clip[1];
      const node &outer_2 = clip[2];
      // calculate gradient
      const bool z01 = g.isAdjacent(inner, outer_1);
      const bool z02 = g.isAdjacent(inner, outer_2);
      const bool z12 = g.isAdjacent(outer_1, outer_2);
      const int grad = (1 - z01 - z02 + z12);

      // update lagmul
      auto &[lagmul, state1, state2] = vals;
      //auto delta = adagrad(grad, state1);
      auto delta = adam(grad, state1, state2, t);
      //auto delta = 0.01f * grad;
      lagmul = std::max(0.0f, lagmul - delta);
      // update lambda
      tie(lambda, edge(inner, outer_1)) -= lagmul;
      tie(lambda, edge(inner, outer_2)) -= lagmul;
      tie(lambda, edge(outer_1, outer_2)) += lagmul;
    }

    for (auto &[cycle, vals] : cycles_)
    {
      // calculate gradient
      int grad = cycle.size() - 1;
      for (uint j = 0; j < cycle.size(); j++)
      {
        const auto [n1, n2] = cycle[j];
        grad -= g.isAdjacent(n1, n2);
      }
      // update lagmul
      auto &[lagmul, state1, state2] = vals;
      //auto delta = adagrad(grad, state1);
      auto delta = adam(grad, state1, state2, t);
      //auto delta = 0.01f * grad;
      lagmul = std::max(0.0f, lagmul - delta);

      // update lambda
      for (uint j = 0; j < cycle.size(); j++)
      {
        tie(lambda, cycle[j]) -= lagmul;
      }
    }
  }

  GradientMethod(float eta)
  {
    eta_ = eta;
  }

private:
  float adagrad(float grad, float &state_1)
  {
    const float eps = 1e-8;
    state_1 += grad * grad;
    const float eta = eta_ / (sqrt(state_1) + eps);
    return eta * grad;
  }

  float adam(float grad, float &state_1, float &state_2, uint t)
  {
    const float alpha = (eta_ == 0.5) ? 0.001 : eta_;
    const float beta = 0.9;
    const float gamma = 0.999;
    const float eps = 1e-8;
    state_1 += (1 - beta) * (grad - state_1);
    state_2 += (1 - gamma) * (grad * grad - state_2);
    const float m = state_1 / (1 - std::pow(beta, t));
    const float v = state_2 / (1 - std::pow(gamma, t));
    return alpha * m / (std::sqrt(v) + eps);
  }

  float &tie(VVVVF &lambda, const edge &e)
  {
    auto [n1, n2] = e;
    if (n1.first > n2.first)
      swap(n1, n2);
    const auto [k1, i1] = n1;
    const auto [k2, i2] = n2;
    return lambda[k1][k2][i1][i2];
  }
};