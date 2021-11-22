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

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "gradient_method.h"

void GradientMethod::
    add_clip(const VVN &clips)
{
  for (const auto &clip : clips)
  {
    clips_.insert({clip, {0.0f, 0.0f, 0.0f}});
  }
}

void GradientMethod::
    add_cycle(const VVE &cycles)
{
  for (const auto &cycle : cycles)
  {
    cycles_.insert({cycle, {0.0f, 0.0f, 0.0f}});
  }
}

auto GradientMethod::
    num_clips() const
{
  return clips_.size();
}

auto GradientMethod::num_cycles() const
{
  return cycles_.size();
}

std::tuple<uint, uint, float>
GradientMethod::
    update(VVVVF &lambda, Align_Graph &g, float lr, uint t)
{
  float g2 = 0.0;
  uint n_violated_clips = 0, n_violated_cycles = 0;
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
    g2 += grad * grad;
    if (grad < 0)
      n_violated_clips++;
  }

  for (auto &[cycle, vals] : cycles_)
  {
    // calculate gradient
    int grad = cycle.size() - 1;
    for (const auto [n1, n2] : cycle)
      grad -= g.isAdjacent(n1, n2);
    g2 += grad * grad;
    if (grad < 0)
      n_violated_cycles++;
  }

  auto eta = eta_ * (lr - lb_) / g2;
  spdlog::debug("LR={}, LB={}, eta={}", lr, lb_, eta);

  std::set<VN> rm_clips;
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
    auto delta = eta * grad;
    lagmul = std::max(0.0f, lagmul - delta);
    // update lambda
    tie(lambda, edge(inner, outer_1)) -= lagmul;
    tie(lambda, edge(inner, outer_2)) -= lagmul;
    tie(lambda, edge(outer_1, outer_2)) += lagmul;

    if (lagmul == 0.0f)
      rm_clips.insert(clip);
  }
  for (const auto clip : rm_clips)
    clips_.erase(clip);

  std::set<VE> rm_cycles;
  for (auto &[cycle, vals] : cycles_)
  {
    // calculate gradient
    int grad = cycle.size() - 1;
    for (const auto [n1, n2] : cycle)
      grad -= g.isAdjacent(n1, n2);
    // update lagmul
    auto &[lagmul, state1, state2] = vals;
    auto delta = eta * grad;
    lagmul = std::max(0.0f, lagmul - delta);

    // update lambda
    for (const auto &e : cycle)
      tie(lambda, e) -= lagmul;

    if (lagmul == 0.0)
      rm_cycles.insert(cycle);
  }
  for (const auto &cycle : rm_cycles)
    cycles_.erase(cycle);

  return {n_violated_clips, n_violated_cycles, eta};
}

GradientMethod::
    GradientMethod(float eta, float lb /*= 0.0*/)
    : eta_(eta), lb_(lb)
{
}

float &GradientMethod::
    tie(VVVVF &lambda, const edge &e)
{
  const auto [n1, n2] = std::minmax(e.first, e.second);
  const auto [k1, i1] = n1;
  const auto [k2, i2] = n2;
  return lambda[k1][k2][i1][i2];
}