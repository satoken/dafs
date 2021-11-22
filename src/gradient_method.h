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
#include <set>

//Align_Graph
#include "alignment_graph.h"
#include "spdlog/spdlog.h"

class GradientMethod
{
private:
  float eta_;
  float lb_;
  std::map<VN, std::tuple<float, float, float>> clips_;
  std::map<VE, std::tuple<float, float, float>> cycles_;

public:
  void add_clip(const VVN &clips);
  void add_cycle(const VVE &cycles);
  auto num_clips() const;
  auto num_cycles() const;
  std::tuple<uint,uint,float> update(VVVVF &lambda, Align_Graph &g, float lr, uint t);
  GradientMethod(float eta, float lb=0.0); 

private:
  float &tie(VVVVF &lambda, const edge &e);
};