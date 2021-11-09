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

#include <typedefs.h>

class Graph{
  typedef std::pair<uint,uint> edge;
  typedef std::vector<edge> VE;
 private:
 public:
  VE  edge_list_;
  void add_edge(uint n1, uint n2);
  void clear();
};

class Undirected_Graph : public Graph{
 private:
 public:
 VU get_components(uint node_num);
 VU get_shortestPath(uint n1, uint n2, uint node_num);
};


class Directed_Graph : public Graph{
 private:
 public:
  VVU get_cycles(uint node_num);
  VU get_topological_order(uint node_num);
};


class Mixed_Graph : public Graph{
 private:
 public:
  VVU get_cycles(uint node_num);
};
