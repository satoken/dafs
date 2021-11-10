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

#include "typedefs.h"
#include "graph_library.h"

typedef std::pair<uint, uint> node;
typedef std::pair<node, node> edge;
typedef std::vector<node> VN;
typedef std::vector<edge> VE;
typedef std::vector<VN> VVN;
typedef std::vector<VE> VVE;

class Align_Graph
{
private:
  uint M_;
  VU seqlen_;
  //
  uint node_num_;
  VN node_list_;
  VVU node_id_;
  uint get_node_id(node n);
  //
  VVU adjacency_list_;
  bool isAdjacent(uint nid_1, uint nid_2);
  //
  Undirected_Graph g_; //Undirected edges
  void set_g();
  //
  uint component_num_;
  VU component_id_;
  VVU component_list_;
  Directed_Graph cog_; //Component Order Graph
  VVU cog_cycles_;
  void set_components();
  void set_cog();
  void set_cog_cycle();
  VVN get_paths(node n1, node n2);
  //
  void detect_clips();
  void detect_mixedCycles();

public:
  VVN clips;
  VVE cycles_1; //MC-I
  VVE cycles_2; //MC-II
  bool isAdjacent(node n1, node n2);
  VVE get_cycleAmongComponents(const VU& compo_cycles);
  VVE get_cycleInComponent(VU compo_nodes);
  uint cog_cycle_num() { return cog_cycles_.size(); }
  //
  void add_edge(node n1, node n2);
  void configure();
  ALN get_alignmentColumns();
  void destroy_wrongEdges();
  Align_Graph(VU &seqlen);
};
