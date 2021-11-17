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

#include "alignment_graph.h"

#include <set>
#include <cassert>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/bellman_ford_shortest_paths.hpp>

#include "spdlog/spdlog.h"

// class: Align_Graph

uint Align_Graph::
    get_node_id(node n) const
{
  const uint k = n.first, i = n.second;
  return node_id_[k][i];
}

bool Align_Graph::
    isAdjacent(uint nid_1, uint nid_2) const
{
  return util::include(adjacency_list_[nid_1], nid_2);
}

bool Align_Graph::
    isAdjacent(node n1, node n2) const
{
  const uint nid_1 = get_node_id(n1);
  const uint nid_2 = get_node_id(n2);
  return isAdjacent(nid_1, nid_2);
}

void Align_Graph::
    set_g()
{
  g_.clear();
  for (uint nid_1 = 0; nid_1 < adjacency_list_.size(); nid_1++)
    for (auto nid_2 : adjacency_list_[nid_1]) 
      g_.add_edge(nid_1, nid_2);
}

void Align_Graph::
    set_components()
{
  component_id_ = g_.get_components(node_num_);
  component_num_ = component_id_[util::max(component_id_)] + 1;

  component_list_.clear();
  component_list_.resize(component_num_, VU());
  for (auto node_id = 0; node_id < component_id_.size(); node_id++)
  {
    auto compo_id = component_id_[node_id];
    component_list_[compo_id].push_back(node_id);
  }
}

void Align_Graph::
    set_cog()
{
  // DFS to detect cycles
  cog_.clear();
  for (uint k = 0; k < node_id_.size(); k++)
  {
    const uint L = node_id_[k].size();
    for (uint i1 = 0; i1 < L; i1++)
      for (uint i2 = i1 + 1; i2 < L; i2++)
      {
        const uint nid_1 = node_id_[k][i1], nid_2 = node_id_[k][i2];
        const uint cid_1 = component_id_[nid_1], cid_2 = component_id_[nid_2];
        // if cid_1 == cid_2, then this is MC-I
        //if (cid_1 != cid_2)
          cog_.add_edge(cid_1, cid_2);
      }
  }
  util::unique(cog_.edge_list_);
}

void Align_Graph::
    set_cog_cycle()
{
  cog_cycles_ = cog_.get_cycles(component_num_);
}

void Align_Graph::
    add_edge(node n1, node n2)
{
  const uint nid_1 = get_node_id(n1);
  const uint nid_2 = get_node_id(n2);
  adjacency_list_[nid_1].push_back(nid_2);
  adjacency_list_[nid_2].push_back(nid_1);
}

VVVU Align_Graph::
    get_all_edges() const
{
  VVVU z(M_, VVU(M_));
  for (uint k1 = 0; k1 != node_id_.size(); ++k1)
    for (uint k2 = 0; k2 != node_id_.size(); ++k2)
      z[k1][k2].resize(node_id_[k1].size(), -1u);

  for (auto nid_1 = 0; nid_1 != adjacency_list_.size(); nid_1++)
  {
    const auto [k1, i1] = node_list_[nid_1];
    for (auto nid_2 : adjacency_list_[nid_1])
    {
      const auto [k2, i2] = node_list_[nid_2];
      z[k1][k2][i1] = i2;
    }
  }

  return z;
}


void Align_Graph::
    configure()
{
  set_g();
  set_components();
  set_cog();
  set_cog_cycle();
  // clips = detect_clips();
  // detect_mixedCycles();
  // cycles_2 = detect_cycles();
#if 0
  spdlog::debug("n_nodes: {}, n_cycles: {}", node_list_.size(), cog_cycles_.size()); 
  for (auto cycle : cog_cycles_)
  {
    uint n_node = 0;
    for (auto cid : cycle) 
      n_node += component_list_[cid].size();
    spdlog::debug("  n_nodes in cycle: {}", n_node);
  }
#endif
}

VVN Align_Graph::
    detect_clips() const
{
  VVN clips;
  for (uint nid = 0; nid < node_list_.size(); nid++)
  {
    const auto [k, i] = node_list_[nid];
    const VU &adj = adjacency_list_[nid];
    if (adj.size() >= 2)
    {
      for (auto ii1 = 0; ii1 != adj.size(); ii1++)
      {
        const uint nid_1 = adj[ii1];
        const auto [k1, i1] = node_list_[nid_1];
        for (auto ii2 = ii1 + 1; ii2 != adj.size(); ii2++)
        {
          const uint nid_2 = adj[ii2];
          if (!util::include(adjacency_list_[nid_2], nid_1))
          {
            const auto [k2, i2] = node_list_[nid_2];
            VN clip;
            clip.emplace_back(k, i);
            if (k1 < k2)
            {
              clip.emplace_back(k1, i1);
              clip.emplace_back(k2, i2);
            }
            else
            {
              clip.emplace_back(k2, i2);
              clip.emplace_back(k1, i1);
            }
            clips.push_back(clip);
          }
        }
      }
    }
  }
  return clips;
}

ALN Align_Graph::
    get_alignmentColumns(const VVVVF& p_z)
{
  // Topological_sort
  VU column_order = cog_.get_topological_order(component_num_);
#if 0
  if (column_order.empty())
  {
    spdlog::warn("Cycles are detected. Modifying alignment graph.");
    destroy_wrongEdges(p_z);
    VVU cog_cycles = cog_.get_cycles(component_num_);
    column_order = cog_.get_topological_order(component_num_);
  }
#endif

  // convert to ALN
  const uint M = node_id_.size();
  ALN aln(M);
  for (uint k = 0; k < M; k++)
  {
    aln[k].first = k;
    aln[k].second.resize(column_order.size(), false);
  }
  for (uint aln_col = 0; aln_col < column_order.size(); aln_col++)
  {
    const uint cid = column_order[aln_col];
    for (auto nid : component_list_[cid])
    {
      const auto [k, i] = node_list_[nid];
      aln[k].second[aln_col] = true;
    }
  }

  return aln;
}

VVE Align_Graph::
    detect_cycles() const
{
  struct EdgeProperties
  {
    int weight;
  };
  typedef boost::adjacency_list<boost::vecS,
                                boost::vecS,
                                boost::directedS,
                                boost::no_property,
                                EdgeProperties>
      Graph;

  assert(node_num_ == node_list_.size());
  assert(node_list_.size() == adjacency_list_.size());
  const auto start_node_id = node_list_.size(); // nid of start node = node_list_.size()
  Graph g(node_list_.size() + 1);
  const auto N = boost::num_vertices(g); // # of nodes

  const auto K = node_id_.size(); // # of sequences
  for (auto k = 0; k < K; ++k)
  {
    boost::add_edge(start_node_id, node_id_[k][0], g);
    for (auto i = 1; i < node_id_[k].size(); ++i)
      boost::add_edge(node_id_[k][i - 1], node_id_[k][i], g);
  }

  for (auto nid_1 = 0; nid_1 < adjacency_list_.size(); nid_1++)
    for (auto nid_2 : adjacency_list_[nid_1])
      boost::add_edge(nid_1, nid_2, g);

  auto weight_pmap = boost::get(&EdgeProperties::weight, g);
  for (auto [ei, ei_end] = boost::edges(g); ei != ei_end; ++ei)
  {
    auto s = boost::source(*ei, g);
    auto t = boost::target(*ei, g);
    weight_pmap[*ei] = s == start_node_id || node_list_[s].first == node_list_[t].first ? -(K + 1) : 1;
  }

  std::vector<int> distance(N, std::numeric_limits<short>::max());
  distance[start_node_id] = 0;
  std::vector<std::size_t> parent(N);
  for (auto i = 0; i < parent.size(); ++i)
    parent[i] = i;

  // spdlog::debug("start bellman_ford");
  bool r = boost::bellman_ford_shortest_paths(g, N,
                                              boost::weight_map(weight_pmap)
                                                  .distance_map(&distance[0])
                                                  .predecessor_map(&parent[0]));
  // spdlog::debug("finish bellman_ford: {}", r);

  VVE cycles;
  if (r) return cycles; // success for breaking all loops

  std::set<std::pair<uint,uint>> finished;
  for (auto i=0; i<parent.size(); ++i) 
  {
    std::vector<uint> visited;
    visited.push_back(i);
    for (auto j = i; j != start_node_id; j = parent[j])
    {
      const uint nid_1 = j;         // target
      const uint nid_2 = parent[j]; // source
      if (node_list_[nid_1].first != node_list_[nid_2].first) // alignment edge
      {
        auto r = finished.insert(std::minmax(nid_1, nid_2));
        if (!r.second) // finished edge
          break;
      }

      auto x = std::find(std::begin(visited), std::end(visited), nid_2) - std::begin(visited);
      if (x>=std::size(visited))
      {
        visited.push_back(nid_2);
      }
      else
      {
        // loop detected
        // spdlog::debug("loop size={}", visited.size()-x);
        visited.push_back(nid_2);
        VE cycle;
        for (auto v = x; v+1 < std::size(visited); ++v)
        {
          const auto [k1, i1] = node_list_[visited[v]];
          const auto [k2, i2] = node_list_[visited[v+1]];
          if (k1 != k2) 
          {
            // spdlog::debug("({},{})-({},{})", k1, i1, k2, i2);
            cycle.emplace_back(node(k1, i1), node(k2, i2));
          }
        }
        cycles.push_back(cycle);
        break;
      }
    }
  }
  std::sort(std::begin(cycles), std::end(cycles));

  return cycles;
}

Align_Graph::
    Align_Graph(VU &seqlen)
{
  M_ = seqlen.size();
  util::insert(seqlen_, seqlen);

  node_num_ = util::accumulate(seqlen);

  const uint M = seqlen.size();
  for (uint k = 0; k < M; k++)
    for (uint i = 0; i < seqlen[k]; i++)
      node_list_.push_back(node(k, i));

  node_id_.resize(M);
  for (uint k = 0; k < M; k++)
    node_id_[k].resize(seqlen[k]);
  for (uint k = 0, c = 0; k < M; k++)
    for (uint i = 0; i < node_id_[k].size(); i++)
      node_id_[k][i] = c++;

  adjacency_list_.resize(node_num_);
}
