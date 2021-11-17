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
}

VVN Align_Graph::
    detect_clips() const
{
  VVN clips;
  for (uint nid = 0; nid < node_list_.size(); nid++)
  {
    const node &n0 = node_list_[nid];
    const auto [k, i] = n0;
    const VU &adj = adjacency_list_[nid];
    if (adj.size() >= 2)
    {
      for (auto ii1 = 0; ii1 != adj.size(); ii1++)
        for (auto ii2 = ii1 + 1; ii2 != adj.size(); ii2++)
        {
          const uint nid_1 = adj[ii1], nid_2 = adj[ii2];
          if (!util::include(adjacency_list_[nid_2], nid_1))
          {
            const auto [k1, i1] = node_list_[nid_1];
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
  return clips;
}

void Align_Graph::
    detect_mixedCycles()
{
  cycles_1.clear();
  cycles_2.clear();
  //VVE cycles;

  for (uint i = 0; i < cog_cycles_.size(); i++)
  {
    VU cog_cycle = cog_cycles_[i];
    const uint size = cog_cycle.size();

    // MC-I
#if 0
    if (size == 1)
    {
      VVE part_of_cycles = get_cycleInComponent(component_list_[i]);
      util::insert(cycles_1, part_of_cycles);
    }
#endif

#if 0
    std::cout << i << ": ";
    for (auto c: cog_cycle) std::cout << c << ", ";
    std::cout << std::endl;
#endif

    // MC-II
    if (size >= 2 && size <= M_)
    {
      VVE part_of_cycles = get_cycleAmongComponents(cog_cycle);
      util::insert(cycles_2, part_of_cycles);
    }
  }

  //return cycles;
}

ALN Align_Graph::
    get_alignmentColumns(const VVVVF& p_z)
{
  // Topological_sort
  VU column_order = cog_.get_topological_order(component_num_);

  if (column_order.empty())
  {
    spdlog::warn("Cycles are detected. Modifying alignment graph.");
    destroy_wrongEdges(p_z);
    VVU cog_cycles = cog_.get_cycles(component_num_);
    column_order = cog_.get_topological_order(component_num_);
  }

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

//

VVE Align_Graph::
    get_cycleAmongComponents(const VU &cog_cycle)
{
  const uint compo_num = cog_cycle.size();
  VVVU compo_list(compo_num + 1, VVU(M_));
  for (uint j = 0; j < compo_list.size(); j++)
  {
    const uint cid = cog_cycle[j % cog_cycle.size()];
    for (auto nid : component_list_[cid])
    {
      const auto [k, i] = node_list_[nid];
      compo_list[j][k].push_back(i);
    }
  }

  VVE cycles_axis; //component transition

  struct InnerFunction
  {
    static void
    function(const VVVU &compo_list, uint i, VE mc, VVE &mixed_cycles, const Align_Graph& g)
    {
      if (i == compo_list.size() - 1)
      {
        mixed_cycles.push_back(mc);
        mc.clear();
        return;
      }

      const VVU &compo_1 = compo_list[i];
      const VVU &compo_2 = compo_list[i + 1];
      const uint M = compo_1.size();
      for (uint k = 0; k < M; k++)
      {
        for (uint j1 = 0; j1 < compo_1[k].size(); j1++)
          for (uint j2 = 0; j2 < compo_2[k].size(); j2++)
          {
            if (j1 > 0 || j2 > 0)
              continue;
            const uint i1 = compo_1[k][j1];
            const uint i2 = compo_2[k][j2];
            if (i1 < i2)
            {
#if 0
              auto nid_1 = g.node_id_[k][i1];
              auto nid_2 = g.node_id_[k][i2];
              auto cid_1 = g.component_id_[nid_1];
              auto cid_2 = g.component_id_[nid_2];
              spdlog::debug("make graph: {}: nid={}, cid={}, ({},{}) -> nid={}, cid={}, ({},{})", i, nid_1, cid_1, k, i1, nid_2, cid_2, k, i2);
#endif
              mc.emplace_back(node(k, i1), node(k, i2));
              function(compo_list, i + 1, mc, mixed_cycles, g);
              mc.pop_back();
            }
          }
      }
    }
  };
  //spdlog::debug("get cycle: size={}", cog_cycle.size());
  VE mc;
  InnerFunction::function(compo_list, 0, mc, cycles_axis, *this);
  //spdlog::debug("cycles_axis size={}", cycles_axis.size());

  VVE cycles;
  for (uint i = 0; i < cycles_axis.size(); i++)
  {
    VE a = cycles_axis[i], c;

#if 0
    spdlog::debug("path size={}", a.size());
    for (uint j = 0; j < a.size() - 1; j++)
    {
      node n1 = a[j].second, n2 = a[j + 1].first;
      auto nid_1 = get_node_id(n1);
      auto nid_2 = get_node_id(n2);
      auto cid_1 = component_id_[nid_1];
      auto cid_2 = component_id_[nid_2];
      spdlog::debug("path: {}, {}: nid={}, cid={}, ({},{}) -> nid={}, cid={}, ({},{}), adjacent={}", i, j, nid_1, cid_1, n1.first, n1.second, nid_2, cid_2, n2.first, n2.second, isAdjacent(n1, n2));
    }
#endif

    for (uint j = 0; j < a.size() - 1; j++)
    {
      node n1 = a[j].second, n2 = a[j + 1].first;
      if (n1==n2)
        continue;
      if (isAdjacent(n1, n2))
        c.emplace_back(n1, n2);
      else
      {
#if 1
        c.clear();
#else
        c.clear();
        const uint nid_1 = get_node_id(n1), nid_2 = get_node_id(n2);
        spdlog::debug("nid={}, cid={}: ({},{})-> nid={}, cid={}: ({},{})", nid_1, component_id_[nid_1], n1.first, n1.second, nid_2, component_id_[nid_2], n2.first, n2.second);
        assert(component_id_[nid_1]==component_id_[nid_2]);
        VU path = g_.get_shortestPath(nid_2, nid_1, node_num_);
        for (uint j=0; j<path.size()-1; j++)
        {
          spdlog::debug("shortest path: ({},{})->({},{}): ({},{})->", n1.first, n1.second, n2.first, n2.second, node_list_[path[j]].first, node_list_[path[j]].second);
        }
#if 0
        const node n1 = node_list_[nid_1], n2 = node_list_[nid_2];
        VVN paths = get_paths(n1, n2);
        for (const auto& path : paths)
          for (uint j = 0; j < path.size() - 1; j++)
            c.emplace_back(path[j], path[j+1]);
#endif
#endif
      }
    }
    if (!c.empty())
      cycles.push_back(c);
  }

  // unique cycle_list
  for (uint i = 0; i < cycles.size(); i++)
    util::cycle_sort(cycles[i]);
  util::unique(cycles);

  return cycles;
}

VVE Align_Graph::
    get_cycleInComponent(const VU &compo_nodes)
{
  VVE cycles;

  const uint M = node_id_.size();
  VVU x(M);
  for (auto ii = compo_nodes.cbegin(); ii != compo_nodes.cend(); ii++)
  {
    const uint nid = *ii;
    const node &n = node_list_[nid];
    const uint k = n.first, i = n.second;
    x[k].push_back(nid);
  }
  for (uint k = 0; k < M; k++)
    util::sort(x[k]);

  for (uint k = 0; k < M; k++)
  {
    const uint size = x[k].size();
    for (uint i = 0; i < size; i++)
      for (uint j = i + 1; j < size; j++)
      {
        const uint nid_1 = x[k][i], nid_2 = x[k][j];
        //OBSOLETE: VU path = get_shortestPath( nid_2, nid_1, node_num_ );
        const node n1 = node_list_[nid_1], n2 = node_list_[nid_2];
        VVN paths = get_paths(n1, n2);
        //DEBUG: util::print(path);
        for (uint i = 0; i < paths.size(); i++)
        {
          const VN &path = paths[i];
          VE cycle;
          for (uint j = 0; j < path.size() - 1; j++)
          {
            const node n1 = path[j], n2 = path[j + 1];
            cycle.push_back(edge(n1, n2));
          }
          cycles.push_back(cycle);
        }
      }
  }

  // unique cycle_list
  for (uint i = 0; i < cycles.size(); i++)
    util::cycle_sort(cycles[i]);
  util::unique(cycles);

  return cycles;
}

VVN Align_Graph::
    get_paths(node n1, node n2) const
{
  struct InnerFunction
  {
    static void
    function(node m, node e, VB v, VN p, VVN &paths, const Align_Graph &g)
    {
      p.push_back(m);
      if (m == e)
      {
        paths.push_back(p);
        p.clear();
        return;
      }

      const uint k = m.first;
      if (v[k] == true)
        return;
      else
        v[k] = true;

      const VU &l = g.adjacency_list_[g.get_node_id(m)];
      for (uint i = 0; i < l.size(); i++)
      {
        node l_i = g.node_list_[l[i]];
        function(l_i, e, v, p, paths, g);
      }
    }
  };
  VB v(M_);
  VN p;
  VVN paths;
  InnerFunction::function(n1, n2, v, p, paths, *this);

  return paths;
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

void Align_Graph::
    destroy_wrongEdges(const VVVVF &p_z)
{
  //DEBUG:
#if 0
  for(uint i=0; i<cog_cycles_.size(); i++)
    util::print(cog_cycles_[i]);
#endif

  assert(node_num_ == node_list_.size());

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

  for (uint t=0; t<1000; ++t)
  {
    spdlog::debug("remove loops: {}", t);
    assert(node_list_.size() == adjacency_list_.size());
    const auto start_node_id = node_list_.size(); // nid of start node = node_list_.size()
    Graph g(node_list_.size() + 1);
    const auto N = boost::num_vertices(g);

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

    spdlog::debug("start bellman_ford");
    bool r = boost::bellman_ford_shortest_paths(g, N,
                                                boost::weight_map(weight_pmap)
                                                    .distance_map(&distance[0])
                                                    .predecessor_map(&parent[0]));

    spdlog::debug("finish bellman_ford: {}", r);
    if (r) break; // success for breaking all loops

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
          spdlog::debug("loop size={}", visited.size()-x);
          visited.push_back(nid_2);
          // loop detected
          float min_p = 1.;
          auto min_p_idx = visited.size();
          for (auto v = x; v+1 < std::size(visited); ++v)
          {
            const auto [k1, i1] = node_list_[visited[v]];
            const auto [k2, i2] = node_list_[visited[v+1]];
            if (k1 != k2) 
            {
              auto p = p_z[k1][k2][i1][i2];
              spdlog::debug("p(({},{})-({},{}))={}", k1, i1, k2, i2, p);
              if (min_p_idx == visited.size() || min_p > p) 
              {
                min_p_idx = v;
                min_p = p;
              }
            }

          }
          assert(min_p_idx < visited.size());
          const auto [k1, i1] = node_list_[visited[min_p_idx]];
          const auto [k2, i2] = node_list_[visited[min_p_idx+1]];
          spdlog::debug("remove ({},{})-({},{})", k1, i1, k2, i2);
          util::remove(adjacency_list_[visited[min_p_idx]], visited[min_p_idx+1]);
          util::remove(adjacency_list_[visited[min_p_idx+1]], visited[min_p_idx]);
          break;
        }
      }
    }
  }

#if 0
  // Mixed Cycle
  VU components; //list of cid
  for (uint i = 0; i < cog_cycles_.size(); i++)
    util::insert(components, cog_cycles_[i]);
  util::unique(components);

  for (auto cid : components)
  {
    const VU &compo_list = component_list_[cid];
    for (auto nid_1 : component_list_[cid])
      for (auto nid_2 : component_list_[cid])
      {
        const uint k1 = node_list_[nid_1].first, k2 = node_list_[nid_2].first;
        if ((k1 - k2) * (k1 - k2) != 1)
        {
          util::remove(adjacency_list_[nid_1], nid_2);
          util::remove(adjacency_list_[nid_2], nid_1);
        }
      }
  }
#endif

#if 0
  VVE cycles;
  util::insert(cycles, cycles_1);
  util::insert(cycles, cycles_2);
  for (uint i = 0; i < cycles.size(); i++)
  {
    VE cycle = cycles[i];
    VF vp;
    for (uint j = 0; j < cycle.size(); j++)
    {
      const edge e = cycle[j];
      node n1 = e.first, n2 = e.second;
      if (n1.first > n2.first)
        swap(n1, n2);
      const uint k1 = n1.first, i1 = n1.second;
      const uint k2 = n2.first, i2 = n2.second;
      vp.push_back(p[k1][k2][i1][i2]);
    }
    const uint j = util::min(vp);
    const edge e = cycle[j];
    const node n1 = e.first, n2 = e.second;
    const uint nid_1 = get_node_id(n1), nid_2 = get_node_id(n2);
    util::remove(adjacency_list_[nid_1], nid_2);
    util::remove(adjacency_list_[nid_2], nid_1);
  }

  //DEBUG:
  for (uint i = 0; i < cycles.size(); i++)
  {
    VE cycle = cycles[i];
    for (uint j = 0; j < cycle.size(); j++)
    {
      edge e = cycle[j];
      node n1 = e.first, n2 = e.second;
      util::print(n1);
      std::cout << "-";
      util::print(n2);
      std::cout << " -> ";
    }
    std::cout << std::endl;
  }
#endif

  configure();
  //std::cout << "DEBUG: COGcycle.size(): " << cog_cycles_.size() << std::endl;
}

//

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
