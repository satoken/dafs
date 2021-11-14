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

#include "graph_library.h"

//Boost Graph Library
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/graph/depth_first_search.hpp>
#include <boost/graph/topological_sort.hpp>
#include <boost/range/algorithm/for_each.hpp>
#include <boost/range/adaptor/reversed.hpp>
#include <boost/graph/graphviz.hpp>
#include <iterator>
#include <utility>
#include <cassert>
#include "spdlog/spdlog.h"

//class: "Graph"
void Graph::
    add_edge(uint n1, uint n2)
{
  edge_list_.emplace_back(n1, n2);
}

void Graph::
    clear()
{
  edge_list_.clear();
}

//class: "Undirected_Graph"
VU Undirected_Graph::
    get_components(uint node_num) const
{
  typedef boost::adjacency_list<boost::listS,
                                boost::vecS,
                                boost::undirectedS>
      Graph;
  Graph g(edge_list_.begin(), edge_list_.end(), node_num);

  VU component(boost::num_vertices(g));
  const uint component_num = boost::connected_components(g, &component[0]);
  return component;
}

VU Undirected_Graph::
    get_shortestPath(uint n1, uint n2, uint node_num) const
{
  typedef boost::adjacency_list<boost::listS,
                                boost::vecS,
                                boost::undirectedS,
                                boost::no_property,
                                boost::property<boost::edge_weight_t, int>>
      Graph;
  Graph g(edge_list_.begin(), edge_list_.end(), node_num);

  typedef boost::graph_traits<Graph>::vertex_descriptor Vertex;
  const Vertex from = n1;
  const Vertex to = n2;

  std::vector<Vertex> parents(num_vertices(g), -1u); // initialized by unreachable flags
  boost::static_property_map<int> weight(1);
  boost::dijkstra_shortest_paths(g, from, 
    boost::visitor(
      boost::make_dijkstra_visitor(
        boost::record_predecessors(&parents[0], boost::on_edge_relaxed())))
    .weight_map(weight));

  VU path;
  for (Vertex v = to; v != from; v = parents[v])
  {
    assert(v != -1u);
    path.push_back(v);
  }
  path.push_back(from);
  util::reverse(path);
  return path;
}

VVU Directed_Graph::
    get_cycles(uint node_num) const
{
  typedef boost::adjacency_list<boost::vecS,
                                boost::vecS,
                                boost::directedS,
                                boost::no_property,
                                boost::property<boost::edge_color_t, boost::default_color_type>>
      Graph;
  typedef boost::graph_traits<Graph>::vertex_descriptor Vertex;
  typedef boost::graph_traits<Graph>::edge_descriptor Edge;

  Graph g(edge_list_.begin(), edge_list_.end(), node_num);
  
  VVU list_cycle;
  struct detect_loops : public boost::dfs_visitor<>
  {
    VU st;
    VVU& list_cycle;

    detect_loops(VVU& lc) : list_cycle(lc) {}

    void start_vertex(Vertex v, const Graph &g)
    {
      st.push_back(v);
    }

    void tree_edge(Edge e, const Graph &g)
    {
      while (st.back() != source(e, g))
      {
        st.pop_back();
      }
      st.push_back(target(e, g));
    }

    void back_edge(Edge e, const Graph &g)
    {
      while (st.back() != source(e, g))
      {
        st.pop_back();
      }
      const uint t = target(e, g);
      const uint i = util::rfind(st, t);
      VU cycle(st.begin() + i, st.end());
      list_cycle.push_back(cycle);
    }
  } vis(list_cycle);

  boost::depth_first_search(g, boost::visitor(vis));

  //DEBUG: boost::write_graphviz(std::cout, g);

  return list_cycle;
}

VU Directed_Graph::
    get_topological_order(uint node_num) const
{
  typedef boost::adjacency_list<boost::vecS,
                                boost::vecS,
                                boost::directedS,
                                boost::no_property,
                                boost::property<boost::edge_color_t, boost::default_color_type>>
      Graph;

  Graph g(edge_list_.begin(), edge_list_.end(), node_num);

  VU result;
  try
  {
    boost::topological_sort(g, std::back_inserter(result));
  }
  catch (boost::not_a_dag &e)
  {
    std::cout << "This alignment graph isn't DAG." << std::endl;
    //std::cout << e.what() << std::endl;
    //boost::write_graphviz(std::cout, g);
  }
  util::reverse(result);
  return result;
}

//class: "Undirected_Graph"
#if 0
VVU
Mixed_Graph::
get_cycles(uint node_num){

}
#endif