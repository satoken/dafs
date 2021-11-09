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

#include <iostream> //DEBUG:



// class: Align_Graph

uint
Align_Graph::
get_node_id(node n){
  const uint k = n.first, i = n.second;
  return node_id_[k][i];
}


bool
Align_Graph::
isAdjacent(uint nid_1, uint nid_2){
  return util::include( adjacency_list_[nid_1], nid_2 );
}


bool
Align_Graph::
isAdjacent(node n1, node n2){
  const uint nid_1 = get_node_id(n1);
  const uint nid_2 = get_node_id(n2);
  return isAdjacent(nid_1, nid_2);
}


void
Align_Graph::
set_g(){
  g_.clear();
  for(uint nid_1=0; nid_1<adjacency_list_.size(); nid_1++){
    VU x = adjacency_list_[nid_1];
    for(VU::iterator ii=x.begin(); ii!=x.end(); ii++){
      const uint nid_2 = *ii;
      g_.add_edge(nid_1, nid_2);
    }
  }
}


void
Align_Graph::
set_components(){
  component_id_ = g_.get_components(node_num_);
  component_num_ = component_id_[util::max(component_id_)] +1;

  component_list_.clear();
  component_list_.resize(component_num_, VU());
  for(VU::iterator ii=component_id_.begin(); ii!=component_id_.end(); ii++){
    const uint node_id = ii-component_id_.begin();
    const uint compo_id = *ii;
    component_list_[compo_id].push_back(node_id);
  }
}


void
Align_Graph::
set_cog(){
  // DFS to detect cycles
  cog_.clear();
  for(uint k=0; k<node_id_.size(); k++){
    const uint L = node_id_[k].size();
    for(uint i1=0; i1<L; i1++)
      for(uint i2=i1 +1; i2<L; i2++){
	const uint nid_1 = node_id_[k][i1], nid_2 = node_id_[k][i2];
	const uint cid_1 = component_id_[nid_1], cid_2 = component_id_[nid_2];
	//if(cid_1 != cid_2)
	  cog_.add_edge(cid_1, cid_2);
      }
  }
}


void
Align_Graph::
set_cog_cycle(){
  cog_cycles_.clear();
  VVU cog_cycles = cog_.get_cycles(component_num_);
  for(uint i=0; i<cog_cycles.size(); i++){
    VU c = cog_cycles[i];
    //if(c.size() <= M_)
    cog_cycles_.push_back(c);
  }
}


void
Align_Graph::
add_edge(node n1, node n2){
  const uint nid_1 = get_node_id(n1);
  const uint nid_2 = get_node_id(n2);
  adjacency_list_[nid_1].push_back(nid_2);
  adjacency_list_[nid_2].push_back(nid_1);
}


void
Align_Graph::
configure(){
  set_g();
  set_components();
  set_cog();
  set_cog_cycle();
  detect_clips();
  detect_mixedCycles();
}


void
Align_Graph::
detect_clips(){
  clips.clear();
  //VVN clips;
  for(uint nid=0; nid<node_list_.size(); nid++){
    const node n0 = node_list_[nid];
    const uint k = n0.first, i = n0.second;
    VU adj = adjacency_list_[nid];
    const uint degree = adj.size();
    if(degree>=2)
      for(VU::iterator ii1=adj.begin(); ii1!=adj.end(); ii1++)
	for(VU::iterator ii2=ii1 +1; ii2!=adj.end(); ii2++){
	  const uint nid_1 = *ii1, nid_2 = *ii2;
	  node n1 = node_list_[nid_1];
	  node n2 = node_list_[nid_2];
	  if(n1.first > n2.first) swap(n1,n2);
	  const uint k1 = n1.first, i1 = n1.second;
	  const uint k2 = n2.first, i2 = n2.second;
	  if( !util::include(adjacency_list_[nid_2], nid_1) ){
	    VN clip;
	    clip.push_back( node(k,i) );
	    clip.push_back( node(k1,i1) );
	    clip.push_back( node(k2,i2) );
	    clips.push_back( clip );
	  }
	}
  }
  //return clips;
}


void
Align_Graph::
detect_mixedCycles(){
  cycles_1.clear();
  cycles_2.clear();
  //VVE cycles;
  
  for(uint i=0; i<cog_cycles_.size(); i++){
    VU cog_cycle = cog_cycles_[i];
    const uint size = cog_cycle.size();
    // MC-I
    /*
    if(size==1){
      VVE part_of_cycles = get_cycleInComponent( component_list_[i] );
      util::insert( cycles_1, part_of_cycles);
    }
    */
    // MC-II
    if(size >= 2 && size <= M_){
      VVE part_of_cycles = get_cycleAmongComponents( cog_cycle );
      util::insert(cycles_2, part_of_cycles);
    }
  }
    
  //return cycles;
}


ALN
Align_Graph::
get_alignmentColumns(){
  // Topological_sort
  VU column_order = cog_.get_topological_order(component_num_);

  if(column_order.empty()){
    std::cout << "modifying alignment graph .." << std::endl;
    destroy_wrongEdges();
    VVU cog_cycles = cog_.get_cycles(component_num_);
    column_order = cog_.get_topological_order(component_num_);
  }
  else std::cout << "topo_success" << std::endl;
 
 // convert to ALN
 const uint M = node_id_.size();
 ALN aln(M);
 for(uint k=0; k<M; k++){
   aln[k].first = k;
   aln[k].second.resize( column_order.size(), false);
 }
 for(uint aln_col=0; aln_col<column_order.size(); aln_col++){
   const uint cid = column_order[aln_col];
   VU c = component_list_[cid];
   for(VU::iterator jj=c.begin(); jj!=c.end(); jj++){
     const node n = node_list_[*jj];
     const uint k = n.first, i = n.second;
     aln[k].second[aln_col] = true;
   }
 }
 
 return aln;
}


//

VVE
Align_Graph::
get_cycleAmongComponents(VU cog_cycle){
  const uint compo_num = cog_cycle.size();
  
  VVVU compo_list(compo_num, VVU(M_));
  for(uint j=0; j<cog_cycle.size(); j++){
    const uint cid = cog_cycle[j];
    VU cl = component_list_[cid];
    for(VU::iterator jj=cl.begin(); jj!=cl.end(); jj++){
      const uint nid = *jj;
      const node n = node_list_[nid];
      const uint k = n.first, i = n.second;
      compo_list[j][k].push_back(i);
    }
  }    
  
  VVE cycles_axis; //component transition
    
  struct InnerFunction{
    static void
    function(VVVU& compo_list, uint i, VE mc, VVE& mixed_cycles){
      if(i==compo_list.size()-1){
	mixed_cycles.push_back(mc);
	mc.clear();
	return;
      }

      VVU compo_1 = compo_list[i], compo_2 = compo_list[i+1];
      const uint M=compo_1.size();
      for(uint k=0; k<M; k++)
	for(uint j1=0; j1<compo_1[k].size(); j1++)
	  for(uint j2=0; j2<compo_2[k].size(); j2++){
	    if(j1>0 || j2>0)continue;
	    const uint i1 = compo_1[k][j1];
	    const uint i2 = compo_2[k][j2];
	    if(i1<i2){
	      mc.push_back( edge( node(k,i1), node(k,i2) ) );	    
	      function(compo_list, i+1, mc, mixed_cycles);
	    }
	  }
    }
  };
  VE mc;
  InnerFunction::function(compo_list, 0, mc, cycles_axis);
   
  VVE cycles;
  for(uint i=0; i<cycles_axis.size(); i++){
    VE a = cycles_axis[i], c;
    a.push_back( a.front() );

    for(uint j=0; j<a.size()-1; j++){
      node n1 = a[j].second, n2 = a[j+1].first;
      if(isAdjacent(n1, n2))
	c.push_back( edge(n1,n2) );
      else{
	c.clear();
	/*
	const uint nid_1 = get_node_id(n1), nid_2 = get_node_id(n2);
	//OBSOLETE: VU path = get_shortestPath( nid_2, nid_1, node_num_ );
	const node n1 = node_list_[nid_1], n2 = node_list_[nid_2];
	VVN paths = get_paths(n1, n2);
	for(uint i=0; i<paths.size(); i++){
	  VN path = paths[i];
	  VE cycle;
	  for(uint j=0; j<path.size()-1; j++){
	    const node m1 = path[j], m2 = path[j+1];
	    c.push_back( edge(m1,m2) ); 
	  }
	}
	*/
      }
    }
    if(!c.empty())
      cycles.push_back(c);
  }

  // unique cycle_list
  for(uint i=0; i<cycles.size(); i++)
    util::cycle_sort(cycles[i]);
  util::unique(cycles);
    
  return cycles;
}


VVE
Align_Graph::
get_cycleInComponent(VU compo_nodes){
  VVE cycles;

  const uint M = node_id_.size();
  VVU x(M);
  for(VU::iterator ii=compo_nodes.begin(); ii!=compo_nodes.end(); ii++){
    const uint nid = *ii;
    const node n = node_list_[nid];
    const uint k = n.first, i = n.second;
    x[k].push_back(nid);
  }
  for(uint k=0; k<M; k++) util::sort(x[k]);
  
  for(uint k=0; k<M; k++){
    const uint size = x[k].size();
    for(uint i=0; i<size; i++)
      for(uint j=i+1; j<size; j++){
	const uint nid_1 = x[k][i], nid_2 = x[k][j];
	//OBSOLETE: VU path = get_shortestPath( nid_2, nid_1, node_num_ );
	const node n1 = node_list_[nid_1], n2 = node_list_[nid_2];
	VVN paths = get_paths(n1, n2);
	//DEBUG: util::print(path);
	for(uint i=0; i<paths.size(); i++){
	  VN path = paths[i];
	  VE cycle;
	  for(uint j=0; j<path.size()-1; j++){
	    const node n1 = path[j], n2 = path[j+1];
	    cycle.push_back( edge(n1,n2) );
	  }
	  cycles.push_back( cycle );	  
	}
      }
  }
  
  // unique cycle_list
  /*
  for(uint i=0; i<cycles.size(); i++)
    util::cycle_sort(cycles[i]);
  util::unique(cycles);
  */

  return cycles;
}


VVN
Align_Graph::
get_paths(node n1, node n2){
  
 struct InnerFunction{
    static void
    function(node m, node e, VB v, VN p, VVN& paths, Align_Graph& g){
      p.push_back(m);
      if(m==e){
	paths.push_back(p);
	p.clear();
	return;
      }
      
      const uint k = m.first;
      if(v[k]==true) return;
      else v[k] = true;
            
      VU l = g.adjacency_list_[ g.get_node_id(m) ];
      for(uint i=0; i<l.size(); i++){
	node l_i = g.node_list_[ l[i] ];
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


void
Align_Graph::
destroy_wrongEdges(){
  //DEBUG:
  /*
  for(uint i=0; i<cog_cycles_.size(); i++)
    util::print(cog_cycles_[i]);
  */

  // Mixed Cycle
  // /*
  VU components; //list of cid
  for(uint i=0; i<cog_cycles_.size(); i++)
    util::insert(components, cog_cycles_[i]);
  util::unique(components);
  
  for(uint i=0; i<components.size(); i++){
    const uint cid = components[i];
    const VU compo_list = component_list_[cid];
    const uint compo_size = compo_list.size();
    for(uint l1=0; l1<compo_size; l1++)
      for(uint l2=0; l2<compo_size; l2++){
	const uint nid_1 = compo_list[l1], nid_2 = compo_list[l2];
	const uint k1 = node_list_[nid_1].first,k2 = node_list_[nid_2].first;
	if( (k1-k2)*(k1-k2)!=1 ){
	  util::remove( adjacency_list_[nid_1], nid_2 );
	  util::remove( adjacency_list_[nid_2], nid_1 );
	}
      }
  }
  // */
  
  /*
  VVE cycles;
  util::insert(cycles, cycles_1);
  util::insert(cycles, cycles_2);
  for(uint i=0; i<cycles.size(); i++){
    VE cycle = cycles[i];
    VF vp;
    for(uint j=0; j<cycle.size(); j++){
      const edge e = cycle[j];
      node n1 = e.first, n2 = e.second;
      if(n1.first > n2.first) swap(n1,n2);
      const uint k1 = n1.first, i1 = n1.second;
      const uint k2 = n2.first, i2 = n2.second;
      vp.push_back( p[k1][k2][i1][i2] );
    }
    const uint j = util::min(vp);
    const edge e = cycle[j];
    const node n1 = e.first, n2 = e.second;
    const uint nid_1 = get_node_id(n1), nid_2 = get_node_id(n2);
    util::remove( adjacency_list_[nid_1], nid_2 );
    util::remove( adjacency_list_[nid_2], nid_1 );
  }

  //DEBUG:                                                                    
  for(uint i=0; i<cycles.size(); i++){
    VE cycle = cycles[i];
    for(uint j=0; j<cycle.size(); j++){
      edge e = cycle[j];
      node n1 = e.first, n2 = e.second;
      util::print(n1);
      std::cout << "-";
      util::print(n2);
      std::cout << " -> ";
    }
    std::cout << std::endl;
  }
  */

  configure();
  //std::cout << "DEBUG: COGcycle.size(): " << cog_cycles_.size() << std::endl;
    
}


//


Align_Graph::
Align_Graph(VU& seqlen){
  M_ = seqlen.size();
  util::insert(seqlen_ ,seqlen);
  
  node_num_ = util::accumulate(seqlen);
 
  const uint M = seqlen.size();
  for(uint k=0; k<M; k++)
    for(uint i=0; i<seqlen[k]; i++)
      node_list_.push_back(node(k,i));

  node_id_.resize(M);
  for(uint k=0; k<M; k++) node_id_[k].resize( seqlen[k] );
  for(uint k=0,c=0; k<M; k++)
    for(uint i=0; i<node_id_[k].size(); i++)
      node_id_[k][i]=c++;

  adjacency_list_.resize(node_num_);
}