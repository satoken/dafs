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


/////////////////////////////////////////////////////////////////
// FOR DEBUG
/////////////////////////////////////////////////////////////////

  /*
  uint sl[] = {3,3,3};
  VU seqlen2(sl, sl+3);
  Align_Graph g(seqlen2);
  
  // graph_1 (non-DAG)
  g.add_edge( node(0,0), node(1,0) );
  g.add_edge( node(0,0), node(2,0) );
  g.add_edge( node(1,0), node(2,0) );
  
  g.add_edge( node(0,1), node(1,1) );
  g.add_edge( node(0,1), node(2,1) );
  g.add_edge( node(1,1), node(2,2) );

  g.add_edge( node(0,2), node(1,2) );  
  g.add_edge( node(0,2), node(2,2) );
  
  
  // graph_2 (DAG)
  g.add_edge( node(0,0), node(1,0) );
  g.add_edge( node(0,0), node(2,0) );
  g.add_edge( node(1,0), node(2,0) );
  
  g.add_edge( node(0,1), node(1,1) );
  g.add_edge( node(0,1), node(2,1) );
  g.add_edge( node(1,1), node(2,1) );

  g.add_edge( node(0,2), node(1,2) );
  g.add_edge( node(0,2), node(2,2) );
  g.add_edge( node(1,2), node(2,2) );
  */

  /*
  // graph_3 (non-DAG)
  uint sl[] = {1,1,4,2,2};
  VU seqlen2(sl, sl+5);
  Align_Graph g(seqlen2);

  g.add_edge( node(0,0), node(1,0) );
  g.add_edge( node(0,0), node(2,1) );
  g.add_edge( node(0,0), node(3,0) );
  g.add_edge( node(2,0), node(3,1) );
  g.add_edge( node(2,0), node(4,1) );
  g.add_edge( node(3,1), node(4,1) );
  //g.add_edge( node(0,0), node(4,0) );
  g.add_edge( node(1,0), node(2,3) );
  */
  
  
  /////////////////////////////////////////////////////////////////
  // 
  /////////////////////////////////////////////////////////////////
