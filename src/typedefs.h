/*
 * Copyright (C) 2016 Kengo Sato, Takuro Anamizu
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

#ifndef __INC_TYPEDEFS_H__
#define __INC_TYPEDEFS_H__

#include <vector>

typedef unsigned int uint;

typedef std::vector<float> VF;
typedef std::vector<VF> VVF;
typedef std::vector<VVF> VVVF;
typedef std::vector<VVVF> VVVVF;

typedef std::vector<int> VI;
typedef std::vector<VI> VVI;
typedef std::vector<VVI> VVVI;
typedef std::vector<VVVI> VVVVI;

typedef std::vector<uint> VU;
typedef std::vector<VU> VVU;
typedef std::vector<VVU> VVVU;

typedef std::vector<char> VC;
typedef std::vector<VC> VVC;

typedef std::vector<bool> VB;
  
typedef std::vector<std::pair<uint,float> > SV; // sparse vectors
typedef std::vector<SV> MP;
typedef std::vector<SV> BP;


typedef std::vector< std::pair<uint, std::vector<bool> > > ALN; // alignments


// Utility

#include <algorithm>
#include <iostream>
#include <numeric>
namespace util{

  template <class X>
  void sort(std::vector<X>& v){
    std::sort(v.begin(), v.end());
  }
 
 template <class X>
  void unique(std::vector<X>& v){
    util::sort(v);
    v.erase( unique(v.begin(),v.end()), v.end());
  }
  
  template <class X>
  void reverse(std::vector<X>& v){
    std::reverse(v.begin(), v.end());
  }

  template <class X>
  const X accumulate(const std::vector<X>& v){
    return accumulate(v.cbegin(), v.cend(), 0);
  }

  template <class X>
  uint max(const std::vector<X>& v){
    return max_element(v.cbegin(), v.cend()) - v.cbegin();
  }

  template <class X>
  uint min(const std::vector<X>& v){
    return min_element(v.cbegin(), v.cend()) - v.cbegin();
  }

  template <class X>
  void print(const std::vector<X>& v){
    for(uint i=0; i<v.size(); i++)
      std::cout << v[i] << ", ";
    std::cout << std::endl;
  }

  template <class X>
  void insert(std::vector<X>& v1, std::vector<X>& v2){
    v1.insert(v1.end(), v2.begin(), v2.end());
  }

  template <class X>
  void copy(std::vector<X>& v1, std::vector<X>& v2){
    v2.reserve( v1.size() );
    std::copy( v1.begin(), v1.end(), std::back_inserter(v2) );
  }

  template <class X>
  void cycle_sort(std::vector<X>& v){
    uint i = min_element(v.begin(), v.end()) - v.begin();
    if(i==0) return;
    std::vector<X> w1(v.begin()+i, v.end());
    std::vector<X> w2(v.begin(), v.begin()+i);
    util::insert(w1,w2);
    swap(v, w1);
  }
 
  template <class X>
  void remove(std::vector<X>& v, X value){
    v.erase( std::remove(v.begin(),v.end(),value) , v.end());
  }
  
  template <class X>
  bool include(const std::vector<X>& v, X value){
    return (std::find(v.cbegin(), v.cend(), value) != v.cend());
  }
  
  template <class X>
  uint find(const std::vector<X>& v, X value){
    return std::find(v.cbegin(), v.cend(), value) - v.cbegin();
  }

  template <class X>
  uint rfind(const std::vector<X>& v, X value){
    return v.crend() - std::find(v.crbegin(), v.crend(), value) -1;
  }

  template <class X>
  void fill(std::vector<X>& v, X value){
    std::fill(v.begin(), v.end(), value);
  }

  /*
  //template <class X>
  void print(std::string str){
    std::cout << str << std::endl;
  }
  */
  template <class X1, class X2>
  void print(const std::pair<X1,X2>& p){
    std::cout << '(' << p.first << ',' << p.second << ')';
  }
  
}


#endif

// Local Variables:
// mode: C++
// End:
