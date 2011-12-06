// $Id:$

#ifndef __INC_TYPEDEFS_H__
#define __INC_TYPEDEFS_H__

#include <vector>

typedef unsigned int uint;

typedef std::vector<float> VF;
typedef std::vector<VF> VVF;

typedef std::vector<int> VI;
typedef std::vector<VI> VVI;

typedef std::vector<uint> VU;
typedef std::vector<VU> VVU;

typedef std::vector<char> VC;
typedef std::vector<VC> VVC;
  
typedef std::vector<std::pair<uint,float> > SV; // sparse vectors
typedef std::vector<SV> MP;
typedef std::vector<SV> BP;

#endif
