// FILE EDGE_RELATIONS.H: Declaration of class edge_relations

#if     !defined(_EDGE_RELATIONS_H)
#define _EDGE_RELATIONS_H      1       //flags that this file has been included

#include "geometry.h"
#include "symb.h"

class edge_relations {
public:
  edge_relations() {;}
  edge_relations(symbdata*, int verb=0);
  int coords(int i) const {return coordindex[i];}
  int gen(int i) const {return gens[i];}  // indexed from 1 not 0
  int get_ngens() const {return ngens;}

private:
  symbdata* sd; // provides plusflag, nsymb, symbol(i), symbops
  int verbose, plusflag, nsymb, nsymbx, ngens;
  vector<int> coordindex, gens;

  void edge_relations_1();    // basic edge relations for alpha = 0
  void edge_relations_2();    // extra edge relations for alphas with denom 2
  void edge_pairing(int i);   // edge relation pair, alpha=r/s with r^2=-1 (s)
  void edge_pairing_double(int i); // edge relation double pairing

public:
};


#endif