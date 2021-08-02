// FILE EDGE_RELATIONS.H: Declaration of class edge_relations

#if     !defined(_EDGE_RELATIONS_H)
#define _EDGE_RELATIONS_H      1       //flags that this file has been included

#include "geometry.h"
#include "P1N.h"

class face_relations;

class edge_relations {
public:
  edge_relations() {;}
  edge_relations(P1N*, int plus, int verb=0);
  int coords(int i) const {return coordindex[i];}
  int gen(int i) const {return gens[i];}  // indexed from 1 not 0
  int get_ngens() const {return ngens;}

protected:
  P1N* P1; // provides nsymb, symbol(i), symbops
  int plusflag, verbose, nsymb, nsymbx, ngens;
  vector<int> coordindex, gens;

private:
  void edge_relations_1();    // basic edge relations for alpha = 0

  void edge_relations_2();    // edge relations for alphas &sigmas with denom 2
  // previous calls one of the following four:
  void edge_relations_2_d3mod8(); // 2 inert, 2 alphas, 0 sigmas
  void edge_relations_2_d7mod8(); // 2 split, 0 alphas, 2 sigmas
  void edge_relations_2_d12mod4(); // 2 ramified, alpha=w/2, sigma=(w+1)/2 or vice versa

  void edge_pairing(int i);   // edge relation pair, alpha=r/s with r^2=-1 (s)
  void edge_pairing_double(int i); // edge relation double pairing
  void report();
  friend class face_relations;
};


#endif
