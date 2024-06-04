// FILE EDGE_RELATIONS.H: Declaration of class edge_relations

#if     !defined(_EDGE_RELATIONS_H)
#define _EDGE_RELATIONS_H      1       //flags that this file has been included

#include "geometry.h"
#include "P1N.h"
#include "swan.h"

class face_relations;

class edge_relations {
public:
  edge_relations() {;} //{cout<<"In default edge_relations constructor"<<endl;}
  edge_relations(P1N*, int plus, int verb=0, long ch=0);
  long coords(int i) const {return coordindex[i];}
  long coords(int i, int t) const {return coordindex[i+offset(t)];}
  long gen(int i) const {return gens[i];}  // indexed from 1 not 0
  long get_ngens() const {return ngens;}
  long offset(int t) const // offset into coordindex for type t
  {
    return nsymb * (t>=0? t: Quad::SD.n_alph()-t-1);
  }
  pair<long, int> symbol_number_and_type(long i) const
  {
    std::ldiv_t st = ldiv(i, nsymb);
    long s = st.rem;  // remainder gives (c:d) symbol number
    int t = st.quot;  // quotient gives symbol type
    if (t >= Quad::SD.n_alph())  // convert singular type to negative
      t = Quad::SD.n_alph()-t-1;
    return {s, t};
  }
  action act_with(const mat22& M);
  action act_with(const Quad& a, const Quad& b, const Quad& c, const Quad& d);

protected:
  P1N* P1; // provides nsymb, symbol(i), symbops
  int plusflag, verbose;
  long characteristic; // =0 or prime p
  int n_alphas, n_sigmas;
  long nsymb, ngens;
  vector<int> coordindex, gens;
  mat22 M_alpha(int j) const {return Quad::SD.Mlist[j];}

private:
  void edge_relations_1();    // basic edge relations for alpha = 0
  void edge_relations_2();    // edge relations for alphas &sigmas with denom 2
  // previous calls one of the following three:
  void edge_relations_2_d3mod8(); // 2 inert, 2 alphas, 0 sigmas
  void edge_relations_2_d7mod8(); // 2 split, 0 alphas, 2 sigmas
  void edge_relations_2_d12mod4(); // 2 ramified, alpha=w/2, sigma=(w+1)/2 or vice versa
  void sigma_relations();          // for sigma with 2*sigma not integral
  void edge_pairing_minus(int i);   // edge relation pair, alpha=r/s with r^2=-1 (s)
  void edge_pairing_plus(int i);    // edge relation pair, alpha=r/s with r^2=+1 (s)
  void edge_pairing_double(int i);  // edge relation double pairing
  void report();
  friend class face_relations;
};

#endif
