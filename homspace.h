// FILE HOMSPACE.H: Declaration of class homspace

#if     !defined(_HOMSPACE_H)
#define _HOMSPACE_H      1       //flags that this file has been included

#include <eclib/arith.h>
#include <eclib/method.h>
#include <eclib/subspace.h>
#include <eclib/smatrix.h>
#ifdef USE_SMATS
#include <eclib/smatrix_elim.h>
#endif
#include "moddata.h"
#include "cusp.h"
#include "symb.h"

class homspace :public symbdata {
friend class newforms;
public:
  int verbose;
  int cuspidal;  // if 1 then compute cuspidal homology
  vector<int> coordindex, gens, needed, freegens;
  long rk, denom1, denom2, dimension, denom3, ncusps;
  ssubspace kern;
  smat tkernbas; // transpose of kernel(delta) basis
  vector<modsym> freemods;
  mat coord, projcoord;
  long hmod; // if >0, failed to lift from modular linear algebra
             // so coord is modulo this

#ifdef USE_SMATS
  smat relmat;
#else
  mat relmat;
#endif
  long ngens, numrel, maxnumrel;
  void add_face_rel(const vector<int>& rel, const vector<int>& types);

  homspace(const Quad& n, int hp, int cuspid, int verb);

  // The next several methods are called only in the constructor, but
  // are separted out for clarity and for ease of separating thec ode
  // for different fields.
  void edge_relations();      // computes coordindex, gens
  void edge_relations_1();    // basic edge relations for alpha = 0
  void edge_relations_2();    // extra edge relations for alphas with denom 2
  void edge_pairing(int i);   // edge relation pair, alpha=r/s with r^2=-1 (s)
  void edge_pairing_double(int i); // edge relation double pairing

  void face_relations();    // computes face relations, fills relmat
  void triangle_relation_0();   // triangle relation for all fields
  void triangle_relation_1_3();   // extra triangle relation for fields 1, 3
  void triangle_relation_2();   // extra triangle relation(s) for fields 19+
  void cyclic_triangle_relation(int i); // generic cyclic triangle relation
  void general_triangle_relation(const vector<int>& tri);  // generic triangle relation
  void square_relation_2();   // extra square relation for field 2
  void rectangle_relation_7();   // extra rectangle relation for field 7
  void hexagon_relation_11();   // extra hexagon relation for field 11
  void square_relation_19();   // extra square relation for field 19
  void square_relation_43();   // extra square relations for field 43
  void square_relation_67();   // extra square relations for field 67

  void solve_relations();       // computes kernel of relmat and sets rk, denom1, coord, freegens
  void kernel_delta();          // computes ker(delta) for cuspidal homology
  void make_freemods();         // computes freemods and needed

   mat opmat(int, int dual=1, int verb=0);
   vec opmat_col(int, int j, int verb=0);
   mat opmat_cols(int, const vec& jlist, int verb=0);
   mat opmat_restricted(int i,const subspace& s, int dual, int verb=0);
  // versions returning an smat:
   smat s_opmat(int i,int dual, int verb=0);
   svec s_opmat_col(int i, int j, int verb=0);
   smat s_opmat_cols(int i, const vec& jlist, int verb=0);
   smat s_opmat_restricted(int i, const ssubspace& s, int dual,int verb=0);
   long matdim(void) {return dimension;}
   vector<long> eigrange(long i);

   long h1cuspdim() const {return dim(kern);}
   long h1dim() const {return dimension;}  // No confusion with subspace::dim
   long h1denom() const {return denom1;}
   long h1cdenom() const {return denom3;}
   long h1ncusps() const {return ncusps;}
   long h1hmod() const {return hmod;}

  vec chaincd(const Quad& c, const Quad& d, int type=0, int proj=0) const;
  vec chain(const symb& s, int type=0, int proj=0) const
  {return chaincd(s.cee(), s.dee(), type, proj);}
  vec projchaincd(const Quad& c, const Quad& d, int type=0) const
  {return chaincd(c, d, type, 1);}
  vec chain(const Quad& a, const Quad& b, int proj=0, const Quad& c=Quad::zero, const Quad& d=Quad::one) const;
  vec projchain(const Quad& a, const Quad& b, const Quad& c=Quad::zero, const Quad& d=Quad::one) const
  {return chain(a, b, 1, c, d);}
  vec chain(const RatQuad& r, int proj=0, const Quad& c=Quad::zero, const Quad& d=Quad::one) const
  {return chain(num(r),den(r), proj, c, d);}
  vec chain(const RatQuad& alpha, const RatQuad& beta, int proj=0) const;
  vec chain(const modsym& m, int proj=0) const
  {return chain(m.alpha(), m.beta(), proj);}

  vec kernelpart(const vec& v) const
  {return v[pivots(kern)];}
  vec cycle(const symb& s, int type=0) const
  {return kernelpart(chain(s, type));}
  vec cycle(const Quad& n, const Quad& d) const
  {return kernelpart(chain(n,d));}
  vec cycle(const RatQuad& r) const
  {return kernelpart(chain(num(r),den(r)));}
  vec cycle(const modsym& m) const
  {return cycle(m.beta())-cycle(m.alpha());}

  vec applyop(const matop& mlist, const RatQuad& m, int proj=0) const;
  vec applyop(const matop& mlist, const modsym& m, int proj=0) const;

  mat calcop(const string opname, const Quad& p, const matop& mlist, int dual=1, int display=0) const;
  vec calcop_col(const string opname, const Quad& p, const matop& mlist, int j, int display=0) const;
  mat calcop_cols(const string opname, const Quad& p, const matop& mlist, const vec& jlist, int display=0) const;
  smat s_calcop(const string  opname, const Quad& p, const matop& mlist,
                int dual, int display) const;
  svec s_calcop_col(const string  opname, const Quad& p, const matop& mlist,
                    int j, int display) const;
  smat s_calcop_cols(const string  opname, const Quad& p, const matop& mlist,
                     const vec& jlist, int display) const;
  mat calcop_restricted(const string opname, const Quad& p, const matop& mlist, const subspace& s,
                        int dual, int display) const;
  smat s_calcop_restricted(const string opname, const Quad& p, const matop& mlist, const ssubspace& s,
                           int dual, int display) const;

public:
   mat heckeop(const Quad& p, int dual=1, int display=0) const;
   vec heckeop_col(const Quad& p, int j, int display=0) const;
   mat heckeop_cols(const Quad& p, const vec& jlist, int display=0) const;
   smat s_heckeop(const Quad& p, int dual, int display) const;
   svec s_heckeop_col(const Quad& p, int j, int display) const;
   smat s_heckeop_cols(const Quad& p, const vec& jlist, int display) const;
   mat heckeop_restricted(const Quad& p, const subspace& s, int dual, int display) const;
   smat s_heckeop_restricted(const Quad& p, const ssubspace& s, int dual, int display) const;
   mat wop(const Quad& q, int dual=1, int display=0) const;
   mat fricke(int dual=1, int display=0) const;
//   mat conj(int display=0) const;
   vec maninvector(const Quad& p) const;
   vec projmaninvector(const Quad& p) const;
  vec manintwist(const Quad& lambda, const vector<Quad>& res, vector<int> chitable) const;
   vec newhecke(const Quad& p, const Quad& n, const Quad& d) const;
};

vec reduce_modp(const vec& v, const scalar& p=DEFAULT_MODULUS);
mat reduce_modp(const mat& m, const scalar& p=DEFAULT_MODULUS);

int check_face_rel(const vector<mat22>& mats, const vector<int>& types);


#endif
