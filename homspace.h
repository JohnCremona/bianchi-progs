// FILE HOMSPACE.H: Declaration of class homspace

#if     !defined(_HOMSPACE_H)
#define _HOMSPACE_H      1       //flags that this file has been included

#include "cusp.h"
#include "face_relations.h"
#include "hecke.h"

class homspace {
friend class newforms;
public:
  int verbose;
  scalar denom1, denom2, denom3;
  int dimension, cuspidal_dimension, ncusps;
  edge_relations ER;
  face_relations FR;
  int plusflag;
  Qideal N; // the level
  P1N P1;
  int ngens, nsymb, nap, nwq;

  ssubspace kern;
  mat deltamat;
  smat tkernbas, bigtkernbas; // transpose of kernel(delta) basis in terms of freegens, gens
  vector<modsym> freemods;
  vector<int> freegens;
  mat projcoord;
  scalar characteristic; // =0 or prime p
  scalar modulus; // prime for linear algebra (=characteristic if that is >0)
  scalar hmod; // if >0, did not lift from modular linear algebra so coord is modulo this

  homspace(const Qideal& I, scalar mod, int hp, int verb=0, scalar ch=scalar(0));

  void kernel_delta();          // computes ker(delta) for cuspidal homology
  void make_freemods();         // computes freemods and needed
  modsym edge_generator(int i); // for i in 0..nsymb*(#types)-1, the i'th edge as a modsym
  modsym generator(int j)      // for j in 1..ngens, the j'th generating modsym
  {return edge_generator(ER.gen(j));}

  int check_conjugate(int verb=0);     // function to check consistency between this and conjugate level

  vec coords(int i) {return FR.coords(i);}
  int coords(const Quad& c, const Quad& d);
  int index(const Quad& c, const Quad& d) {return P1.index(c, d);}

  vec cuspidalpart(const vec& v) const {return v[pivots(kern)];}
  int is_cuspidal(const subspace& s) const; // test for cupidality (of a non-dual subspace only)

  int h1cuspdim() const {return cuspidal_dimension;}
  int h1dim() const {return dimension;}  // No confusion with subspace::dim
  int h1ncusps() const {return ncusps;}
  scalar h1denom() const {return denom1;}
  scalar h1cdenom() const {return denom3;}
  scalar h1hmod() const {return hmod;}

  vec chaincd(const Quad& c, const Quad& d, int type=0, int proj=0);
  vec chain(const Quad& a, const Quad& b, int proj=0, const Quad& c=Quad::zero, const Quad& d=Quad::one);
  vec chain(const RatQuad& r, int proj=0, const Quad& c=Quad::zero, const Quad& d=Quad::one)
  {return chain(r.num(), r.den(), proj, c, d);}
  vec chain(const RatQuad& alpha, const RatQuad& beta, int proj=0);
  vec chain(const modsym& m, int proj=0)
  {return chain(m.alpha(), m.beta(), proj);}

  vec applyop(const matop& T, const RatQuad& m, int proj=0);
  vec applyop(const matop& T, const modsym& m, int proj=0);

  mat calcop(const matop& T, int cuspidal=0, int dual=1, int display=0);
  ZZX charpoly(const matop& T, int cuspidal=0);
  mat calcop(const gmatop& T, int cuspidal=0, int dual=1, int display=0);
  ZZX charpoly(const gmatop& T, int cuspidal=0);
  vec calcop_col(const matop& T, int j, int)  {return applyop(T,freemods[j-1]);}
  mat calcop_cols(const matop& T, const vec_i& jlist, int);
  mat calcop_restricted(const matop& T, const subspace& s, int dual, int display);
  smat s_calcop(const matop& T, int cuspidal, int dual, int display);
  smat s_calcop_cols(const matop& T, const vec_i& jlist, int verb=0);
  smat s_calcop_restricted(const matop& T, const ssubspace& s, int dual, int display);

  vec maninvector(Quadprime& P, int proj=0);
  vec manintwist(const Quad& lambda, const vector<Quad>& res, vector<int> chitable, int proj=0);
  // no longer used, implementation commented out:
  // vec newhecke(const Quadprime& P, const Quad& n, const Quad& d);
  // The (dual) subspace cut out by the given eigenvalues for the
  // basis of Quad::class_group_2_rank unramified characters.
  ssubspace unramified_character_subspace(const vector<int>& eigs);

  // The dimension and cuspidal of the previous space (for when we do
  // not need the subspace itself).
  pair<int,int> unramified_character_subspace_dimensions(const vector<int>& eigs);

  // The dimension or cuspidal dimension of the previous space (for
  // when we do not need the subspace itself).
  int unramified_character_subspace_dimension(const vector<int>& eigs, int cuspidal)
  {
    pair<int,int> dims = unramified_character_subspace_dimensions(eigs);
    return (cuspidal? dims.second: dims.first);
  }

  // The (dual) subspace with eigenvalue +1 for all quadratic
  // unramified characters.
  ssubspace trivial_character_subspace()
  {
    return unramified_character_subspace(vector<int>(Quad::class_group_2_rank, +1));
  }

  // total and cuspidal dimensions of subspace on which all T(A,A) act trivially
  pair<int,int> trivial_character_subspace_dimensions()
  {
    return unramified_character_subspace_dimensions(vector<int>(Quad::class_group_2_rank, +1));
  }

  // total or cuspidal dimension of subspace on which all T(A,A) act trivially
  int trivial_character_subspace_dimension(int cuspidal)
  {
    pair<int,int> dims = trivial_character_subspace_dimensions();
    return (cuspidal? dims.second: dims.first);
  }

  // list of (total,cuspidal) dimensions of subspaces on which all T(A,A)
  // act trivially with self-twist by unramified quadratic char D for
  // each D (including D=1, meaning no self-twist)
  vector<pair<int,int>> trivial_character_subspace_dimensions_by_twist(int use_lower_bounds, int use_cuspidal_lower_bounds, vector<int> lower_bounds, vector<int> cuspidal_lower_bounds);

  // list of total or cuspidal dimensions of subspaces on which all T(A,A)
  // act trivially with self-twist by unramified quadratic char D for
  // each D (including D=1, meaning no self-twist)
  vector<int> trivial_character_subspace_dimensions_by_twist(int cuspidal, int use_lower_bounds, vector<int> lower_bounds={});
};

// functions for caching homspaces, full Hecke matrices, full and new Hecke polynomials

// Key is ideal_label(N)
extern map<string,homspace*> H1_dict;
homspace* get_homspace(const Qideal& N, scalar mod);

// Key is ideal_label(N)-T.name()
extern map<string, ZZX> full_poly_dict;
ZZX get_full_poly(const Qideal& N,  const matop& T, const scalar& mod);
ZZX get_full_poly(const Qideal& N,  const gmatop& T, const scalar& mod);

extern map<string, mat> full_mat_dict;
mat get_full_mat(const Qideal& N,  const matop& T, const scalar& mod);
mat get_full_mat(const Qideal& N,  const gmatop& T, const scalar& mod);

extern map<string, ZZX> new_poly_dict;
ZZX get_new_poly(const Qideal& N, const matop& T, const scalar& mod);
ZZX get_new_poly(const Qideal& N, const gmatop& T, const scalar& mod);

// Key is ideal_label(N)-mod-p
extern map<string,homspace*> H1_modp_dict;
homspace* get_homspace_modp(const Qideal& N, scalar p);

// Key is ideal_label(N)-ideal_label(P)
extern map<string, ZZ_pX> full_poly_modp_dict;
ZZ_pX get_full_poly_modp(const Qideal& N,  const Quadprime& P, scalar p);

// Key is ideal_label(N)-ideal_label(P)
extern map<string, ZZ_pX> new_poly_modp_dict;
ZZ_pX get_new_poly_modp(const Qideal& N, const Quadprime& P, scalar p);

#endif
