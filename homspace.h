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

  // Return the (dual) subspace of the full space cut out by the given
  // eigenvalues for the basis of Quad::class_group_2_rank unramified
  // characters.  If dual=0 (only) and cuspidal=1, then return the
  // subspace of the cuspidal subspace similarly cut out.
  ssubspace unramified_character_subspace(const vector<int>& eigs, int cuspidal, int dual);

  // The dimension and cuspidal dimension of the previous space (for
  // when we do not need the subspace itself).
  pair<int,int> unramified_character_subspace_dimensions(const vector<int>& eigs);

  // The dimension or cuspidal dimension of the previous space (for
  // when we do not need the subspace itself).
  int unramified_character_subspace_dimension(const vector<int>& eigs, int cuspidal)
  {
    pair<int,int> dims = unramified_character_subspace_dimensions(eigs);
    return (cuspidal? dims.second: dims.first);
  }

  // The (dual) subspace with eigenvalue +1 for all quadratic
  // unramified characters.  Value cached in triv_char_subspace.
  int triv_char_subdim;      // dimension of trivial char subspace, or -1 if not yet computed
  int triv_char_dual_subdim; // dimension of dual trivial char subspace, or -1 if not yet computed
  int triv_char_cuspidal_subdim; // dimension of cuspidal trivial char subspace, or -1 if not yet computed
  ssubspace triv_char_subspace;
  ssubspace triv_char_dual_subspace;
  ssubspace triv_char_cuspidal_subspace;
  // return triv_char_subspace or triv_char_cuspidal_subspace, after computing if necessary
  // NB cannot have cuspidal=dual=1
  ssubspace trivial_character_subspace(int cuspidal, int dual);

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

// Return true iff T's new poly is squarefree and coprime to its old poly
int test_splitting_operator(const Qideal& N, const gmatop& T, const scalar& mod, int verbose=0);

///////////////////////////////////////////////////////////////////////////////////////////
// functions for caching homspaces, full Hecke matrices, full and new Hecke polynomials
///////////////////////////////////////////////////////////////////////////////////////////

// Utilities for creating keys

string Nkey(Qideal& N);
string NPkey(Qideal& N, Qideal& P);
string NTkey(Qideal& N, const matop& T);
string NTkey(Qideal& N, const gmatop& T);
string Nmodpkey(Qideal& N, const scalar p);
string NPmodpkey(Qideal& N, Quadprime& P, scalar p);

// Dicts (i.e. maps) holding various caches (declared in homspace.cc)
// Keys are strings encoding level, operator label created from above utilities

extern map<string,homspace*> H1_dict;  // Key: ideal_label(N)
extern map<string, mat> full_mat_dict; // Key: ideal_label(N)-T.name(). Value: matrix of T
                                       // on full space
extern map<string, ZZX> poly_dict;          // char polys on full space
extern map<string, ZZX> cuspidal_poly_dict; // char polys of restriction to cuspidal subspace
extern map<string, ZZX> new_poly_dict;      // char polys of restriction to new subspace
extern map<string, ZZX> new_cuspidal_poly_dict; // char polys of restriction to new cuspidal subspace
extern map<string, ZZX> tc_poly_dict;           // char polys on trivial char subspace
extern map<string, ZZX> tc_cuspidal_poly_dict;  // char polys on trivial char cuspidal subspace
extern map<string, ZZX> tc_new_poly_dict;       // char polys on trivial char new subspace
extern map<string, ZZX> tc_new_cuspidal_poly_dict;  // char polys on trivial char new cuspidal subspace

// mod p version of some of the above

extern map<string,homspace*> H1_modp_dict; // Key is ideal_label(N)-mod-p
extern map<string, ZZ_pX> full_poly_modp_dict;
extern map<string, ZZ_pX> new_poly_modp_dict;

// Functions to retrieve a value from one of these dicts given its
// key, computing and storing it if the key is not already there:

// from H1_dict
homspace* get_homspace(const Qideal& N, scalar mod);

// from full_mat_dict
mat get_full_mat(const Qideal& N,  const matop& T, const scalar& mod);
mat get_full_mat(const Qideal& N,  const gmatop& T, const scalar& mod);

// from one of poly_dict, tc_poly_dict, cuspidal_poly_dict, tc_cuspidal_poly_dict
// depending on flags cuspidal & triv_char
ZZX get_poly(const Qideal& N,  const gmatop& T, int cuspidal, int triv_char, const scalar& mod);
inline ZZX get_poly(const Qideal& N,  const matop& T, int cuspidal, int triv_char, const scalar& mod)
{return get_poly(N, gmatop(T), cuspidal, triv_char, mod);}

// from one of new_poly_dict, tc_new_poly_dict, new_cuspidal_poly_dict, tc_new_cuspidal_poly_dict
// depending on flags cuspidal & triv_char

// NB In even class number this may raise an error when newspaces at
// lower levels have self-twist since the effect of this on lowering
// oldspace dimensions is ignored.
ZZX get_new_poly(const Qideal& N, const gmatop& T, int cuspidal, int triv_char, const scalar& mod);
inline ZZX get_new_poly(const Qideal& N, const matop& T, int cuspidal, int triv_char, const scalar& mod)
{return get_new_poly(N, gmatop(T), cuspidal,  triv_char, mod);}

// Functions to output and re-input poly dicts
void output_poly_dict(ostream& os, map<string, ZZX> D);
map<string, ZZX> input_poly_dict(istream& is);

// mod p cached functions for homspace, full and new polys
homspace* get_homspace_modp(const Qideal& N, scalar p);
ZZ_pX get_full_poly_modp(const Qideal& N,  const Quadprime& P, scalar p);
ZZ_pX get_new_poly_modp(const Qideal& N, const Quadprime& P, scalar p);

// Functions to output and re-input poly dicts (mod p version)
void output_poly_dict(ostream& os, map<string, ZZ_pX> D);
map<string, ZZ_pX> input_poly_dict(istream& is,  const ZZ& p);


#endif
