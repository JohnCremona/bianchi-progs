// FILE SWAN.H: declaration of Swan's algorithm functions

#if     !defined(_SWAN_H)
#define _SWAN_H      1       //flags that this file has been included

#include <iostream>
#include <set>

#include "swan_utils.h"

class SwanData {
public:

  CuspList alist;               // Main list of alphas
  std::set<Quad> a_denoms; // Denominators of alphas
  map<vector<RAT>, int> a_ind;  // Index of an alpha in alist (from coords as key)
  vector<int> a_inv;            // Permutation of order 2 swapping a to a' where M_a(oo)=a'
  vector<int> a_flip;           // Permutation of order 2 swapping alpha to -alpha mod 1
  vector<mat22> Mlist;          // List of M_a with det(M_a)=1 such that M_a(a)=oo and M_a(oo) in alphas
  vector<int> edge_pairs_plus;  // indices of first of a pair (r/s, -r/s) with r^2=+1 (mod s)
  vector<int> edge_pairs_minus; // indices of first of a pair (r/s, -r/s) with r^2=-1 (mod s)
  vector<int> edge_fours;       // indices of first of a 4-tuple (r1,-r1,r2,-r2) of alphas with r1*r2=-1 (mod s)
  vector<vector<Quad>> alpha_sets; // list of all {s, r1, r2} with r1*r2=-1(s)

  CuspList slist;
  map<vector<RAT>, int> s_ind; // Index of a sigma in slist (from coords as key)
  vector<int> s_flip;          // Permutation of order 2 swapping sigma to -sigma mod 1

  H3pointList corners;
  vector<POLYHEDRON> singular_polyhedra, principal_polyhedra, all_polyhedra;

  // all oriented faces up to GL2-equivalence, rotation and reflection, a single mixed list of:
  // - principal triangles T {a1, oo, a2} with a1 reduced fundamental
  // - principal squares   Q {a1, oo, a2, a3} with a1 reduced fundamental
  // - principal hexagons  H {a1, oo, a2, a3, a4, a5, a6} with a1 reduced fundamental
  // - singular triangles  U {a, oo, s} with a reduced fundamental, s singular
  vector<CuspList> all_faces;
  // Separate lists of faces of type T, U, Q, H excluding redundants
  vector<CuspList> aaa, aas, sqs, hexs;
  // Encoded faces as POLYGONS of type T, U, Q, H excluding redundants
  vector<POLYGON> T_faces, U_faces, Q_faces, H_faces;

  vector<vector<int>> M32;  // matrix (encoded as vector<vector<int>>) with one row per polyhedron
                            // giving its boundary as a Z-linear combination of oriented faces
                            // (created by make_all_faces())

  SwanData(int s=0) :showtimes(s), maxn(0) {;} // constructor: does no work

  // clears all alphas-sigma data, polyhedra, faces
  void clear();
  // clears all and recomputes alphas-sigma data, polyhedra, faces
  void create(int verbose=0);
  // read from geodata file, return 1 if successful or 0 if file absent
  int read(string subdir="", int verbose=0)
  {
    clear();
    return read_geodata(0, subdir, verbose>1);
  }
  // read from geodata file, or create from scratch and store if not successful (file absent)
  void read_or_create(string subdir="", int verbose=0);

  void make_sigmas();
  CuspList get_sigmas() {
    make_sigmas();
    return slist;
  }
  int n_sig(int finite_only=0) const {
    int n = slist.size();
    return finite_only? n-1: n;
  }
  // Return i such that slist[i]=s, else -1
  int sigma_index(const RatQuad& s) const;
  // Return i and set shift such that slist[i]+shift=a/b, else -1
  int sigma_index_with_translation(const Quad& a, const Quad& b, Quad& shift) const;
  // Return i and set shift such that slist[i]+shift=s, else -1
  int sigma_index_with_translation(const RatQuad& s, Quad& shift) const
  {return sigma_index_with_translation(s.num(), s.den(), shift);}
  // Return i such that slist[i]+shift=s, else -1
  // (same as above for when shift is not needed)
  int sigma_index_upto_translation(const Quad& a, const Quad& b) const
  {return sigma_index_upto_translation(RatQuad(a,b));}
  // Return i such that slist[i]+shift=s, else -1
  // (same as above for when shift is not needed)
  int sigma_index_upto_translation(const RatQuad& s) const;

  void make_alphas(int verbose=0);
  CuspList get_alphas(int verbose=0) {
    make_alphas(verbose);
    return alist;
  }
  int n_alph() const {return alist.size();}
  // Return i such that alist[i]=a, else -1
  int alpha_index(const RatQuad& a) const;
  // Return i and set t such that alist[i]+t=a, else -1
  int alpha_index_with_translation(const RatQuad& a, Quad& t) const;
  // Return i such that alist[i]+t=a, else -1
  // (same as above for when t is not needed)
  int alpha_index_upto_translation(const RatQuad& a) const;

  RatQuad base_point(int t) const {
    return t>=0? alist[t] : slist[-t];
  }

  H3pointList get_corners() const {
    return corners;
  }

  void make_principal_polyhedra(int verbose=0);
  void make_singular_polyhedra(int verbose=0);
  void make_all_polyhedra(int verbose=0);

  void make_all_faces(int verbose=0);
  // Encode all faces found as POLYGONs in T_faces, U_faces, H_faces,
  // Q_faces; report if verbose; check their encodings/decodings for
  // consistency if check
  int encode_all_faces(int check=1, int verbose=0);

  // For use after reading face data from a geodata file.  From
  // T_faces, U_fces, Q_faces, H_faces reconstruct all_faces.  Must
  // include the standard faces not (for historical rasons) included
  // in geodata files, namely: (1) the universal triangle {0,oo,1};
  // (2) one standard triangle for d=19, 43, 67, 163; (3) one extra Q
  // or H face for d=2, 7, 11.
  void decode_all_faces();

  // Write alpha data to subdir/geodata_d.dat (overwriting if
  // existing) in A-lines.  If include_small_denoms, output all alpha
  // orbits, otherwise omit those of denominator 1, 2, 3.
  void output_alphas(int include_small_denoms=0, string subdir="");
  // Append sigma data to subdir/geodata_d.dat in S-lines.  If
  // include_small_denoms, output all finite sigma orbits (up to
  // sign), otherwise omit oo and those of denominator 2, 3.
  void output_sigmas(int include_small_denoms=0, string subdir="");
  // Write alpha and sigma data to subdir/geodata_d.dat (overwriting
  // if existing) in A- and S-lines. If include_small_denoms, output
  // all alpha and and sigma orbits , otherwise omit those of
  // denominator 0,1,2,3.
  void output_alphas_and_sigmas(int include_small_denoms=0, string subdir="")
  {
    output_alphas(include_small_denoms, subdir);
    output_sigmas(include_small_denoms, subdir);
  }

  // Read from subdir/geodata_d.dat all A- and S-lines and fill alist,
  // slist, Mlist, a_denoms, a_ind, a_inv, a_flip, edge_pairs_plus,
  // edge_pairs_minus, edge_fours.  If include_small_denoms, the file
  // contains all alpha and sigma orbits, otherwise construct those of
  // denominator 0, 1, 2, 3 on the fly.

  // Return 1 if the geodata file exists, else 0.
  int read_alphas_and_sigmas(int include_small_denoms=0, string subdir="", int verbose=0);

  // Append to subdir/geodata_d.dat all TUQH lines from
  // T_faces, U_faces, Q_faces, H_faces
  void output_face_data(string subdir="", int verbose=0);

  //
  void output_geodata(string subdir="", int include_small_denoms=0, int verbose=0)
  {
    output_alphas_and_sigmas(include_small_denoms, subdir);
    output_face_data(subdir, verbose);
  }

  // Read from subdir/geodata_d.dat all TUQH lines and fill T_faces,
  // U_faces, Q_faces, H_faces.  Return 1 if Euclidean, or file exists
  // and has at least 1 TUQH line.
  int read_face_data(string subdir="", int verbose=0);

  // Read from subdir/geodata_d.dat all ASTUQH lines and fill
  // everything needed for homspace computations.  Return 2 if geodata
  // file exists with all data, 1 if it exists with alpha-sigma data
  // only, or 0 if it does not exist.

  // (1) alist, slist, Mlist, a_denoms, a_ind, a_inv, a_flip,
  // edge_pairs_plus, edge_pairs_minus, edge_fours.  If
  // include_small_denoms, the file contains all alpha and sigma
  // orbits, otherwise construct those of denominator 0, 1, 2, 3 on
  // the fly.

  // (2) T_faces, U_faces, Q_faces, H_faces

  int read_geodata(int include_small_denoms=0, string subdir="", int verbose=0)
  {
    int res = read_alphas_and_sigmas(include_small_denoms, subdir, verbose);
    if (!res) return 0;
    res = read_face_data(subdir, verbose);
    return 1+res;
  }

  // return the invariants of H_1 as a Z-module for either GL2
  // (group=1) or SL2 (group=2) or both (group=3)
  vector<vector<long>> integral_homology(int group, int debug=0);

private:
  timer SwanTimer;
  int showtimes;
  CuspList alistx; // list of alphas + 8 integer translates
  CuspList alistF4; // sublist of those in quarter rectangle
  INT maxn; // max denom norm of alphas considered systematically
  std::set<RatQuad> alist_ok, alist_open; // partition of current alphas (ok=surrounded, open=not yet)
  map<RatQuad, CuspList> nbrs, nbrs_ok, nbrs_open;
  CuspList slistx; // list of sigmas + 8 integer translates

  // add one alpha; use covered=1 after finding covering alphas and saturating with more
  // (called by add_new_alphas() and by saturate_aphas())
  int add_one_alpha(const RatQuad& a, int covered=0, int verbose=0);

  // add next batch of alphas from denom_looper, return number added
  // (called only by find_covering_alphas())
  int add_new_alphas(Quadlooper& denom_looper, int verbose=0);

  // (called only by make_alphas())
  void find_covering_alphas(int verbose=0);
  void saturate_alphas(int verbose=0);

  // list of singular corners [s,0] on S_a (s in slist or a translate)
  // (called only by find_corners_from_one())
  H3pointList singular_corners(const RatQuad& a);

  // Find and fill corners list, replacing alist/alistF4 with sublist of alphas (/in F4) on >=3 corners
  // (called only by saturate_alphas())
  void find_corners(int verbose=0);
#if(0)
  void old_find_corners(int verbose=0);
#endif
  // After an unsuccessful saturation loop which produces extra alphas
  // properly covering some old corners, use these to compute more corners
  // (called only by saturate_alphas())
  H3pointList find_extra_corners(const CuspList& extra_alphas);
  // Find corners from one alpha.  The new corners are not added to
  // the class list points, but are returned.  The parameter redundant
  // is set to 1 if a has <3 corners (including singular ones).
  // (called only by find_corners() and find_extra_corners())
  H3pointList find_corners_from_one(const RatQuad& a, int& redundant, int verbose=0);

  // (called only by are_alphas_surrounded())
  int are_sigmas_surrounded(int verbose=0);
  // (called only by are_sigmas_surrounded())
  int is_sigma_surrounded(const RatQuad& s, int verbose=0);

  // (called only by find_covering_alphas())
  int are_alphas_surrounded(int verbose=0);
  // (called only by are_alphas_surrounded())
  int is_alpha_redundant(const RatQuad& a, int verbose=0);
  // (called only by are_alphas_surrounded())
  int is_alpha_surrounded(const RatQuad& a, int verbose=0);
  // (called only by are_alphas_surrounded())
  int are_intersection_points_covered(const RatQuad& a, const RatQuad& b, int verbose=0);
  // test if a is singular by reducing to rectangle and comparing with slist
  // (called only by are_intersection_points_covered())
  int is_singular(const RatQuad& a);

  void process_sigma_orbit(const Quad& r, const Quad& s);
  void process_sigma_orbit(const vector<Quad>& rs) {
    process_sigma_orbit(rs[0], rs[1]);
  }
  void process_sigma_orbit(const RatQuad& sig) {
    process_sigma_orbit(sig.num(), sig.den());
  }

  void process_alpha_orbit(const Quad& s, const Quad& r1, const Quad& r2, int verbose=0);
  void process_alpha_orbit(const vector<Quad>& sr1r2, int verbose=0) {
    process_alpha_orbit(sr1r2[0], sr1r2[1], sr1r2[2], verbose);
  }

  void make_alpha_orbits();
  POLYHEDRON make_principal_polyhedron(const H3point& P, std::set<int>& orbit, int verbose=0);

  // Return the image under delta of the face, as a vector of length #alist+#slist-1
  vector<int> face_boundary_vector(const CuspList& face);

  // Return the index of an edge {a,b} in the range 0..#alphas+#sigmas-2
  // For e={a,b} with a,b principal return i (>=0), where [oo,alphas[i]]
  // is SL(2,Ok)-equivalent to e.

  // If a is principal and b singular, return i+len(alphas)-1 (>=1),
  // where [oo,sigmas[i]] is SL(2,Ok)-equivalent to e.

  // If a is singular and b principal, return -(i+len(alphas)-1) (<0),
  // where [sigmas[i],oo] is SL(2,Ok)-equivalent to e.

  // Raise an error if both a and b are singular.

  // Note that we assume sigmas[0]=oo which is not a singular point;
  // the number of singular points is len(sigmas)-1.

  int edge_index(const EDGE& e);

  // Return the edge boundary matrix M10 (matrix of delta: 1-chains -> 0-chains).
  // The edge basis consists of first the [a,oo] for a in alist (whose
  // boundary is trivial), then the [s,oo] for s in slist[1:].  The
  // cusp basis is indexed by ideal classes. Same for SL2 and GL2.
  vector<vector<int>> edge_boundary_matrix();

  // Use alist, slist, edge_pairs_plus, edge_pairs_minus, fours to
  // return a matrix with one row per pair of glued oriented edges:
  vector<vector<int>> edge_pairings(int GL2, int debug=0);
  // Use all_faces to return a matrix with one row per face giving its
  // boundary as a Z-linear combination of oriented edges
  vector<vector<int>> face_boundaries(int GL2, int debug=0);
  // Row concatenation of previous 2, giving the matrix M21 of the
  // boundary map from 2-cells to 1-cells
  vector<vector<int>> face_boundary_matrix(int GL2, int debug=0);

  // functions for pseudo-euclidean algorithm

  // return type t if translation by -q followed by Mlist[t] reduces
  // (a,b) to (a',b') with N(b')<N(b), or -1 if none exists, which
  // will be if and only if a/b is singular.  If lucky is nonzero,
  // return the first t which reduces N(b), otherwise try all and
  // return the best.

  int find_best_alpha(const Quad& a, const Quad& b, Quad& shift, int lucky=0) const;

/******************************************************************

One pseudo-Euclidean step applies a translation and *if possible* an
M_alpha inversion to a/b (or column vector [a;b]) reducing b, also
multiplying row vector [c,d] my M_alpha on the right.  In the
Euclidean case, the shift is -q where q=a/b (rounded) and the
inversion is via S=[0,-1;1,0], with t=0.  In general if t>=0 then
the t'th inversion was applied.

If the class number is >1 and the ideal (a,b) is non-principal, then
possibly after translation we have that a/b is a singular point s, in
which case no inversion is done (as none can reduce N(b)); then t<0
where s is the |t|'th singular point.  (The singular points list
effectively starts at index 1.)

a,b,c1,d1,c2,d2 are changed in place, though if either (c2,d2) or
both (c1,d1), (c2,d2) are left as defaults they are not updated.

When applied repeatedly, there are two possible stopping
conditions; note that the ideal (a,b) is unchanged throughout since
we only apply SL(2,O_K)-transformations.  The process is guaranteed
to stop after a finite number of steps since either N(b) is reduced
or the second stopping conditions is reached.

(1) When the ideal (a0,b0) is principal, the stopping condition is
b==0.  Then a = (a0,b0) = a0*d1-b0*d2, and c2/c1=a0/b0 reduced to
lowest terms.

(2) When (a0,b0) is not principal, the stopping condiction is
t<0. Then a/b is the |t|'th singular point, represented as a
fraction with ideal (a,b)=(a0,b0), which may not be the "standard"
representation of the singular point r/s.  Since a/b=r/s, we have
a/r=b/s=lambda, say, where lambda*r=a and lambda*b=s, but *lambda
is not integral* in general.

The quantities c1*a+d1*b andc2*a+d2*b are invariant, i.e. M*v is
invariant where M=[c1,d1;c2,d2] and v=[a;b], since we multuply v on
the left by unimodular matrices at the same time as multiplyin M on
the right by the inverse.

We could use a mat22 [c2,d2;c1,d1] (note the numbering) instead of
c1,d1,c2,d2, but we do not always need one or both pairs: for a gcd
computation, we need neither, for a bezout we need c1,d1 and for
continued fractions we need both.

 This function is crucial in reducing modular symbols, expressing
 each as a linear combination of generalised M-symbols.
*********************************************************************/
public:
  void pseudo_euclidean_step(Quad& a, Quad& b, int& t,
                             Quad& c1=Quad::zero, Quad& d1=Quad::zero,
                             Quad& c2=Quad::zero, Quad& d2=Quad::zero) const;

  // Generalization of extended Euclidean algorithm.

  // Given a,b, returns M=[d1,-d2;-c1,c2] (the inverse of [c2,d2;c1,d1])
  // such that g/h = M(a/b), i.e. g=d1*a-d2*b, h = -c1*a+c2*b, where
  //
  // (1) if (a,b) is principal then (a,b)=(g) and h=0, and s=0;
  // (2) otherwise (a,b)=(g,h) and g/h is  the s'th singular point (s>=1).
  //
  // So in Case (1) we get essentially the same information as
  // quadbezout_psea(aa, bb, xx, yy) with xx,yy the top row of M.
  //
  // Note that in case 2, g/h is equal to the singular point sigma=g0/h0
  // as an element of the field, but not as a fraction, since the ideal
  // (g,h)=(a,b) is in the same (non-principal) ideal class as
  // (g0,h0). but is not the same ideal.  In fact, (g,h)=lambda*(g0,h0)
  // with lambda = g/g0 = h/h0 (since g*h0=h*g0), but in general lambda
  // is not integral.

  mat22 generalised_extended_euclid(const Quad& aa, const Quad& bb, int& s) const;

  // Each relation is a signed sum of edges (M)_alpha = {M(alpha},
  // M(oo)} for M in the list mats and alpha=alphas[t] (when t>=0) or
  // sigmas[-t] (when t<0), for t in the list types.  Here we check that such a
  // relation holds identically in H_3 (not just modulo the congruence
  // subgroup!)

  // General case:
  int check_rel(const vector<mat22>& mats, const vector<int>& types, const vector<int>& signs) const;
  // Special case: all signs +1
  int check_rel(const vector<mat22>& mats, const vector<int>& types) const;

  int check_aaa_triangle(const POLYGON& T, int verbose=0) const;
  int check_aas_triangle(const POLYGON& T, int verbose=0) const;
  int check_triangles(int verbose=0) const;
  int check_square(const POLYGON& squ, int verbose=0) const;
  int check_squares(int verbose=0) const;
  int check_hexagon(const POLYGON& hex, int verbose=0) const;
  int check_hexagons(int verbose=0) const;
  int check_faces(int verbose=0) const
  {
    // We want to check all even if some fail
    int TU_ok = check_triangles(verbose);
    int Q_ok = check_squares(verbose);
    int H_ok = check_hexagons(verbose);
    return TU_ok && Q_ok && H_ok;
  }
};

#endif
