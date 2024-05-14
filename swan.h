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

  vector<vector<int>> M32;  // matrix (encoded as vector<vector<int>>) with one row per polyhedron
                            // giving its boundary as a Z-linear combination of oriented faces

  SwanData(int s=0) :showtimes(s), maxn(0) {;}

  void make_sigmas();
  CuspList get_sigmas() {
    make_sigmas();
    return slist;
  }
  void output_sigmas(int include_small_denoms=0, string subdir="");

  void make_alphas(int verbose=0);
  CuspList get_alphas(int verbose=0) {
    make_alphas(verbose);
    return alist;
  }
  void output_alphas(int include_small_denoms=0, string subdir="");
  void output_alphas_and_sigmas(int include_small_denoms=0, string subdir="")
  {
    output_alphas(include_small_denoms, subdir);
    output_sigmas(include_small_denoms, subdir);
  }
  void read_alphas_and_sigmas(int include_small_denoms=0, string subdir="");

  H3pointList get_corners() const {
    return corners;
  }

  void make_principal_polyhedra(int verbose=0);
  void make_singular_polyhedra(int verbose=0);
  void make_all_polyhedra(int verbose=0);

  void make_all_faces(int verbose=0);
  // Report on faces found if verbose; check their encodings/decodings for consistency:
  int check_all_faces(int verbose=0);

  void output_face_data(string subdir="", int verbose=0);

private:
  timer SwanTimer;
  int showtimes;
  Quadlooper denom_looper; // default init: norms from 1 to oo, both conjugates
  CuspList alistx; // list of alphas + 8 integer translates
  CuspList alistF4; // sublist of those in quarter rectangle
  INT maxn; // max denom norm of alphas considered systematically
  CuspList alist_ok, alist_open; // partition of current alphas (ok=surrounded, open=not yet)
  map<RatQuad, CuspList> nbrs, nbrs_ok, nbrs_open;
  CuspList slistx; // list of sigmas + 8 integer translates

  // add one alpha; use covered=1 after finding covering alphas and saturating with more
  // (called by add_new_alphas() and by saturate_aphas())
  int add_one_alpha(const RatQuad& a, int covered=0, int verbose=0);

  // add next batch of alphas from denom_looper, return number added
  // (called only by find_covering_alphas())
  int add_new_alphas(int verbose=0);

  // (called only by make_alphas())
  void find_covering_alphas(int verbose=0);
  void saturate_alphas(int verbose=0);

  // list of singular corners [s,0] on S_a (s in slist or a translate)
  // (called only by find_corners_from_one())
  H3pointList singular_corners(const RatQuad& a);

  // Find and fill corners list, replacing alist/alistF4 with sublist of alphas (/in F4) on >=3 corners
  // (called only by saturate_alphas())
  void find_corners(int verbose=0);
  void old_find_corners(int verbose=0);

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

  void process_alpha_orbit(const Quad& s, const Quad& r1, const Quad& r2);
  void process_alpha_orbit(const vector<Quad>& sr1r2) {
    process_alpha_orbit(sr1r2[0], sr1r2[1], sr1r2[2]);
  }

  void make_alpha_orbits();
  POLYHEDRON make_principal_polyhedron(const H3point& P, std::set<int>& orbit, int verbose=0);
};

#endif
