// FILE SWAN_ALPHAS.H: declaration of alpha-finding functions for Swan's algorithm

#if     !defined(_SWAN_ALPHAS_H)
#define _SWAN_ALPHAS_H      1       //flags that this file has been included

#include "swan_utils.h"

// Return sorted list of (saturated covering) alphas (denom 1, denom 2, denom 3, larger denoms in pairs or fours)
// + pluspairs, minuspairs, fours
CuspList sort_alphas(const CuspList& A,
                     vector<vector<Quad>>& pluspairs, vector<vector<Quad>>& minuspairs, vector<vector<Quad>>& fours,
                     int verbose=0, int debug=0);

// Return sorted list of (saturated covering) alphas (denom 1, denom
// 2, denom 3, then larger denoms in orbit pairs or fours) and fill triples
// with {s,r1,r2} where r1*r2=-1 (mod s)
CuspList alpha_orbits(const CuspList& alist,
                      vector<vector<Quad>>& triples,
                      int verbose=0, int debug=0);

// Output sorted list of alphas (denom > 3 in pairs or fours)
void output_alphas(const vector<vector<Quad>>& pluspairs,
                   const vector<vector<Quad>>& minuspairs,
                   const vector<vector<Quad>>& fours,
                   int to_file=1, int to_screen=0);

// direct lists of alphas of denominator 2 or 3:
CuspList denom_2_alphas();
CuspList denom_3_alphas();

// Given principal cusps a1, a2, a such that the circles S_a1 and
// S_a2 intersect in distinct points, test whether S_a covers either
// or both these points.

// Returns 0 if neither, 2 if both, +1 or -1 if just one.  The signs
// are consistent so that if a returns +1 and a' returns -1 then each
// intersection point is covered by either S_a or S_a'.

int are_intersection_points_covered_by_one(const RatQuad& a1, const RatQuad& a2, const RatQuad& a,
                                           int debug=0);

// Given principal cusps a0, a1 whose circles S_a0, S_a1 intersect,
// and a list of principal cusps alist each of whose circles S_a also
// intersects S_a0, test whether each of the two intersection points
// of S_a0 and S_a1 is either singular or strictly inside one of the
// S_a.  We treat as a special case when the two intersection points
// are in k.  If not, the code still uses exact arithmetic.

int are_intersection_points_covered(const RatQuad& a0, const RatQuad& a1, const CuspList& alist,
                                    const CuspList& slist, int debug=0);


// Given a principal cusp a0, a candidate list of principal cusps
// alist, tests whether the boundary of the disc S_a0 is contained in
// the union of the S_b for b in alist, apart from any singular
// points on the boundary.  It suffices to consider all b such that
// S_b intersects S_a in two points and check that each of the
// points is either singular or contained in some other S_b.  This
// is simplest when the intersection points are in k; if not then the
// method still uses exact arithmetic in k throughout.

// pairs_ok is a list of {a1,a2} whose intersection points are
// known to be covered, so can be skipped.

// Returns 0 or 1, and pairs_ok is updated to a list of pairs
// whose intersections have now been shown to be covered.

int is_alpha_surrounded(const RatQuad& a0, const CuspList& alist, const CuspList& slist,
                        vector<CuspPair>& pairs_ok, int debug=0);


// Given alist_ok and alist_open, lists of principal cusps, and slist,
// a complete list of singular points, tests whether the boundary of
// every disc S_a for a in alist_open is contained in the union of
// the translates of the S_b for b in alist_ok+alist_open, apart from
// any singular points on the boundary.

// Any a which pass are added to alist_ok, and removed from
// alist_open, so success means that the latter is empty.  This allows
// for incremental testing by adding more a to alist_open.

// pairs_ok is list of pairs {a1,a2} whose intersection points are
// known to be covered, which will be added to.

// Returns 0/1

// NB All a in alist_open will be tested, i.e. we carry on after a
// failure.

int are_alphas_surrounded(CuspList& alist_ok, CuspList& alist_open,
                          const CuspList& slist, vector<CuspPair>& pairs_ok,
                          int verbose=0, int debug=0);
// Returns a finite list of principal cusps a such that the S_{a+t}
// for all integral t cover CC apart from singular points.

// For n>=1 successively, we test as a candidate set all a=r/s with
// r,s coprime, r reduced mod s, N(s)<=n (omitting any for which S_a
// is contained in any earlier S_a') until we succeed.

// slist can be set to a list of singular points (up to translation),
// otherwise these will be computed.

// Other functions will then (1) saturate the set, (2) discard redundancies.
CuspList covering_alphas(const CuspList& slist, int verbose=0);

// Given a list of principal cusps alpha (all reduced mod O_k) return
// a list of "corners" P = [z,tsq] each the intersection of an S_a
// with at least two other S_{b+t} with z in the fundamental
// rectangle and tsq>0.
H3pointList triple_intersections(const CuspList& alist, int debug=0);

H3pointList new_triple_intersections(const CuspList& alist, int debug=0);
H3pointList old_triple_intersections(const CuspList& alist, int debug=0);

// count how many P in points are on S_a
int nverts(const RatQuad& a, const H3pointList& points);

// return sublist of a in alist which have at least 3 vertices in points
CuspList remove_redundants(const CuspList& alist, const H3pointList& points);

// Given a list of covering alphas with max denom norm maxn, return a saturated irredundant list
CuspList saturate_covering_alphas(const CuspList& alist, const CuspList& slist, INT maxn, int debug=0, int verbose=0);

// return  a saturated irredundant list of alphas in the fundamental rectangle
CuspList find_alphas(const CuspList& slist, int debug=0, int verbose=0);

// return list of alphas (or translates) which pass through a finite cusp
CuspList neighbours(const RatQuad& a, const CuspList& alist);

// return sorted list of alphas (or translates) which pass through a finite cusp,
// i.e. angle_under_pi(sigma, a[i-1], a[i]) for all 0<=i<n and angle_under_pi(sigma, a[n-1], a[0]).
CuspList sorted_neighbours(const RatQuad& sigma, const CuspList& alist);

// test if all singular points (slist) are surrounded by alpha circles:
int are_sigmas_surrounded(const CuspList& slist, const CuspList& alist, int debug=0);
// test if one singular point (sigma) is surrounded by alpha circles:
int is_sigma_surrounded(const RatQuad& sigma, const CuspList& alist, int debug=0);

// Output string for one alpha orbit
// Returns "A sr si r1r r1i r2r r2i" for alphas r1/s, r2/s with r1*r2=-1 (mod s)
string make_A_line(const Quad& s, const Quad& r1, const Quad& r2);
// same with parameter {s,r1,r2};
inline string make_A_line(const vector<Quad>& sr1r2)
{
  return make_A_line(sr1r2[0], sr1r2[1], sr1r2[2]);
}
#endif
