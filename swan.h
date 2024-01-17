// FILE SWAN.H: declaration of Swan's algorithm functions

#if     !defined(_SWAN_H)
#define _SWAN_H      1       //flags that this file has been included

#include <iostream>
#include <set>

#include "ratquads.h"

typedef vector<RatQuad> CuspList;  // may refactor using sets later
typedef std::set<RatQuad, Cusp_comparison> CuspPair;
typedef pair<RatQuad,RAT> H3point;
// Comparison function (based only height, highest first)
struct H3_comparison {
  bool operator() (const H3point& lhs, const H3point& rhs) const
  {return lhs.second > rhs.second;}
};
extern H3_comparison H3_cmp;
typedef vector<H3point> H3pointList;

ostream& operator<<(ostream& s, const H3point& P);

// Given an ideal I, return a list of singular points of class [I]
// (one representative for each orbit under integral translations).
CuspList singular_points_in_class(Qideal I, int verbose=0);

// Return a list of lists of singular points in each ideal class.
vector<CuspList> singular_points_by_class();

// Return one list of all singular points.
CuspList singular_points();

// Return sorted list of singular points (oo, denom 2, denom 3, larger denoms in +/- pairs)
CuspList sort_singular_points(const CuspList& S, int verbose=0);

// Output sorted list of singular points (oo, denom 2, denom 3, larger denoms in +/- pairs)
void output_singular_points(const CuspList& S, int to_file=1, int to_screen=0);

// Return sorted list of (saturated covering) alphas (denom 1, denom 2, denom 3, larger denoms in pairs or fours)
// + pluspairs, minuspairs, fours
CuspList sort_alphas(const CuspList& A,
                     vector<vector<Quad>>& pluspairs, vector<vector<Quad>>& minuspairs, vector<vector<Quad>>& fours,
                     int verbose=0, int debug=0);

// Output sorted list of alphas (denom > 3 in pairs or fours)
void output_alphas(vector<vector<Quad>>& pluspairs, vector<vector<Quad>>& minuspairs, vector<vector<Quad>>& fours,
                   int to_file=1, int to_screen=0);

// direct lists of alphas of denominator 2 or 3:
CuspList denom_2_alphas();
CuspList denom_3_alphas();

// Square radius for principal cusp
RAT radius_squared(const RatQuad& a);

// For a1, a2 normalised principal cusps with circles S_ai,
// return +2, +1, 0, -1, -2:
// +2 if they do not intersect and are external to each other
// +1 if they are externally tangent
// 0  if they intersect in two distinct points
// -1 if they are internally tangent (or equal)
// -2 if they do not intersect and one is inside the other
int tau(const RatQuad& a1, const RatQuad& a2);

// return 1 iff the circles S_Ai intersect in two distinct points (not necessarily k-rational)
inline int circles_intersect(const RatQuad& a1, const RatQuad& a2) {return tau(a1,a2)==0;}

// return 1 iff the circle S_A1 is inside S_a2
int circle_inside_circle(const RatQuad& a1, const RatQuad& a2, int strict=1);

// return 1 iff the circle S_a is inside S_b for any b in blist
int circle_inside_any_circle(const RatQuad& a, const CuspList& blist, int strict=1);

// return a list of up to 2 k-rational cusps where the S_ai intersect
CuspList intersection_points_in_k(const RatQuad& a1, const RatQuad& a2);

// return 1 iff a is [strictly] inside S_b
int is_inside(const RatQuad& a, const RatQuad& b, int strict=0);

// return 1 iff a is [strictly] inside S_b for at least one b in blist
int is_inside_one(const RatQuad& a, const CuspList& blist, int strict=0);

// Return the height of S_a above z, or 0 if S_a does not cover z
RAT height_above(const RatQuad& a, const RatQuad& z);

// return -1,0,+1 according as P is over, on, under S_a (a principal)
int is_under(const H3point& P, const RatQuad& a);

// return +1 iff P is under at least one S_a for a in sliat
int is_under_any(const H3point& P, const CuspList& alist);

// Return [P] where P is the triple intersection point of the
// hemispheres S_a_i, where a0, a1, a2 are principal cusps, if there
// is one, else [].

H3pointList tri_inter_points(const RatQuad& a0, const RatQuad& a1, const RatQuad& a2);

// list of principal cusps with given denominator norm
CuspList principal_cusps_of_dnorm(const INT& n);

// list of principal cusps with denominator norm up to given bound,
// omitting any whose circles are contained in an earlier one.
CuspList principal_cusps_of_dnorm_up_to(const INT& maxn);

// list of principal cusps with given denominator
CuspList principal_cusps_with_denominator(const Quad& s);

// list of principal cusps with denominator in given list
CuspList principal_cusps_with_denominators(const vector<Quad>& slist);

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
                                    const CuspList& sigmas, int debug=0);


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

int is_alpha_surrounded(const RatQuad& a0, const CuspList& alist, const CuspList& sigmas,
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

// sigmas can be set to a list of singular points (up to translation),
// otherwise these will be computed.

// Other functions will then (1) saturate the set, (2) discard redundancies.
CuspList covering_alphas(const CuspList& sigmas, int verbose=0);

// Given a list of principal cusps alpha (all reduced mod O_k) return
// a list of "corners" P = [z,tsq] each the intersection of an S_a
// with at least two other S_{b+t} with z in the fundamental
// rectangle and tsq>0.
H3pointList triple_intersections(const CuspList& alphas, int debug=0);

// count how many P in points are on S_a
int nverts(const RatQuad& a, const H3pointList& points);

// return sublist of a in alist which have at least 3 vertices in points
CuspList remove_redundants(const CuspList& alist, const H3pointList& points);

// For P=[z,t2] in H_3, returns a list of principal cusps alpha =r/s
// such that P lies on or under S_alpha, and N(s)>=norm_s_lb.

// If option is +1 ('exact') only returns alpha for which P is on S_alpha exactly.
// If option is -1 ('strict') only returns alpha for which P is strictly under S_alpha.
// Otherwise (default), returns alpha for which P is under or on S_alpha.

CuspList covering_hemispheres(const H3point& P, int option=0, long norm_s_lb=1, int debug=0);

CuspList properly_covering_hemispheres(const H3point& P, long norm_s_lb=1, int debug=0);

// Of the properly_covering_hemispheres(P) extract the subset for which the covering height is maximal
CuspList best_covering_hemispheres(const H3point& P, long norm_s_lb=1, int debug=0);

// return max denominator norm of a list of principal cusps
INT max_dnorm(const CuspList& alphas);

// Given a list of covering alphas with max denom norm maxn, return a saturated irredundant list
CuspList saturate_covering_alphas(const CuspList& alphas, const CuspList& sigmas, INT maxn, int debug=0, int verbose=0);

// return  a saturated irredundant list of alphas in the fundamental rectangle
CuspList find_alphas(const CuspList& sigmas, int debug=0, int verbose=0);

// return  a saturated irredundant list of alphas, and list of sigmas, in the fundamental rectangle
pair<CuspList,CuspList> find_alphas_and_sigmas(int debug=0, int verbose=0);

// return list of alphas (or translates) which pass through a singular point
CuspList neighbours(const CuspList& sigmas, const CuspList& alphas);

// test if all singular points (sigmas) are surrounded by alpha circles:
int are_sigmas_surrounded(const CuspList& sigmas, const CuspList& alphas, int debug=0);
// test if one singular point (sigma) is surrounded by alpha circles:
int is_sigma_surrounded(const RatQuad& sigma, const CuspList& alphas, int debug=0);

int compare_CuspLists_as_sets(const CuspList& A, const CuspList& B);
int compare_CuspLists_as_sets_mod_translation(const CuspList& A, const CuspList& B);
#endif
