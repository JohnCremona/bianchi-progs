// FILE SWAN.H: declaration of Swan's algorithm functions

#if     !defined(_SWAN_H)
#define _SWAN_H      1       //flags that this file has been included

#include <iostream>
#include <set>

#include "ratquads.h"
#include "swan_utils.h"

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

// Return [P] where P is the triple intersection point of the
// hemispheres S_a_i, where a0, a1, a2 are principal cusps, if there
// is one, else [].

H3pointList tri_inter_points(const RatQuad& a0, const RatQuad& a1, const RatQuad& a2);

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

// test whether angle between s-->a1 and s-->a2 is <180 degrees
int angle_under_pi(const RatQuad& s, const RatQuad& a1, const RatQuad& a2);

// test whether a,b,c,d are coplanar
int coplanar(const RatQuad& a, const RatQuad& b, const RatQuad& c, const RatQuad& d);

// return list of alphas (or translates) which pass through a finite cusp
CuspList neighbours(const RatQuad& sigma, const CuspList& alphas);

// Given a base cusp s and a list of cusps alist, return a sorted
// alist with respect to the circular ordering around s
CuspList circular_sort(const RatQuad& s, const CuspList& alist);

// return sorted list of alphas (or translates) which pass through a finite cusp,
// i.e. angle_under_pi(sigma, a[i-1], a[i]) for all 0<=i<n and angle_under_pi(sigma, a[n-1], a[0]).
CuspList sorted_neighbours(const RatQuad& sigma, const CuspList& alphas);

// test if all singular points (sigmas) are surrounded by alpha circles:
int are_sigmas_surrounded(const CuspList& sigmas, const CuspList& alphas, int debug=0);
// test if one singular point (sigma) is surrounded by alpha circles:
int is_sigma_surrounded(const RatQuad& sigma, const CuspList& alphas, int debug=0);

int compare_CuspLists_as_sets(const CuspList& A, const CuspList& B);
int compare_CuspLists_as_sets_mod_translation(const CuspList& A, const CuspList& B);

// Functions which for a principal cusp alpha return a matrix M with M(alpha)=oo

// Basic version, for principal alpha: returns any matrix with M(alpha)=oo
mat22 Malpha(const RatQuad& alpha);
// Version which also ensures M(oo) is in the list alist; sets j so that M(oo)=alist[j]
mat22 Malpha(const RatQuad& alpha, const CuspList& alist, int& j);
// Version which also ensures M(s) is in the list slist; sets j so that M(s)=slist[j]
mat22 Malpha(const RatQuad& alpha, const RatQuad& s, const CuspList& slist, int& j);
// Version which also ensures M(P) is in the list Plist; sets j so that M(P)=Plist[j]
mat22 Malpha(const RatQuad& alpha, const H3point& P, const H3pointList& Plist, int& j);

inline int nverts(const POLYHEDRON& P) {return P.vertices.size();}
inline int nedges(const POLYHEDRON& P) {int n = P.edges.size(); return n/2;}
inline int nfaces(const POLYHEDRON& P) {return P.faces.size();}
// return triple (V,E,F) for a polyhedron
inline vector<int> VEF(const POLYHEDRON& P)
{
  return {nverts(P), nedges(P), nfaces(P)};
}
// return triple (V,E,F3,F4,F6) for a polyhedron
vector<int> VEFx(const POLYHEDRON& P);
// return name (e.g. "tetrahedron") of a polyhedron
string poly_name(const POLYHEDRON& P);

// given vertices and edges, fill in faces:
void fill_faces(POLYHEDRON& P, int verbose=0);

// return a tetrahedron from a list of its vertices
POLYHEDRON tetrahedron(const CuspList& V);

// return a list of tetrahedra (i.e. lists of 4 cusps (oo, sigmas[j],
// a1, a2) with a1,a2 fundamental); as a side-effect set flags[j]=1
// and flags[j']=1 where M_a_i(sigma[j])=sigma[j'] for i=1,2, for
// each.
vector<POLYHEDRON>
singular_tetrahedra(int j, const CuspList& sigmas, const CuspList& alphas, vector<int>& flags, int verbose=0);

// return a list of all singular tetrahedra
vector<POLYHEDRON>
singular_tetrahedra(const CuspList& sigmas, const CuspList& alphas, int verbose=0);

// return a polyhedron (as a list of EDGEs), from the j'th corner Plist[j]
POLYHEDRON
principal_polyhedron(int j, const CuspList& alphas, const H3pointList& Plist,
                     vector<int>& flags, int verbose=0);

// return a list of all principal polyhedra
vector<POLYHEDRON>
principal_polyhedra(const CuspList& alphas, int verbose=0);

#endif
