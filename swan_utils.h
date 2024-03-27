// FILE SWAN_UTILS.H: declaration of Swan's algorithm utility functions

#if     !defined(_SWAN_UTILS_H)
#define _SWAN_UTILS_H      1       //flags that this file has been included

#include <iostream>
#include <set>

#include "mat22.h"
#include "geometry.h"
#include "ratquads.h"
#include "looper.h"

typedef modsym EDGE;
struct POLYHEDRON {
  CuspList vertices;
  vector<EDGE> edges;
  vector<CuspList> faces;
};

ostream& operator<<(ostream& s, const POLYHEDRON& P);

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

// return the point on the intersection of the hemispheres S_a_i
// (where a1, a2 are principal cusps) which is on the line from a1 to
// a2.  Used when both S_a_i pass through a singular point.
H3point bi_inter(const RatQuad& a1, const RatQuad& a2);

// return 1 iff a is on (the boundary of) S_b  (b principal, a arbitrary)
int is_on(const RatQuad& a, const RatQuad& b);

// return 1 iff a is [strictly] inside S_b  (b principal, a arbitrary)
int is_inside(const RatQuad& a, const RatQuad& b, int strict=0);

// return 1 iff a is [strictly] inside S_b for at least one b in blist
int is_inside_one(const RatQuad& a, const CuspList& blist, int strict=0);

// Return the height of S_a above z, or 0 if S_a does not cover z
RAT height_above(const RatQuad& a, const RatQuad& z);

// return -1,0,+1 according as P is over, on, under S_a (a principal)
int is_under(const H3point& P, const RatQuad& a);

// return +1 iff P is under at least one S_a for a in slist
int is_under_any(const H3point& P, const CuspList& alist);

// multiply a point by fundamental unit (usually -1, hence the name here)
H3point negate(const H3point& P);

// translate a point
H3point translate(const H3point& P, const Quad& t);

// return 1 iff P is an integer translate of Q, with t=P-Q
int is_translate(const H3point& P, const H3point& Q, Quad& t);

// Return index i of P mod O_K in Plist, with t=P-Plist[i], or -1 if not in list
int point_index_with_translation(const H3point& P, const vector<H3point>& Plist, Quad& t);

// list of principal cusps with given denominator norm
CuspList principal_cusps_of_dnorm(const INT& n);

// list of principal cusps with denominator norm up to given bound,
// omitting any whose circles are contained in an earlier one.
CuspList principal_cusps_of_dnorm_up_to(const INT& maxn);

// list of principal cusps with given denominator
CuspList principal_cusps_with_denominator(const Quad& s);

// list of principal cusps with denominator in given list
CuspList principal_cusps_with_denominators(const vector<Quad>& slist);

// det([[a1,a2,a3],[a1bar,a2bar,a3bar],[1,1,1]])
RatQuad tri_det(const RatQuad& a1, const RatQuad& a2, const RatQuad& a3);

// return max denominator norm of a list of principal cusps
INT max_dnorm(const CuspList& alphas);

// test whether angle between s-->a1 and s-->a2 is <180 degrees
int angle_under_pi(const RatQuad& s, const RatQuad& a1, const RatQuad& a2);

// test whether a,b,c,d are coplanar
int coplanar(const RatQuad& a, const RatQuad& b, const RatQuad& c, const RatQuad& d);

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

// Given a full list of alphas, return the list of all M_alphas (such
// that M_alpha(alpha)=oo and M_alpha(oo)=alpha' is in the list).
// Also sets inv s.t. M_alpha[i](oo) = alpha[inv[i]].
vector<mat22> all_M_alphas(const CuspList& alist, vector<int>& inv);

// test for being a valid edge M{alphas[j],oo}
int valid_edge(const RatQuad& a, const RatQuad& b, const CuspList& alphas,
               mat22& M, int& j);

// polyhedron utilities
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

// Return [P] where P is the triple intersection point of the
// hemispheres S_a_i, where a0, a1, a2 are principal cusps, if there
// is one, else [].

H3pointList tri_inter_points(const RatQuad& a0, const RatQuad& a1, const RatQuad& a2);

// Given a base cusp s and a list of cusps alist, return a sorted
// alist with respect to the circular ordering around s
CuspList circular_sort(const RatQuad& s, const CuspList& alist);

// Return the height of S_a above P, or 0 if S_a does not cover P
RAT height_above(const RatQuad& a, const RatQuad& z);
// return -1,0,+1 according as P is over, on, under S_a (a principal)
int is_under(const H3point& P, const RatQuad& a);
// return +1 iff P is under at least one S_a for a in alist
int is_under_any(const H3point& P, const CuspList& alist);

// For z in F4 (quarter rectangle) return list of z and the 8
// neighbours of z, -z, zbar, -zbar in the 8 surrounding
// quarter-rectangles
CuspList F4nbrs(const RatQuad& z);

// For P=[z,t2] in H_3, returns a list of principal cusps alpha =r/s
// such that P lies on or under S_alpha, and N(s)>=norm_s_lb.

// If option is +1 ('exact') only returns alpha for which P is on S_alpha exactly.
// If option is -1 ('strict') only returns alpha for which P is strictly under S_alpha.
// Otherwise (default), returns alpha for which P is under or on S_alpha.

CuspList covering_hemispheres(const H3point& P, int option=0, long norm_s_lb=1, int debug=0);

CuspList properly_covering_hemispheres(const H3point& P, long norm_s_lb=1, int debug=0);

// Of the properly_covering_hemispheres(P) extract the subset for which the covering height is maximal
CuspList best_covering_hemispheres(const H3point& P, long norm_s_lb=1, int debug=0);

// return whether the cusp is finite singular
int is_cusp_singular(const RatQuad& a, const CuspList& sigmas);

// return number of vertices which are finite singular
int is_face_singular(const CuspList& face, const CuspList& sigmas);

#endif
