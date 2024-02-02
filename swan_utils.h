// FILE SWAN_UTILS.H: declaration of Swan's algorithm utility functions

#if     !defined(_SWAN_UTILS_H)
#define _SWAN_UTILS_H      1       //flags that this file has been included

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
typedef modsym EDGE;
struct POLYHEDRON {
  CuspList vertices;
  vector<EDGE> edges;
  vector<CuspList> faces;
};

ostream& operator<<(ostream& s, const H3point& P);
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

// return max denominator norm of a list of principal cusps
INT max_dnorm(const CuspList& alphas);

#endif
