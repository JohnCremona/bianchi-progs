// FILE SWAN_UTILS.CC: implementation of Swan utility functions decalred in swan_utils.h

#include <iostream>

#include "swan_utils.h"
#include "looper.h"
#include "mat22.h"

ostream& operator<<(ostream& s, const POLYHEDRON& P)
{
  s << "[V:" << P.vertices<<", E:"<<P.edges<<", F:"<<P.faces<<"]";
  return s;
}

// Square radius for principal cusp
RAT radius_squared(const RatQuad& a)
{
  return RAT(INT(1), a.den().norm());
}

// For a1, a2 normalised principal cusps with circles S_ai,
// return +2, +1, 0, -1, -2:
// +2 if they do not intersect and are external to each other
// +1 if they are externally tangent
// 0  if they intersect in two distinct points
// -1 if they are internally tangent (or equal)
// -2 if they do not intersect and one is inside the other

int tau(const RatQuad& a1, const RatQuad& a2)
{
  Quad r1=a1.num(), s1=a1.den(), r2=a2.num(), s2=a2.den();
  INT n1 = s1.norm(), n2 = s2.norm(), n3 = mms(r1,s2,r2,s1).norm();
  INT x = n3-n1-n2;
  int sd1 = sign(x);
  int sd2 = sign(x*x-4*n1*n2);
  int ans = (sd2 < 0? 0: sd1 * (sd2==0? 1: 2));
  return ans;
}

// return 1 iff the circle S_A1 is inside S_a2
int circle_inside_circle(const RatQuad& a1, const RatQuad& a2, int strict)
{
  int s = sign(a1.den().norm() - a2.den().norm());
  if ((strict? s<0: s<=0))
    return 0;
  int t  = tau(a1,a2);
  return (strict? t==-2: t<0);
}

// return 1 iff the circle S_a is inside S_b for any b in blist
int circle_inside_any_circle(const RatQuad& a, const CuspList& blist, int strict)
{
  return std::any_of(blist.begin(), blist.end(),
                     [a, strict](const RatQuad& b) {return circle_inside_circle(a,b,strict);});
}

// return a list of up to 2 k-rational cusps where the S_ai intersect
CuspList intersection_points_in_k(const RatQuad& a1, const RatQuad& a2)
{
  RAT r1sq = radius_squared(a1), r2sq = radius_squared(a2);
  RatQuad delta = a2-a1;
  RAT n = delta.norm();
  RAT d1 = n - (r1sq + r2sq);
  RAT d2 = d1*d1 - 4*r1sq*r2sq;

  CuspList ans;
  if (d2 > 0)
    return ans;
  delta = delta.conj();
  RatQuad z = a1 + a2 + RatQuad(r1sq-r2sq)/delta;
  // find sqrts of d2 in k, if any
  for ( const auto& sqrtd2 : sqrts_in_k(d2))
    {
      RatQuad pt = (z + sqrtd2/delta)/2;
      assert ((pt-a1).norm() == r1sq);
      assert ((pt-a2).norm() == r2sq);
      ans.push_back(pt);
    }
  return ans;
}

// return the point on the intersection of the hemispheres S_a_i
// (where a1, a2 are principal cusps) which is on the line from a1 to
// a2.  Used when both S_a_i pass through a singular point.
H3point bi_inter(const RatQuad& a1, const RatQuad& a2)
{
  RAT r1sq = radius_squared(a1), r2sq = radius_squared(a2);
  RatQuad delta = a2-a1;
  RAT n = delta.norm();
  RAT la = (1+(r2sq-r1sq)/n)/2;
  RAT mu = 1 - la;
  RatQuad z = la*a1 + mu*a2;
  RAT t1 = r1sq-mu*mu*n;
  RAT t2 = r2sq-la*la*n;
  assert (t1==t2);
  assert (t2>0);
  H3point P = {z,t2};
  return P;
}

// return 1 iff a is on (the boundary of) S_b  (b principal, a arbitrary)
int is_on(const RatQuad& a, const RatQuad& b)
{
  int t = sign(a.den().norm() - (a.num()*b.den()-a.den()*b.num()).norm());
  return t==0;
}

// return 1 iff a is [strictly] inside S_b (b principal, a arbitrary)
int is_inside(const RatQuad& a, const RatQuad& b, int strict)
{
  int t = sign(a.den().norm() - (a.num()*b.den()-a.den()*b.num()).norm());
  return (strict ? 0<t : 0<=t);
}

// return 1 iff a is [strictly] inside S_b for at least one b in blist
int is_inside_one(const RatQuad& a, const CuspList& blist, int strict)
{
  return std::any_of(blist.begin(), blist.end(),
                     [a, strict](const RatQuad& b)
                     {return is_inside(a,b,strict);});
}

// return 1 iff a1 and a2 are both [strictly] inside the same S_b for at least one b in blist
int are_both_inside_one(const RatQuad& a1, const RatQuad& a2, const CuspList& blist, int strict)
{
  return std::any_of(blist.begin(), blist.end(),
                     [a1, a2, strict](const RatQuad& b)
                     {return is_inside(a1,b,strict) && is_inside(a2,b,strict);});
}

// list of principal cusps with given denominator norm
CuspList principal_cusps_of_dnorm(const INT& n)
{
  return principal_cusps_with_denominators(quads_of_norm(n));
}

// list of principal cusps with denominator norm up to given bound,
// omitting any whose circles are contained in an earlier one.

CuspList principal_cusps_of_dnorm_up_to(const INT& maxn)
{
  return principal_cusps_with_denominators(quads_of_norm_up_to(maxn, 1, 1));
}

// list of principal cusps with given denominator
CuspList principal_cusps_with_denominator(const Quad& s)
{
  auto invres = invertible_residues(s);
  CuspList alist(invres.size());
  std::transform(invres.begin(), invres.end(), alist.begin(),
                 [s](const Quad& r) {return reduce_to_rectangle(RatQuad(r, s));});
  return alist;
}

// list of principal cusps with denominator in given list
CuspList principal_cusps_with_denominators(const vector<Quad>& slist)
{
  CuspList alist;
  for ( const auto& s : slist)
    {
      // cout << "denominator "<<s<<" (norm "<<s.norm()<<"): "<<flush;
      auto alist1 = principal_cusps_with_denominator(s);
      // cout << alist1 << endl;
      alist.insert(alist.end(), alist1.begin(), alist1.end());
    }
  return alist;
}

void test_principal_cusps(int n1, int n2)
{
  for (int n=1; n<=n1; n++)
    {
      auto alist = principal_cusps_of_dnorm(INT(n));
      cout << "Principal cusps of denominator norm "<<n<<": "<<alist<<endl;
    }
  auto alist = principal_cusps_of_dnorm_up_to(INT(n2));
  cout << "Principal cusps of denominator norm up to "<<n2<<": "<<alist<<endl;
}

// det([[a1,a2,a3],[a1bar,a2bar,a3bar],[1,1,1]])
RatQuad tri_det(const RatQuad& a1, const RatQuad& a2, const RatQuad& a3)
{
  RatQuad a1b=a1.conj(), a2b=a2.conj(), a3b=a3.conj();
  return (a2b - a3b)*a1 + (a3b - a1b)*a2 + (a1b - a2b)*a3;
}

// Return [P] where P is the triple intersection point of the
// hemispheres S_a_i, where a0, a1, a2 are principal cusps, if there
// is one, else [].

H3pointList tri_inter_points(const RatQuad& a0, const RatQuad& a1, const RatQuad& a2)
{
  H3pointList points;
  RatQuad delta = tri_det(a0,a1,a2);
  if (delta.is_zero())
    return points; // empty
  delta = delta.conj(); // to match Sage code
  RAT
    rho0(radius_squared(a0)),
    rho1(radius_squared(a1)),
    rho2(radius_squared(a2)),
    n0(a0.norm()),
    n1(a1.norm()),
    n2(a2.norm());
  RatQuad z = (a1*(n0-n2+rho2-rho0) + a2*(n1-n0+rho0-rho1) + a0*(n2-n1+rho1-rho2)) / delta;
  RatQuad zbar = z.conj();
  RAT znorm = z.norm();
  RAT t2 =  (rho0 - n0 - znorm) + 2*(a0*zbar).x_coord(1);
  RAT t2a = (rho1 - n1 - znorm) + 2*(a1*zbar).x_coord(1);
  RAT t2b = (rho2 - n2 - znorm) + 2*(a2*zbar).x_coord(1);
  assert (t2==t2a);
  assert (t2==t2b);
  H3point P = {z,t2};
  assert (is_on(P,a0));
  assert (is_on(P,a1));
  assert (is_on(P,a2));
  if (t2>=0)
    points.push_back(P);
  return points;
}

// return max denominator norm of a list of principal cusps
INT max_dnorm(const CuspList& alphas)
{
  INT m;
  std::for_each(alphas.begin(), alphas.end(),
                [&m](const RatQuad& a) {INT n = a.den().norm(); if (n>m) m=n;});
  return m;
}

// test whether angle between s-->a1 and s-->a2 is <180 degrees
int angle_under_pi(const RatQuad& s, const RatQuad& a1, const RatQuad& a2)
{
  return sign_im_cr(a2,s,a1)>0;
}

// test whether a,b,c,d are coplanar
int coplanar(const RatQuad& a, const RatQuad& b, const RatQuad& c, const RatQuad& d)
{
  return sign_im_cr(a,b,c,d)==0;
}

int compare_CuspLists_as_sets(const CuspList& A, const CuspList& B)
{
  return (A.size()==B.size()) &&
    std::all_of(A.begin(), A.end(), [B](const RatQuad& a) {return std::count(B.begin(), B.end(), a)==1;});
}

int compare_CuspLists_as_sets_mod_translation(const CuspList& A, const CuspList& B)
{
  Quad t;
  return (A.size()==B.size()) &&
    std::all_of(A.begin(), A.end(),
                [B,&t](const RatQuad& a) {return cusp_index_with_translation(a, B, t) !=-1;});
}

// Functions which for a principal cusp alpha return a matrix M with M(alpha)=oo

// Basic version, for principal alpha: returns any matrix with M(alpha)=oo
mat22 Malpha(const RatQuad& alpha)
{
  Quad a, b, c = alpha.den(), d = alpha.num();
  if (c==Quad::zero)
    return mat22::identity;
  if (d==Quad::zero)
    return mat22::S;
  Quad g = quadbezout(c, d, b, a); // a*d+b*c=g=1
  if (g==Quad::zero)
    cerr<<"Error in Malpha(): cusp "<<alpha<<" is not principal"<<endl;
  assert(a*d+b*c==g);
  if (!g.is_one()) cout<<"Malpha("<<alpha<<"): a="<<a<<", b="<<b<<", g="<<g<<endl;
  assert(g.is_one());
  mat22 M(-a,-b,c,-d);
  assert (M.is_unimodular()); // strict by default, i.e. det=1 not just a unit
  assert (M(alpha)==RatQuad::infinity());
  return M;
}

// Version which also ensures M(oo) is in the list alist; sets j so that M(oo)=alist[j]
mat22 Malpha(const RatQuad& alpha, const CuspList& alist, int& j)
{
  return Malpha(alpha, RatQuad::infinity(), alist, j);
}

//#define DEBUG_M_ALPHA

// Version which also ensures M(s) is in the list slist; sets j so that M(s)=slist[j]
mat22 Malpha(const RatQuad& alpha, const RatQuad& s, const CuspList& slist, int& j)
{
#ifdef DEBUG_M_ALPHA
  cout << "In Malpha("<<alpha<<") with s="<<s<<" and slist="<<slist<<endl;
#endif
  mat22 M = Malpha(alpha);
  Quad x;
  RatQuad sd = M(s);
#ifdef DEBUG_M_ALPHA
  cout << " - first M = "<<M<<" with M(s)="<<sd<<endl;
  cout<<"Differences of this and elements of slist:\n";
  for (const auto& a : slist) cout<<(sd-a)<<endl;
#endif
  j = cusp_index_with_translation(sd, slist, x);
  assert (j>=0);
  M = mat22::Tmat(-x)*M;
  assert (M(alpha).is_infinity());
  assert (M(s) == slist[j]);
  return M;
}

// Version which also ensures M(P) is in the list Plist; sets j so that M(P)=Plist[j]
mat22 Malpha(const RatQuad& alpha, const H3point& P, const H3pointList& Plist, int& j)
{
  mat22 M = Malpha(alpha);
  Quad x(0);
  H3point Q = M(P);
  j = point_index_with_translation(Q, Plist, x);
  assert (j>=0);
  M = mat22::Tmat(-x)*M;
  if (!(M(alpha).is_infinity()))
    {
      cout<<"alpha = "<<alpha<<", M="<<M<<": M(alpha)="<<M(alpha)<<endl;
    }
  assert (M(alpha).is_infinity());
  Q = M(P);
  assert (Q == Plist[j]);
  return M;
}

// Given a full list of alphas, return the list of all M_alphas (such
// that M_alpha(alpha)=oo and M_alpha(oo)=alpha' is in the list).
// Also sets inv s.t. M_alpha[i](oo) = alpha[inv[i]].
vector<mat22> all_M_alphas(const CuspList& alist, vector<int>& inv)
{
  vector<mat22> Mlist; Mlist.reserve(alist.size());
  inv.clear();  inv.reserve(alist.size());
  for (const auto& a : alist)
    {
      int j;
      mat22 M = Malpha(a, alist, j); // i.e. M(a)=oo and M(oo)=alist[j]
      assert (j>=0);
      Mlist.push_back(M);
      inv.push_back(j);
    }
  return Mlist;
}

// test for being a valid edge M{alpha,oo} with alpha=alphas[j]
int valid_edge(const RatQuad& a, const RatQuad& b, const CuspList& alphas,
               mat22& M, int& j)
{
  if (!b.is_principal()) return 0;
  M = Malpha(b);
  Quad u;
  j = cusp_index_with_translation(M(a), alphas, u);
  if (j<0)
    return 0;
  // Now M(a) = u+alphas[j],
  // so {a,b} = M^(-1){u+alphas[j],oo} = M^(-1)*T^u{alphas[j],oo}
  M = M.inverse() * mat22::Tmat(u);
  return 1;
}

// count how many P in points are on S_a
int nverts(const RatQuad& a, const H3pointList& points)
{
  return std::count_if(points.begin(), points.end(),
                       [a](const H3point& P) { return is_on(P,a);} );
}

// POLYHEDRON utilities

// return triple (V,E,F3,F4,F6) for a polyhedron
vector<int> VEFx(const POLYHEDRON& P)
{
  vector<int> ans = {nverts(P), nedges(P), 0, 0, 0};
  for (const auto& f : P.faces)
    {
      int n = f.size();
      int i = (n==3? 2 : (n==4? 3 : (n==6? 4 : -1)));
      if (i>0)
        ans[i]++;
      else
        cout << "Polyhedron has a face "<<f<<" with "<<n<<" faces!"<<endl;
    }
  return ans;
}

static const map<vector<int>, string> poly_names =
  {
    {{4,6,4,0,0}, "tetrahedron"},
    {{8,12,0,6,0}, "cube"},
    {{6,12,8,0,0}, "octahedron"},
    {{6,9,2,3,0}, "triangular prism"},
    {{5,8,4,1,0}, "square pyramid"},
    {{12,18,0,6,2}, "hexagonal prism"},
    {{9,15,4,3,1}, "hexagonal cap"},
    {{12,18,4,0,4}, "truncated tetrahedron"},
    {{12,24,8,6,0}, "cuboctahedron"},
    {{8,14,4,4,0}, "sliced cube"},
    {{18,30,0,12,2}, "double hexagonal prism"},
    {{7,14,8,1,0}, "half star"},
    {{18,33,6,9,2}, "tented hexagonal prism"},
    {{5,9,6,0,0}, "3-dipyramid"},
    {{7,15,10,0,0}, "5-dipyramid"},
    {{8,18,12,0,0}, "6-dipyramid"},
    {{7,13,6,2,0}, "triangular prism plus square pyramid"},
    {{6,11,6,1,0}, "tetrahedron plus square pyramid"},
    {{12,20,0,10,0}, "DoubleCube"},
    {{12,22,4,8,0}, "DoubleTentedCube"},
    {{9,19,10,2,0}, "unnamed #1"},
    {{16,32,8,10,0}, "unnamed #2"},
    {{7,11,2,4,0}, "unnamed #3"}
  };

string poly_name(const POLYHEDRON& P)
{
  auto vef = VEFx(P);
  // return poly_names.at(vef); // will crash if P is not in the above list
  auto search = poly_names.find(vef);
  if ( search!=poly_names.end() )
    return search->second;
  cout << "Polyhedron with (V,E,F3,F4,F6) = " << vef << " not recognised" << endl;
  // int nv = vef[0];
  // int n = nv-2;
  // if (vef[1]==3*n && vef[2]==2*n && vef[3]==0 && vef[4]==0)
  //   cout << " -- looks like a "<<n<<"-dipyramid";
  // cout <<endl;
  return "new type of polyhedron";
}

// Given a base cusp s and a list of cusps alist, return a sorted
// alist with respect to the circular ordering around b
CuspList circular_sort(const RatQuad& s, const CuspList& alist)
{
  int n = alist.size();
  assert (n>=3);
  if (n==3) return alist;
  int debug=0;
  if (debug)
    cout<<"circular sort of "<<alist<<" around "<<s<<endl;
  CuspList sorted_alist;
  RatQuad a = alist[0];
  // can start with any a for circular sort:
  sorted_alist.push_back(a);
  if (debug)
    cout<<" - initial list: "<<sorted_alist<<endl;

  auto is_after = [s](const RatQuad& a, const RatQuad& b) {return angle_under_pi(s,a,b);};
  auto half_afters = [alist, is_after](const RatQuad& a) {
		       CuspList afters;
		       std::copy_if(alist.begin(), alist.end(), std::back_inserter(afters),
				    [a, is_after](const RatQuad& b){return is_after(a,b);});
		       return afters;
		     };
  auto is_first = [is_after](const CuspList& as, const RatQuad& a) {
		    return std::all_of(as.begin(), as.end(), [is_after, a](const RatQuad& b) {return b==a || is_after(a,b);});
		      };
  auto find_first = [is_first](const CuspList& as) {
		      auto it = std::find_if(as.begin(), as.end(), [is_first, as](const RatQuad& a) {return is_first(as,a);});
		      return *it;
		    };
  auto next = [find_first, half_afters](const RatQuad& a) { return find_first(half_afters(a)); };

  // Now repeatedly append the next one:
  while (sorted_alist.size() < alist.size())
    {
      sorted_alist.push_back(next(sorted_alist.back()));
      if (debug)
	cout<<" - next list: "<<sorted_alist<<endl;
    }

  if (debug)
    cout<<" - final list is "<<sorted_alist<<"\n ... now checking..."<<endl;
  // Check:
  for (int i=0; i<n; i++)
    {
      int i1 = (i+1)%n;
      int t = angle_under_pi(s, sorted_alist[i], sorted_alist[i1]);
      if (debug || !t)
	cout<<"angle_under_pi("<<s<<", "<< sorted_alist[i] << ", " << sorted_alist[i1] << ") = " << t << endl;
      assert(t);
    }

  return sorted_alist;
}

// Return the height of S_a above P, or 0 if S_a does not cover P
RAT height_above(const RatQuad& a, const RatQuad& z)
{
  return radius_squared(a)-(z-a).norm(); // positive iff S_a covers z
}

// return -1,0,+1 according as P is over, on, under S_a (a principal)
int is_under(const H3point& P, const RatQuad& a)
{
  return sign(height_above(a, P.z) - P.t2);
}

// return 1 iff P is strictly under S_a
// (same as is_under(P,a)==1 but faster)
int is_strictly_under(const H3point& P, const RatQuad& a)
{
  const Quad& s = a.den(), d = P.z.den();
  INT x = d.norm() - mms(s, P.z.num(), a.num(), d).norm();
  if (x<0) return 0;
  return (s * d).norm() * P.t2.num() < x * P.t2.den();
}

// return 1 iff P is exactly on S_a
// (same as is_under(P,a)==0 but faster)
int is_on(const H3point& P, const RatQuad& a)
{
  const Quad& s = a.den(), d = P.z.den();
  INT quo, rem;
  if (!divrem((s * d).norm() * P.t2.num(), P.t2.den(), quo, rem)) return 0;
  return d.norm() == mms(s, P.z.num(), a.num(), d).norm() + quo;
}


// return +1 iff P is under at least one S_a for a in alist
int is_under_any(const H3point& P, const CuspList& alist)
{
  return std::any_of(alist.begin(), alist.end(),
                     [P](const RatQuad& a) {return is_strictly_under(P,a);});
}

// For z in F4 (quarter rectangle) return list of z and the 8
// neighbours of z, -z, zbar, -zbar in the 8 surrounding
// quarter-rectangles.  These are (with cz=z.conj()): {z, -z, cz, 1-z,
// -cz, 1-cz, w-z, w+cz} and either w+cz-1 or w+1-cz depending on
// Trace(w).

CuspList oldF4nbrs(const RatQuad& z)
{
  RatQuad cz = z.conj();
  int t = Quad::t;
  Quad w = Quad::w, one(1);
  CuspList ans = {z};
  for ( RatQuad a : {-z, cz, one-z, -cz, one-cz, w-z, w+cz, (t? w+cz-one : w+one-z)})
    if (std::find(ans.begin(), ans.end(), a) == ans.end())
      ans.push_back(a);
  return ans;
}

// For principal a in F4 (quarter rectangle) return list of a and up
// to 3 neighbours of a, -a, abar, -abar in the 8 surrounding
// quarter-rectangles, omitting those whose circles do not intersect
// F4 so cannot contribut to a triple intersection.

// Write a=r/s= x+y*sqrt(-d), n=N(s), r*conj(s) = n*a = u+v*sqrt(-d).
// Then
// (L) a is near left   edge iff u^2 <= n;
// (R) a is near right  edge iff (n-2*u)^2 <= 4*n;
// (B) a is near bottom edge iff d*v^2 <= n;
// (T) a is near top    edge iff d*(n-2*v)^2 <=  4*n (t==0)
//                            or d*(n-4*v)^2 <= 16*n (t==1)
//
// if (L) we add -abar;
// if (R) we add 1-abar;
// if (B) we add abar;
// if (T) we add w+abar (t==0) or w-a (t==1);
// if (L) & (B) we also add -a;
// if (B) & (R) we also add 1-a;
// if (R) & (T) we also add 1+w-a (t==0) or w+abar   (t==1);
// if (T) & (L) we also add   w-a (t==0) or w-1+abar (t==1).

// These are (with
// cz=z.conj()): {z, -z, cz, 1-z, -cz, 1-cz, w-z, w+cz} and either
// w+cz-1 or w+1-cz depending on Trace(w).
CuspList F4nbrs(const RatQuad& a)
{
  RatQuad abar = a.conj();
  int t = Quad::t, d = Quad::d;
  Quad w = Quad::w, r = a.num(), s = a.den(), one(1);
  INT i = racb(r, s), j = iacb(r, s), n = s.norm(); // so r*conj(s) = i+j*w
  INT p = (t? 2*i+j: 2*i), q = (t? j: 2*j);          //   = (p*q*sqrt(-d))/2
  INT p2 = n-p, q2 = n-2*q, n4 = 4*n;
  int L = (p*p<=n4), R = p2*p2<=n4, B = d*q*q <= n4, T = d*q2*q2 <= 4*n4;

  std::set<RatQuad> nbrs = {a};
  if (L) nbrs.insert(-abar);
  if (R) nbrs.insert(one-abar);
  if (B) nbrs.insert(abar);
  if (T) nbrs.insert((t? w-a : w+abar));
  if (L&&B) nbrs.insert(-a);
  if (B&&R) nbrs.insert(one-a);
  if (R&&T) nbrs.insert((t? w+abar : w+1-a));
  if (T&&L) nbrs.insert((t? w-1+abar : w-a));

  CuspList ans;
  std::copy(nbrs.begin(), nbrs.end(), std::back_inserter(ans));
  return ans;
}

// For P=[z,t2] in H_3, returns a list of principal cusps alpha =r/s
// such that P lies on or under S_alpha, and N(s)>=norm_s_lb.  For
// each s with norm_s_lb <= N(s) <= 1/t2 the only candidate(s) is r
// such that N(s*z-r)<1, since the inequality to be satisfied is
// N(s*z-r)+t2*N(s)<=1.

// If option is +1 ('exact') only returns alpha for which P is on S_alpha exactly.
// If option is -1 ('strict') only returns alpha for which P is strictly under S_alpha.
// Otherwise (default), returns alpha for which P is under or on S_alpha.

CuspList covering_hemispheres(const H3point& P, int option, long norm_s_lb, int debug)
{
  if (debug)
    {
      cout << "Finding a for which S_a covers P = "<<P<<endl;
      cout << " option = "<<option;
      if (option==1)
        cout<<" (exact, i.e. P is on S_a)";
      else
        if (option==-1)
          cout<<" (strict, i.e. P is strictly under S_a)";
        else cout<<" (default, i.e. P is on or under S_a)";
      cout <<endl;
    }
  CuspList ans;
  RatQuad z = P.z, sz;
  RAT t2 = P.t2;
  Quad temp;
  int ok, test;
  INT norm_s_ub = t2.recip().floor();
  if (debug)
    cout << "t2 = "<<t2<<" so bounds on N(s) are ["<<norm_s_lb<<","<<norm_s_ub<<"]\n";

  // We could do Quadlooper sloop(norm_s_lb, norm_s_ub, 1); but it's
  // more efficient to construct all possible s at once, though that
  // takes more memory
  auto slist = quads_of_norm_between(INT(norm_s_lb), norm_s_ub, 1, 0); // including conjugates, not sorted
  for (const auto& s : slist)
    {
      if (debug)
        cout << " s = "<<s<<" ";
      sz = s*z;
      auto rlist = nearest_quads(sz, 0); // 0 means all (there can be 2), 1 means at most one
      if (debug)
        cout << " possible r: "<<rlist<<endl;
      for (const auto& r : rlist)
        {
          test = sign((sz-r).norm() + t2*s.norm() - ONE);
          //  0 for on, -1 for strictly under, +1 for over
          if (debug)
            cout << " r = "<<r<<"\n test = "<<test<< endl;
          switch (option) {
          case 1:
            ok = (test==0); break;
          case -1:
            ok = (test==-1); break;
          default:
            ok = (test<1);
          }
          if (debug)
            cout << " ok = "<<ok<< endl;
          if (!ok)
            continue;
          ok = coprime(r,s);
          if (debug)
            cout << " coprime(r,s)? "<< ok << endl;
          if (!ok)
            continue;
          if (debug)
            cout << " +++ Success with a = "<<RatQuad(r,s)<<endl;
          ans.push_back(RatQuad(r,s));
        }
    }
  if (debug)
    cout << "Covering hemispheres are S_a for a in "<<ans<<endl;
  return ans;
}

CuspList properly_covering_hemispheres(const H3point& P, long norm_s_lb, int debug)
{
  return covering_hemispheres(P, -1, norm_s_lb, debug);
}

// multiply a point by fundamental unit (usually -1, hence the name here)
H3point negate(const H3point& P)
{
  return {fundunit * P.z, P.t2};
}

H3point translate(const H3point& P, const Quad& t)
{
  return {P.z + t, P.t2};
}

// Of the properly_covering_hemispheres(P) extract the subset for which the covering height is maximal
CuspList best_covering_hemispheres(const H3point& P, long norm_s_lb, int debug)
{
  // Find all covering hemispheres
  CuspList alist = properly_covering_hemispheres(P, norm_s_lb, debug);
  if (alist.empty())
    return alist;

  if(debug)
    cout<<"S_a covering P="<<P<<":"<<endl;

  // Find the max height of all S_a above P
  RAT m(0);
  if(debug)
    std::for_each(alist.begin(), alist.end(),
                  [P](const RatQuad& a) {cout<<"a="<<a<<": height above P is "<<height_above(a,P.z)<<endl;});
  std::for_each(alist.begin(), alist.end(),
                [P,&m](const RatQuad& a) {m = max(m, height_above(a,P.z));});
  if(debug)
    cout<<"max height is "<<m<<endl;
  // Discard those a whose height above P is not maximal
  alist.erase(std::remove_if(alist.begin(), alist.end(),
                              [P,m](const RatQuad& a) { return height_above(a,P.z) < m;}),
               alist.end());
  if(debug)
    {
      cout<<"S_a maximally covering P="<<P<<": "<<alist<<" (ht "<<m<<" above P)"<<endl;
      // for (const auto& a: alist)
      //   cout<<a<<" with height "<<height_above(a,P.z)<<endl;
    }
  return alist;
}

H3point mat22::operator()(const H3point& P) const
{
  RatQuad z = P.z;
  RAT t2 = P.t2;
  RAT  n = (c*z+d).norm() + c.norm()*t2;
  RatQuad new_z = ((a*z+b)*(c*z+d).conj() + a*c.conj()*t2) / n;
  RAT new_t2 = t2 / (n*n);
  return {new_z, new_t2};
}

// return 1 iff P is an integer translate of Q, with t=P-Q
int is_translate(const H3point& P, const H3point& Q, Quad& t)
{
  return ( (P.t2 == Q.t2) &&
           integral_difference(P.z, Q.z, t) );
}

// Return index i of P mod O_K in Plist, with t=P-Plist[i], or -1 if not in list
int point_index_with_translation(const H3point& P, const vector<H3point>& Plist, Quad& t)
{
  int j = 0;
  for ( const auto& Q : Plist)
    {
      if (is_translate(P, Q, t))
        return j;
      j++;
    }
  return -1;
}

// return whether the cusp is finite singular
int is_cusp_singular(const RatQuad& a, const CuspList& sigmas)
{
  Quad x;
  return (a.is_finite()) && cusp_index_with_translation(a, sigmas, x)>=0;
}

// return number of vertices which are finite singular
int is_face_singular(const CuspList& face, const CuspList& sigmas)
{
  return std::count_if(face.begin(), face.end(),
                       [sigmas](const RatQuad& a) {
                         return is_cusp_singular(a, sigmas);
                       });
}

CuspList cusp_shifts(const RatQuad& a, const vector<Quad>& sh)
{
  CuspList blist(sh.size());
  std::transform(sh.begin(), sh.end(), blist.begin(),
                 [a] (const Quad& s) { return a+s; });
  return blist;
}

CuspList cusp_shifts(const CuspList& alist, const vector<Quad>& sh)
{
  CuspList blist;
  for ( const auto& a : alist)
    {
      auto a_sh = cusp_shifts(a, sh);
      blist.insert(blist.end(), a_sh.begin(), a_sh.end());
    }
  return blist;
}

H3pointList H3point_shifts(const H3point& P, const vector<Quad>& sh)
{
  H3pointList Qlist(sh.size());
  std::transform(sh.begin(), sh.end(), Qlist.begin(),
                 [P] (const Quad& s) {return translate(P,s);});
  return Qlist;
}

// Given a principal cusp a0, and a list of others, alist, return the
// list of b in alist which intersect a0 in 2 distinct points.
CuspList intersecting_alphas(const RatQuad& a0, const CuspList& alist)
{
  CuspList blist;
  std::copy_if(alist.begin(), alist.end(), std::back_inserter(blist),
               [a0](const RatQuad& b){return circles_intersect(a0, b);});
  return blist;
}
