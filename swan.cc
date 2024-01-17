// FILE SWAN.CC: implementation of Swan's algorithm

#include <iostream>

#include "swan.h"
#include "geometry.h"
#include "looper.h"

H3_comparison H3_cmp;

// Given an ideal I, return a list of singular points of class [I]
// (one representative for each orbit under integral translations).

CuspList singular_points_in_class(Qideal I, int verbose)
{
  if (I.is_principal())
    return {RatQuad(Quad::one, Quad::zero)};
  INT n = I.norm();
  Quad temp, r, s(n);
  vector<Quad> slist = {s};
  Quad s2 = I.zgen(1);
  if (s2.norm()==n*n)
    slist.push_back(s2);
  else
    {
      s2 = n-s2;
      if (s2.norm()==n*n)
        slist.push_back(s2);
    }
  if (verbose)
    cout<<" I = "<<I<<" has minimal elememts "<<slist<<endl;
  CuspList S;
  for ( const auto& s : slist)
    {
      if (verbose)
        cout<<" - using s = "<<s<<endl;
      auto rlist = Qideal(s).residues();
      if (verbose)
        cout<<"Residues modulo "<<s<<" are "<<rlist<<endl;
      for ( auto r : rlist)
        {
          if (I==Qideal({r,s}))
            {
              RatQuad sig = reduce_to_rectangle(RatQuad(r,s), temp);
              assert (sig.in_rectangle());
              if (verbose)
                cout<<" - using r = "<<r<<" to give "<<sig<<endl;
              S.push_back(sig);
            }
        }
    }
  return S;
}

// Return a list of lists of singular points in each ideal class.

vector<CuspList> singular_points_by_class()
{
  vector<CuspList> sigma_lists;
  for ( const auto& I : Quad::class_group)
    sigma_lists.push_back(singular_points_in_class(I));
  return sigma_lists;
}

// Return one list of all singular points.

CuspList singular_points()
{
  CuspList sigma_list;
  for ( const auto& I : Quad::class_group)
    {
      CuspList S = singular_points_in_class(I);
      sigma_list.insert(sigma_list.end(), S.begin(), S.end());
    }
  return sigma_list;
}

// Return sorted list of singular points (oo, denom 2, denom 3, larger denoms in +/- pairs)

CuspList sort_singular_points(const CuspList& slist, int verbose)
{
  CuspList sorted_slist;
  sorted_slist.push_back(RatQuad(Quad::one, Quad::zero));
  if (Quad::class_number==1)
    return sorted_slist;

  // denom 2 sigmas we construct directly
  Quad w = Quad::w;
  long d = Quad::d;

  switch (d%8) {
  case 1: case 5:
    sorted_slist.push_back(RatQuad(1+w,TWO));
    break;
  case 2: case 6:
    sorted_slist.push_back(RatQuad(w,TWO));
    break;
  case 7:
    sorted_slist.push_back(RatQuad(w,TWO));
    sorted_slist.push_back(RatQuad(1-w,TWO)); // NB not in rectangle: (w-1)/2 is
    break;
  default:
    ;
  }

  // denom 3 sigmas we construct directly
  switch (d%12) {
  case 2: case 5:
    if (d>5)
      {
        sorted_slist.push_back(RatQuad(1+w,THREE));
        sorted_slist.push_back(RatQuad(-1-w,THREE));
        sorted_slist.push_back(RatQuad(1-w,THREE));
        sorted_slist.push_back(RatQuad(-1+w,THREE));
      }
    break;
  case 3:
    if (d>15)
      {
        sorted_slist.push_back(RatQuad(1+w,THREE));
        sorted_slist.push_back(RatQuad(-1-w,THREE)); // NB not in rectangle: (2-w)/3 is
      }
    break;
  case 6: case 9:
    if (d>6)
      {
        sorted_slist.push_back(RatQuad(w,THREE));
        sorted_slist.push_back(RatQuad(-w,THREE));
      }
    break;
  case 11:
    if (d>23)
      {
        sorted_slist.push_back(RatQuad(w,THREE));
        sorted_slist.push_back(RatQuad(-w,THREE));
        sorted_slist.push_back(RatQuad(1-w,THREE));
        sorted_slist.push_back(RatQuad(-1+w,THREE));
      }
    break;
  default:
    ;
  }

  if (verbose)
    cout<<"Sigmas with small denominators: "<<sorted_slist<<endl;

  // Now process the other sigmas (if any)
  Quad temp;
  for ( auto s : slist) // not const or reference as we may change it
    {
      if (verbose)
        cout <<"sigma = "<<s<<endl;
      assert (s.in_rectangle());

      // skip oo and denom 2 or 3 sigmas:
      if (s.is_infinity() or (TWO*s).is_integral() or (THREE*s).is_integral())
        {
          if (verbose)
            cout <<" - skipping (small denominator)"<<endl;
          continue;
        }

      RatQuad ms = reduce_to_rectangle(-s, temp);
      assert (ms.in_rectangle());

      // skip sigmas if we have seen its negative:
      if (std::find(sorted_slist.begin(), sorted_slist.end(), ms) != sorted_slist.end())
        {
          if (verbose)
            cout <<" - skipping (negative seen already)";
          continue;
        }

      if (s.y_coord()<0)
        {
          sorted_slist.push_back(ms);
          sorted_slist.push_back(s);
        }
      else
        {
          sorted_slist.push_back(s);
          sorted_slist.push_back(ms);
        }
    }
  if (verbose)
    cout<<"All sigmas: "<<sorted_slist<<endl;

  return sorted_slist;
}


// Output sorted list of singular points (oo, denom 2, denom 3, larger denoms in +/- pairs)

void output_singular_points(const CuspList& S, int to_file, int to_screen)
{
  vector<Quad> small_denoms = {Quad::zero, Quad::one, 2*Quad::one, 3*Quad::one};
  ofstream geodata;
  stringstream ss;
  if (to_file)
    {
      ss << "geodata_" << Quad::d << ".dat";
      geodata.open(ss.str().c_str(), ios_base::app);
    }
  int nlines=0;
  for ( const auto& s : S)
    {
      Quad sden = s.den(), snum = s.num();
      if ((s.y_coord()>0) && (std::find(small_denoms.begin(), small_denoms.end(), sden) == small_denoms.end()))
        {
          nlines++;
          if (to_file)
            {
              geodata << Quad::d << " S ";
              geodata << snum.re() << " " << snum.im() << " ";
              geodata << sden.re() << " " << sden.im() << endl;
            }
          if (to_screen)
            {
              cout << Quad::d << " S ";
              cout << snum.re() << " " << snum.im() << " ";
              cout << sden.re() << " " << sden.im() << endl;
            }
        }
    }
  if (to_file)
    geodata.close();
  if (to_screen)
    cout << nlines << " S lines output" <<endl;
}

// direct lists of alphas and sigmas of denominator 2 or 3:
CuspList denom_2_alphas()
{
  CuspList alist;
  if (Quad::is_Euclidean) return alist;
  Quad w = Quad::w;
  vector<Quad> numlist;
  int d8 = (Quad::d)%8;
  if (d8==1 || d8==5) numlist = {w};
  if (d8==2 || d8==6) numlist = {1+w};
  if (d8==3) numlist = {w,w-1};
  for (const auto& n : numlist) alist.push_back(RatQuad(n,TWO));
  return alist;
}

CuspList denom_2_sigmas()
{
  CuspList slist;
  if (Quad::class_number == 1) return slist;
  Quad w = Quad::w;
  vector<Quad> numlist;
  int d8 = (Quad::d)%8;
  if (d8==1 || d8==5) numlist = {1+w};
  if (d8==2 || d8==6) numlist = {w};
  if (d8==7) numlist = {w,1-w};
  for (const auto& n : numlist) slist.push_back(RatQuad(n,TWO));
  return slist;
}

CuspList denom_3_alphas()
{
  CuspList alist;
  if (Quad::is_Euclidean) return alist;
  int d = Quad::d, d12 = (Quad::d)%12;
  if (d==5 || d==6 || d==15 || d==19 || d==23) return alist;
  Quad w = Quad::w;
  vector<Quad> numlist;
  switch (d12) {
  case 3:
    numlist = {w, w-1}; break;
  case 7:
    if (d>31)
      numlist = {w, 1-w, 1+w};
    else
      numlist = {1+w};
    break;
  case 11:
    numlist = {w+1}; break;
  case 1: case 10:
    numlist = {w, 1+w, 1-w}; break;
  case 2: case 5:
    numlist = {w}; break;
  case 6: case 9:
    numlist = {w+1, w-1}; break;
  }
  for (const auto& n : numlist)
    {
      alist.push_back(RatQuad(n,THREE));
      alist.push_back(RatQuad(-n,THREE));
    }
  return alist;
}

CuspList denom_3_sigmas()
{
  CuspList slist;
  if (Quad::class_number == 1) return slist;
  int d = Quad::d, d12 = (Quad::d)%12;
  if (d==5 || d==6 || d==23 || d==15) return slist;
  Quad w = Quad::w;
  vector<Quad> numlist;
  switch (d12) {
  case 2: case 5:
    numlist = {1+w, 1-w}; break;
  case 3:
    numlist = {w+1}; break;
  case 6: case 9:
    numlist = {w}; break;
  case 11:
    numlist = {w, w-1}; break;
  }
  for (const auto& n : numlist)
    {
      slist.push_back(RatQuad(n,THREE));
      slist.push_back(RatQuad(-n,THREE));
    }
  return slist;
}

// Square radius for principal cusp
RAT radius_squared(const RatQuad& a)
{
  RAT rsq(a.den().norm());
  return rsq.recip();
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
  RAT r1sq = radius_squared(a1), r2sq = radius_squared(a2);
  RAT d1 = (a1-a2).norm() - (r1sq + r2sq);
  RAT d2 = d1*d1 - 4*r1sq*r2sq;
  return (d2 < 0? 0: d1.sign() * (d2==0? 1: 2));
}

// return 1 iff the circle S_A1 is inside S_a2
int circle_inside_circle(const RatQuad& a1, const RatQuad& a2, int strict)
{
  RAT r1sq = radius_squared(a1), r2sq = radius_squared(a2);
  if (! (strict? r1sq<r2sq: r1sq<=r2sq))
    {
      // cout<<"is "<<a1<<" inside "<<a2<<"? no (radius test)\n";
      return 0;
    }
  int t  = tau(a1,a2);
  int ans = (strict? t==-2: t<0);
  // cout<<"is "<<a1<<" inside "<<a2<<"? no (tau="<<t<<")\n";
  return ans;
}

// return 1 iff the circle S_a is inside S_b for any b in blist
int circle_inside_any_circle(const RatQuad& a, const CuspList& blist, int strict)
{
  return std::any_of(blist.begin(), blist.end(),
                     [a, strict](RatQuad b) {return circle_inside_circle(a,b,strict);});
}

// Return list of 0, 1 or 2 sqrts of a rational r in k
vector<RatQuad> sqrts(const RAT& r)
{
  int s = r.sign();
  if (s>0)
    return vector<RatQuad>(); // empty
  if (s==0)
    return {RatQuad()};       // 0
  RAT root;
  if (r.is_square(root))      // rational square
    {
      RatQuad rt(root);
      assert (rt*rt==RatQuad(r));
      return {rt, -rt};
    }
  RAT rd = -r*INT(Quad::d);
  if (rd.is_square(root))      // -d * rational square
    {
      // now root^2 = -d*r, so (rt/sqrt(-d))^2 = r
      Quad root_minus_d = (Quad::t ? Quad(-1,2) : Quad::w);
      RatQuad rt(root);
      rt /= root_minus_d;
      assert (rt*rt==RatQuad(r));
      return {rt, -rt};
    }
  else
    {
      return vector<RatQuad>(); // empty
    }
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
  RatQuad z = a1 + a2 + (r1sq-r2sq)/delta;
  // find sqrts of d2 in k, if any
  vector<RatQuad> d2sqrts = sqrts(d2);
  for ( const auto& sqrtd2 : d2sqrts)
    {
      RatQuad pt = z/TWO + sqrtd2/(TWO*delta);
      assert ((pt-a1).norm() == r1sq);
      assert ((pt-a2).norm() == r2sq);
      ans.push_back(pt);
    }
  return ans;
}


// return 1 iff a is [strictly] inside S_b
int is_inside(const RatQuad& a, const RatQuad& b, int strict)
{
  RAT t = radius_squared(b) - (a-b).norm();
  if (0)
    {
      cout<<"Testing whether "<<a<<" is strictly inside S_{"<<b<<"}\n";
      cout<<" square radius = "<<radius_squared(b);
      cout<<" square distance = "<<(a-b).norm();
      cout<<" result: "<<(t>0)<<endl;
    }
  return (strict ? 0<t : 0<=t);
}

// return 1 iff a is [strictly] inside S_b for at least one b in blist
int is_inside_one(const RatQuad& a, const CuspList& blist, int strict)
{
  return std::any_of(blist.begin(), blist.end(),
                     [a, strict](RatQuad b) {return is_inside(a,b,strict);});
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
  return principal_cusps_with_denominators(quads_of_norm_up_to(maxn));
}

// list of principal cusps with given denominator
CuspList principal_cusps_with_denominator(const Quad& s)
{
  CuspList alist;
  Quad temp;
  auto rlist = invertible_residues(s);
  for ( const auto& r : rlist)
    {
      auto a = reduce_to_rectangle(RatQuad(r, s), temp);
      if (!circle_inside_any_circle(a, alist))
        alist.push_back(a);
    }
  return alist;
}

// list of principal cusps with denominator in given list
CuspList principal_cusps_with_denominators(const vector<Quad>& slist)
{
  CuspList alist;
  for ( const auto& s : slist)
    {
      auto alist1 = principal_cusps_with_denominator(s);
      alist.insert(alist.end(), alist1.begin(), alist1.end());
    }
  return alist;
}

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
  RAT
    rho0(radius_squared(a0)),
    rho1(radius_squared(a1)),
    rho2(radius_squared(a2));
  RatQuad delta = tri_det(a0,a1,a2);
  if (delta.is_zero())
    return points; // empty
  delta = delta.conj(); // to match Sage code
  RAT
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
  assert (is_under(P,a0)==0);
  assert (is_under(P,a1)==0);
  assert (is_under(P,a2)==0);
  if (t2>=0)
    points.push_back(P);
  return points;
}

// Given principal cusps a1, a2, a such that the circles S_a1 and
// S_a2 intersect in distinct points, test whether S_a covers either
// or both these points.

// Returns 0 if neither, 2 if both, +1 or -1 if just one.  The signs
// are consistent so that if a returns +1 and a' returns -1 then each
// intersection point is covered by either S_a or S_a'.

//#define DEBUG_ARE_INTERSECTION_POINTS_COVERED_BY_ONE

int are_intersection_points_covered_by_one(const RatQuad& a1, const RatQuad& a2, const RatQuad& a, int debug)
{
#ifdef DEBUG_ARE_INTERSECTION_POINTS_COVERED_BY_ONE
  debug=1;
#endif
  if (debug)
    cout << "\tTesting if intersection points of "<<a1<<" and "<<a2<<" are covered by "<<a;

  // Check the cusps are principal, not infinity, and with unit ideal
  assert (a.is_finite() && a.is_principal());
  assert (a1.is_finite() && a1.is_principal());
  assert (a2.is_finite() && a2.is_principal());

  // Define the square radii and centres
  RatQuad delta = (a2-a1).conj();
  RAT
    r1sq = radius_squared(a1),
    r2sq = radius_squared(a2),
    rsq = radius_squared(a),
    n1 = a1.norm(),
    n2 = a2.norm(),
    n = delta.norm(),
    d1 = n - (r1sq + r2sq),
    d2 = d1*d1 - 4*r1sq*r2sq;
  assert (d2 < 0);

  RatQuad z0 = ((a1+a2) + (r1sq-r2sq)/delta)/TWO;
  RAT T = 2 * n * (rsq - (z0-a).norm()) + d2/TWO;
  RAT T2 = T*T;
  RatQuad D = tri_det(a, a2, a1); // pure imaginary
  RAT D2 = (D*D).x_coord(1);      // negative rational
  RAT d2D2 = d2*D2;               // positive rational

  // the covering condition is \pm sqrt(d2)*D < T

  int code = 0;
  if (d2D2 < T2)
    {
      code = (T>0 ? 2 : (T<0? 0: 99));
    }
  if (d2D2 > T2)
    {
      Quad w = Quad::w;
      RAT u = (D*(w-w.conj())).x_coord(1);
      code = u.sign();
    }
  if (debug)
    cout << " - test returns code "<<code<<endl;
  return code;
}

// Given principal cusps a0, a1 whose circles S_a0, S_a1 intersect,
// and a list of principal cusps alist each of whose circles S_a also
// intersects S_a0, test whether each of the two intersection points
// of S_a0 and S_a1 is either singular or strictly inside one of the
// S_a.  We treat as a special case when the two intersection points
// are in k.  If not, the code still uses exact arithmetic.

int are_intersection_points_covered(const RatQuad& a0, const RatQuad& a1, const CuspList& alist,
                                    const CuspList& sigmas, int debug)
{
#ifdef DEBUG_ARE_INTERSECTION_POINTS_COVERED_BY_ONE
  debug = 1;
#endif
  if (debug)
    cout << "Testing if intersection points of "<<a0<<" and "<<a1<<" are covered"<<endl;// by "<<alist<<endl;

  CuspList zlist = intersection_points_in_k(a0,a1);
  if (!zlist.empty())
    {
      if (debug)
        cout << " intersection points are k-rational: "<<zlist<<endl;
      for ( const auto& z : zlist)
        {
          Quad t;
          if (cusp_index_with_translation(z, sigmas, t)!=-1) // z is singular: OK
            {
              if (debug)
                cout << " ok, "<<z<<" is singular"<<endl;
              continue;
            }
          if (!is_inside_one(z, alist, 1))  // z is not covered: not OK
            {
              if (debug)
                cout << " returning no, "<<z<<" is not covered"<<endl;
              return 0;
            }
        }
      // we reach here if every z is either singular or covered: OK
      if (debug)
        cout << "+++returning yes, both are covered"<<endl;
      return 1;
    }

  // Now the intersection points are not in k. Check that either one
  // S_a covers both, or two cover one each:

  int t = 0; // will hold +1 or -1 if we have covered only one of the two
  for ( const auto& a2 : alist)
    {
      if ((a2 == a1) || (a2==a0))
        continue;
      if (tau(a0,a2)>0 || tau(a1,a2)>0)
        continue;
      int t2 = are_intersection_points_covered_by_one(a0, a1, a2);
      if (t2==2) // both are covered by a2
        {
      if (debug)
        cout << "+++returning yes, both are covered"<<endl;
      return 1;
        }
      if (t2==0) // neither is covered by a2
        continue;
      // Now t2 is +1 or -1; we win if t is its negative
      if (t==-t2) // then they are (1,-1) or (-1,1) so we have covered both points
        {
          if (debug)
            cout << " returning yes, both are now covered"<<endl;
          return 1;
        }
      if (debug)
        cout << "   not yet, only one is covered so far"<<endl;
      t = t2;     // = +1 or -1: we have covered one of the points, so remember which
    }
  // If we reach here then none of the S_a covers both points
  if (debug)
    cout << "---returning no"<<endl;
  return 0;
}

// Given a principal cusp a0, a candidate list of principal cusps
// alist, tests whether the boundary of the disc S_a0 is contained in
// the union of the S_b for b in alist, apart from any singular
// points on the boundary.  It suffices to consider all b such that
// S_b intersects S_a in two points and check that each of the
// points is either singular or contained in some other S_b.  This
// is simplest when the intersection points are in k; if not then the
// method still uses exact arithmetic in k throughout.

// pairs_ok is list of pairs {a1,a2} whose intersection points are
// known to be covered, which will be added to.

// Returns 0 or 1, and pairs_ok is updated to a list of pairs
// whose intersections have now been shown to be covered.

int is_alpha_surrounded(const RatQuad& a0, const CuspList& alist, const CuspList& sigmas,
                        vector<CuspPair>& pairs_ok, int debug)
{
  if (debug)
    cout<<"Testing if S_{"<<a0<<"} is surrounded"<<endl;
  assert (a0.in_rectangle());
  // check if S_a0 is strictly entirely contained in one S_alpha:
  if (circle_inside_any_circle(a0, alist, 1))
    {
      if (debug)
        cout<<" - yes, S_{"<<a0<<"} is contained within another S_a"<<endl;
      return 1;
    }
  // extract the relevant alphas, if any, namely those for which
  // S_alpha and S_a0 properly intersect:

  CuspList a0list(alist.size()); // upper bound on size
  auto it1 = std::copy_if(alist.begin(), alist.end(), a0list.begin(),
                         [a0](RatQuad a){return circles_intersect(a0, a);});
  a0list.resize(std::distance(a0list.begin(),it1));  // shrink to new size
  if (debug)
    cout<<" - intersecting neighbours: "<<a0list<<endl;

  // extract the pairs in pairs_ok which contain a0:

  vector<CuspPair> a0_pairs_ok(pairs_ok.size());
  auto it2 = std::copy_if(pairs_ok.begin(), pairs_ok.end(), a0_pairs_ok.begin(),
                         [a0](CuspPair pr){return (pr.find(a0)!=pr.end());});
  a0_pairs_ok.resize(std::distance(a0_pairs_ok.begin(),it2));  // shrink to new size

  int all_ok = 1;
  for ( const auto& a1 : a0list)
    {
      CuspPair pr = {a0,a1};
      if (std::find(a0_pairs_ok.begin(), a0_pairs_ok.end(), pr) != a0_pairs_ok.end())
        continue; // this pair is already ok
      if (are_intersection_points_covered(a0, a1, alist, sigmas))
        {
          pairs_ok.push_back(pr); // record that this pair is ok
          if (debug)
            cout<<" - intersection points of S_{"<<a0<<"} and S_{"<<a1<<"} are surrounded"<<endl;
        }
      else
        {
          if (debug)
            cout<<" - intersection points of S_{"<<a0<<"} and S_{"<<a1<<"} are NOT surrounded"<<endl;
          all_ok = 0; // but continue checking the other a1s
        }
    }
  if (debug)
    {
      if (all_ok)
        cout << " S_{"<<a0<<"} IS surrounded"<<endl;
      else
        cout << " S_{"<<a0<<"} is NOT surrounded"<<endl;
    }
  return all_ok;
}

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

// We assume that all a in alist_open are in the quarter rectangle

// Returns 0/1

// NB All a in alist_open will be tested, i.e. we carry on after a
// failure.

int are_alphas_surrounded(CuspList& alist_ok, CuspList& alist_open,
                          const CuspList& slist, vector<CuspPair>& pairs_ok,
                          int verbose, int debug)
{
  if (debug)
    {
      cout << "At start of are_alphas_surrounded()" << endl;
      cout << "alist_ok = "<<alist_ok<<endl;
      cout << "alist_open = "<<alist_open<<endl;
    }
  // concatenate the two alists
  CuspList alist;
  alist.insert(alist.end(), alist_ok.begin(), alist_ok.end());
  alist.insert(alist.end(), alist_open.begin(), alist_open.end());

  // add translates
  CuspList alistx;
  for ( auto& a : alist)
    for (auto x : {-1,0,1})
      for (auto y : {-1,0,1})
        {
          alistx.push_back(a+Quad(x,y));
        }
  int i=0, n_open = alist_open.size();

  // We make a copy of alist_open to loop over, so we can delete elements from the original as we go
  CuspList new_alist_open = alist_open;
  for ( const auto& a : new_alist_open)
    {
      assert (a.in_quarter_rectangle());
      i++;
      if (verbose) cout <<"Testing alpha #"<<i<<"/"<<n_open<<" = "<<a<<"...";
      if (is_alpha_surrounded(a, alistx, slist, pairs_ok))
        {
          if (verbose) cout << " ok! surrounded\n";
          // add this alpha to the ok list end remove from the open list
          alist_ok.push_back(a);
          alist_open.erase(std::find(alist_open.begin(), alist_open.end(), a));
        }
      else
        {
          if (verbose) cout << " no, not surrounded" << endl;
          return 0;
        }
    }

  // Test that the singular points are surrounded:
  int ok = are_sigmas_surrounded(slist, alist, debug);
  if (verbose)
    {
      if (ok)
        cout << " and singular points are surrounded!" <<endl;
      else
        cout << " but singular points are not yet surrounded" <<endl;
    }
  return ok;
}

// Returns a finite list of principal cusps a such that the S_{a+t}
// for all integral t cover CC apart from singular points.

// For n>=1 successively, we test as a candidate set all a=r/s with
// r,s coprime, r reduced mod s, N(s)<=n (omitting any for which S_a
// is contained in any earlier S_a') until we succeed.

// sigmas can be set to a list of singular points (up to translation),
// otherwise these will be computed.

// Other functions will then (1) saturate the set, (2) discard redundancies.

CuspList covering_alphas(const CuspList& sigmas, int verbose)
{
  CuspList alphas_ok;  // list that will be returned
  INT maxn = Quad::absdisc/4; // we first consider all alphas with dnorm up to this
  if (maxn==0) maxn=1; // just for discriminant -3
  int first = 1;
  string s = "of norm up to ";
  CuspList alist, alphas_open;
  vector<CuspPair>pairs_ok;
  Quadlooper looper(1, 0, 1);  // norms from 1 to oo, both conjugates (up to units)
  while (1)
    {
      // Get the next batch new_alphas, either of all dnorms up to maxn, or those with the next dnorm

      CuspList new_alphas = (first
                             ?
                             principal_cusps_with_denominators(looper.values_with_norm_up_to(maxn))
                             :
                             principal_cusps_with_denominators(looper.values_with_current_norm())
                             );
      //cout << "new_alphas = " << new_alphas << endl;
      maxn = new_alphas.back().den().norm();

      if (verbose)
        cout << "----------------------------------------------------------\n"
             << "Considering "<<new_alphas.size()<<" extra principal cusps " << s << maxn
          //<< ": "<<new_alphas
             << endl;
      first = 0;
      s = "of norm ";

      // Of these, check for redundancy; if not redundant, put into
      // alist and also either into new_alphas_open (if in quarter
      // rectangle) or into alphas_ok.  We only use new_alphas_open
      // for ease of reporting output.

      CuspList new_alphas_open;

      for (auto a : new_alphas)
        {
          assert (a.in_rectangle());
          if (circle_inside_any_circle(a, alist, 0)) // strict=0
            {
              //cout<<"alpha = "<<a<<" is weakly inside one of "<<alist<<endl;
              continue;
            }
          // cout<<"alpha = "<<a<<" is useful and "
          //     <<(a.in_quarter_rectangle()?"is":"is not")<<" in the quarter rectangle"<<endl;
          alist.push_back(a);
          if (a.in_quarter_rectangle())
            new_alphas_open.push_back(a);
          else
            alphas_ok.push_back(a);
        }
      int nc = new_alphas_open.size();

      // report on which new alphas we will be testing, if any:

      if (verbose)
        {
          if (nc)
            {
              cout << "Adding "<<nc<<" alphas " << s << maxn << ": " << new_alphas_open;
              cout << " (plus symmetrics); ";
              cout << "#alphas="<<alphas_ok.size()+alphas_open.size()+new_alphas_open.size();
              // cout << " of which "<<alphas_ok.size()<<" are proved surrounded so far";
              cout << endl;
            }
          else
            {
              cout << "All are redundant\n";
            }
        }
      if (nc==0)
        continue;

      // Append new_alphas_open to alphas_open:

      alphas_open.insert(alphas_open.end(), new_alphas_open.begin(), new_alphas_open.end());

      // Test whether all alphas_open are now surrounded (by
      // translates of alphas_ok+alphas_open).  As a side effect, some
      // alphas will be moved from alphas_open to alphas_ok, and
      // pairs_ok (which holds a list of pairs of alphas whose
      // intersections are known to be covered) will be updated.

      if (are_alphas_surrounded(alphas_ok, alphas_open, sigmas, pairs_ok, verbose, verbose>1))
        {
          if (verbose)
            cout << "Success using "<<alphas_ok.size()<<" alphas of with max norm "<<maxn<<"\n";
          return alphas_ok;
        }
      else
        {
          if (verbose)
            cout << "Some alphas are not surrounded, continuing...\n";
        }
    }
}

// Return the height of S_a above P, or 0 if S_a does not cover P
RAT height_above(const RatQuad& a, const RatQuad& z)
{
  return radius_squared(a)-(z-a).norm(); // positive iff S_a covers z
}

// return -1,0,+1 according as P is over, on, under S_a (a principal)
int is_under(const H3point& P, const RatQuad& a)
{
  return sign(height_above(a, P.first) - P.second);
}

// return +1 iff P is under at least one S_a for a in sliat
int is_under_any(const H3point& P, const CuspList& alist)
{
  return std::any_of(alist.begin(), alist.end(),
                     [P](RatQuad a) {return is_under(P,a)==1;});
}

// Given a list of principal cusps alpha (all reduced mod O_k) return
// a list of "corners" P = [z,tsq] each the intersection of an S_a
// with at least two other S_{b+t} with z in the fundamental
// rectangle and tsq>0.

// Let u = (w-wbar)/2.  The fundamental rectangle F has TR corner at
// (u+1)/2 and BL corner minus this.  Using symmetries (negation and
// conjugation) we can work with the quarter-rectangle F4 with the
// same TR and BL=0.  To recover F from F4 take the union of
// z,-z,zbar,-zbar for z in F4.

// The 9 quarter-rectangles adjacent to F4 consist of
//  -z, zbar, 1-zbar; -zbar, z, 1-zbar; w-z, w+zbar, and either w+1-z or w+zbar-1
// for z in F4.

CuspList nbrs(const RatQuad& z)
{
  RatQuad cz = z.conj();
  int t = Quad::t;
  Quad w = Quad::w;
  return {z,-z, cz, ONE-z, -cz, ONE-cz, w-z, w+cz, (t? w+cz-ONE : w+ONE-z)};
}

ostream& operator<<(ostream& s, const H3point& P)
{
  s << "[" << P.first<<","<<P.second<<"]";
  return s;
}

H3pointList triple_intersections(const CuspList& alphas, int debug)
{
    if (debug)
      cout << "Finding triple intersections..."<<endl;

    // Extract the alphas in F4:
    CuspList alphasF4(alphas.size());
    auto it1 = std::copy_if(alphas.begin(), alphas.end(), alphasF4.begin(),
                            [](RatQuad a){return a.in_quarter_rectangle();});
    alphasF4.resize(std::distance(alphasF4.begin(),it1));  // shrink to new size
    if (debug)
      cout << alphasF4.size() <<" alphas are in the quarter rectangle" << endl;

    // Extend these by 8 translations:
    CuspList alphasF4X;
    for ( auto& z : alphasF4)
      {
        CuspList z_nbrs = nbrs(z);
        alphasF4X.insert(alphasF4X.end(), z_nbrs.begin(), z_nbrs.end());
      }
    if (debug)
      cout << alphasF4X.size() <<" neighbours of these" << endl;

    // convert each cusp z to a point P = [z,tsq] with tsq the square
    // radius of S_z:

    // from utils import frac
    // Alist = [[a,1/frac(a)[1].norm()] for a in XA4]

    // For each i get a list of j>i for which S_ai and S_aj intersect properly
    map <int, vector<int> > i2j;
    int i=0, n=alphasF4X.size();
    for ( const auto& a : alphasF4X)
      {
        vector<int> i2j_i;
        for (int j = i+1; j<n; j++)
          if (circles_intersect(a, alphasF4X[j]))
            i2j_i.push_back(j);
        i2j[i] = i2j_i;
        i++;
      }
    if (debug)
      cout << " finished making i2j" <<endl;

    // Hence make a list of triples (i,j,k) with i<j<k with pairwise proper intersections

    vector<vector<int>> ijk_list;
    for (const auto& i_j_list : i2j)
      {
        int i = i_j_list.first;
        vector<int> j_list = i_j_list.second;
        for (const auto& j : j_list)
          for (const auto& k : j_list)
            if (std::find(i2j[j].begin(), i2j[j].end(), k) != i2j[j].end()) // k is in i2j[j]
              ijk_list.push_back({i,j,k});
      }

    if (debug)
      cout <<" finished making ijk_list: "<<ijk_list.size()<<" triples\n";

    H3pointList points;

    for ( const auto& ijk : ijk_list)
      {
        int i=ijk[0], j=ijk[1], k=ijk[2];
        H3pointList points1 = tri_inter_points(alphasF4X[i], alphasF4X[j], alphasF4X[k]);
        if (points1.empty())
          continue;
        H3point P = *points1.begin();
        RatQuad z = P.first;
        RAT t2 = P.second;
        if (t2.sign()==0)
          continue;
        if (debug>2)
          cout << " found P = "<<P<<"\n";
        if (!z.in_quarter_rectangle())
          continue;
        if (std::find(points.begin(), points.end(), P) != points.end())
          continue;
        if (is_under_any(P, alphasF4X))
          continue;
        // These corners are in F4, so we apply symmetries to get all those in F:
        RatQuad zbar = z.conj();
        for (const auto& z2 : {z, -z, zbar, -zbar})
          {
            if (!z2.in_rectangle())
              continue;
            H3point P2 = {z2, t2};
            if (std::find(points.begin(), points.end(), P2) != points.end())
              continue;
            if (debug)
              cout << " adding P2 = "<<P2<<" from (i,j,k) = "<< ijk <<endl;
            points.push_back(P2);
          }
      }
    if (debug)
      cout << " returning "<<points.size() <<" corners" <<endl;
    return points;
}


// count how many P in points are on S_a
int nverts(const RatQuad& a, const H3pointList& points)
{
  return std::count_if(points.begin(), points.end(),
                       [a](H3point P) {return is_under(P,a)==0;});
}


// return sublist of a in alist which have t least 3 vertices in points
CuspList remove_redundants(const CuspList& alist, const H3pointList& points)
{
  CuspList new_alist;
  std::copy_if(alist.begin(), alist.end(), std::back_inserter(new_alist),
               [points](RatQuad a) {return nverts(a, points) >= 3;});
  return new_alist;
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
        else cout<<" (feault, i.e. P is on or under S_a)";
      cout <<endl;
    }
  CuspList ans;
  RatQuad z = P.first, sz;
  RAT t2 = P.second;
  Quad r, temp;
  int ok, test;
  long norm_s_ub = I2long((1/t2).floor());
  if (debug)
    cout << "t2 = "<<t2<<" so bounds on N(s) are ["<<norm_s_lb<<","<<norm_s_ub<<"]\n";

  // We could do Quadlooper sloop(norm_s_lb, norm_s_ub, 1); but it's
  // more efficient to construct all possible s at once, though that
  // takes more memory
  auto slist = quads_of_norm_between(norm_s_lb, norm_s_ub, 1, 0); // including conjugates, not sorted
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

H3point translate(const H3point& P, const Quad& t)
{
  return {P.first + t, P.second};
}

// return max denomaintor norm of a list of principal cusps
INT max_dnorm(const CuspList& alphas)
{
  INT m;
  std::for_each(alphas.begin(), alphas.end(),
                [&m](RatQuad a) {INT n = a.den().norm(); if (n>m) m=n;});
  return m;
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
                  [P](RatQuad a) {cout<<"a="<<a<<": height above P is "<<height_above(a,P.first)<<endl;});
  std::for_each(alist.begin(), alist.end(),
                [P,&m](RatQuad a) {m = max(m, height_above(a,P.first));});
  if(debug)
    cout<<"max height is "<<m<<endl;
  // Discard those a whose height above P is not maximal
  alist.erase(std::remove_if(alist.begin(), alist.end(),
                              [P,m](RatQuad a) { return height_above(a,P.first) < m;}),
               alist.end());
  if(debug)
    {
      cout<<"S_a maximally covering P="<<P<<": "<<alist<<" (ht "<<m<<" above P)"<<endl;
      // for (const auto& a: alist)
      //   cout<<a<<" with height "<<height_above(a,P.first)<<endl;
    }
  return alist;
}

// Given a covering set of alphas as produced by
// find_covering_alphas(), add extras if necessary so that they are
// "saturated", i.e. define the extended fundamental domain.

// By Swan, we need to find the points P in H^3 with positive height
// where at least 3 hemispheres S_a intersect, and for each P check
// whether P is properly covered by an S_a for a not in the set of
// alphas (up to translation).  If so, we need to add a to the set of
// alphas.  If none, then we have the fundamental region (and can go
// on to discard any redundant alphas).

// We assume that we have already considered all alpha=r/s with N(s)<=maxn.

// At the end we discard any alphas with <3 vertices (including
// translates and singular points), and return the new set of alphas

CuspList saturate_covering_alphas(const CuspList& alphas, const CuspList& sigmas, INT maxn, int debug, int verbose)
{
  INT m;
  if (verbose)
    {
      m = max_dnorm(alphas);
      cout << "Saturating "<<alphas.size()<<" alphas with max dnorm "<< m <<endl;
    }
  // First delete any alphas with <3 vertices, allowing for translates
  H3pointList points = triple_intersections(alphas, debug);
  if (verbose)
    cout << "Found "<<points.size() << " triple intersection points" <<endl;

  // add translates of these and singular points
  vector<Quad> translates = { Quad::one, Quad::w, Quad::one+Quad::w, Quad::one-Quad::w };
  H3pointList pointsx;
  for ( const auto& P : points)
    {
      pointsx.push_back(P);
      for ( const auto& t : translates)
        {
          pointsx.push_back(translate(P,t));
          pointsx.push_back(translate(P,-t));
        }
    }
  for ( const auto& s : sigmas)
    {
      H3point P = {s, ZERO};
      pointsx.push_back(P);
      for ( const auto& t : translates)
        {
          pointsx.push_back(translate(P,t));
          pointsx.push_back(translate(P,-t));
        }
    }
  CuspList new_alphas = remove_redundants(alphas, pointsx);
  m = max_dnorm(new_alphas);
  if (verbose)
    cout << "After removing alphas which go through <3 vertices, we now have "
         <<new_alphas.size()<<" alphas with max norm "<< m <<endl;

  int sat = 0, first_run=1;
  H3pointList checked_points;
  while (!sat)
    {
      if (!first_run)
        points = triple_intersections(new_alphas, debug);
      first_run = 0;
      if (verbose)
        {
          cout << "Found "<<points.size()<<" potential vertices"<<endl;
          // cout << points <<endl;
        }
      if(debug)
        cout << "Extracting those of square height less than 1/"<<maxn<<endl;
      // Remove points which cannot be better covered by an alpha with dnorm>maxn
      points.erase(std::remove_if(points.begin(), points.end(),
                                  [maxn](H3point P) { return maxn*P.second>=ONE;}),
                   points.end());
      if (debug)
        {
          cout << " -- of which "<<points.size()<<" are low enough to be properly covered by a new alpha"<<endl;
          // for (const auto& P : points)
          //   cout << "P = " << P << " is in first quadrant? "<< P.first.in_quarter_rectangle() << endl;
        }
      points.erase(std::remove_if(points.begin(), points.end(),
                                  [](H3point P) { return !P.first.in_quarter_rectangle();}),
                   points.end());
      if (debug)
        {
          cout << " -- of which "<<points.size()<<" lie in the first quadrant" << endl;
          // cout << points <<endl;
        }

      points.erase(std::remove_if(points.begin(), points.end(),
                                  [checked_points](H3point P)
                                  {return std::find(checked_points.begin(), checked_points.end(), P) != checked_points.end();}),
                   points.end());
      if (debug)
        {
          cout << " -- of which "<<points.size()<<" have not already been checked" << endl;
          // cout << points <<endl;
        }

      sat = 1;        // will be set to 0 if we find out that the alphas are not already saturated
      CuspList extra_alphas; // will be filled with any extra alphas needed on this pass
      int iP = 0, nP = points.size();
      for (const auto& P : points)
        {
          iP++;
          if (debug)
            cout << " - checking corner #"<<iP<<"/"<<nP<<": "<<P<<"...";
          CuspList extras = best_covering_hemispheres(P, I2long(maxn)+1, debug>1);
          if (debug)
            cout << " done...";
          if (extras.empty())
            {
              if (debug)
                cout << "   - no properly covering alphas found" <<endl;
              checked_points.push_back(P);
              continue; // on to the next P
            }
          sat = 0; // we now know the alphas are not saturated
          if (debug)
            {
              cout << "   - found " <<extras.size()<<" maximally properly covering hemispheres";
              vector<INT> norms;
              norms.resize(extras.size());
              std::transform(extras.begin(), extras.end(), norms.begin(),
                             [](RatQuad a) {return a.den().norm();});
              cout << " with denominator norms " << norms;
              cout << ";  heights above P (height "<<P.second<<") are ";
              vector<RAT> hts;
              hts.resize(extras.size());
              std::transform(extras.begin(), extras.end(), hts.begin(),
                             [P](RatQuad a) {return height_above(a,P.first);});
              cout <<hts <<endl;
            }
          for ( auto& a : extras)
            {
              Quad temp;
              RatQuad ca = a.conj();
              CuspList blist = {a,-a,ca,-ca};
              for ( auto& b : blist)
                {
                  b = reduce_to_rectangle(b, temp);
                  if (debug)
                    cout << " - testing potential new alpha " << b << endl;
                  if (b.in_rectangle() &&
                      std::find(extra_alphas.begin(), extra_alphas.end(), b) == extra_alphas.end())
                    {
                      if (debug)
                        cout << " - adding alpha " << b << endl;
                      extra_alphas.push_back(b);
                    }
                  else
                      if (debug)
                        {
                          cout << " - not adding it: ";
                          if (b.in_rectangle())
                            cout << "it's a repeat" <<endl;
                          else
                            cout << "it's not in the rectangle" <<endl;
                        }
                } // loop over 4 flips of a
            } // loop over a in extras
        } // checked all P in points
      if (verbose)
        {
          if (sat)
            {
              m = max_dnorm(new_alphas);
              cout << " alphas are saturated! "<<new_alphas.size()<<" alphas with max norm "<<m<<endl;
            }
          else
            {
              m = max_dnorm(extra_alphas);
              cout << " alphas not saturated, "<<extra_alphas.size()<<" extras needed: "<<extra_alphas<<" (with norms at most "<<m<<")"<<endl;
            }
        }
      new_alphas.insert(new_alphas.end(), extra_alphas.begin(), extra_alphas.end());
    } // ends while(!sat)

  m = max_dnorm(new_alphas);
  if (verbose)
    cout << "After saturation we now have "<<new_alphas.size()<<" alphas with max norm "<<m<<endl;

  // Now again delete any alphas with <3 vertices, allowing for translates
  points = triple_intersections(new_alphas, debug);
  pointsx.clear();
  for ( const auto& P : points)
    {
      pointsx.push_back(P);
      for ( const auto& t : translates)
        {
          pointsx.push_back(translate(P,t));
          pointsx.push_back(translate(P,-t));
        }
    }
  for ( const auto& s : sigmas)
    {
      H3point P = {s, ZERO};
      pointsx.push_back(P);
      for ( const auto& t : translates)
        {
          pointsx.push_back(translate(P,t));
          pointsx.push_back(translate(P,-t));
        }
    }
  new_alphas = remove_redundants(new_alphas, pointsx);
  m = max_dnorm(new_alphas);
  if (verbose)
    cout << "After removing alphas which go through <3 vertices, we now have "
         <<new_alphas.size()<<" alphas with max norm "<< m <<endl;
  std::sort(new_alphas.begin(), new_alphas.end(), Cusp_cmp);
  return new_alphas;
}

// return  a saturated irredundant list of alphas in the fundamental rectangle
CuspList find_alphas(const CuspList& sigmas, int debug, int verbose)
{
  auto alphas = covering_alphas(sigmas, verbose);
  INT maxn = max_dnorm(alphas);
  return saturate_covering_alphas(alphas, sigmas, maxn, debug, verbose);
}

// return  a saturated irredundant list of alphas, and list of sigmas, in the fundamental rectangle
pair<CuspList,CuspList> find_alphas_and_sigmas(int debug, int verbose)
{
  auto sigmas = singular_points();
  auto alphas = find_alphas(sigmas, debug, verbose);
  return {alphas, sigmas};
}


// test whether angle between s-->a1 and s-->a2 is <180 degrees
int angle_under_pi(const RatQuad& s, const RatQuad& a1, const RatQuad& a2)
{
  return ((s-a1)*(s-a2).conj()).y_coord() > 0;
}

// return list of alphas (or translates) which pass through a finite singular point
CuspList neighbours(const RatQuad& s, const CuspList& alphas)
{
  CuspList nbrs;
  for ( const auto& a : alphas)
    {
      RAT rsq = radius_squared(a); // same for integral translates
      for (auto x : {-1,0,1})
        for (auto y : {-1,0,1})
          {
            RatQuad b = a + Quad(x,y);
            if ((s-b).norm()==rsq)
              nbrs.push_back(b);
          }
    }
  return nbrs;
}

// test if one singular point (sigma) is surrounded by alpha circles:
int is_sigma_surrounded(const RatQuad& sigma, const CuspList& alphas, int debug)
{
  if (sigma.is_infinity())
    return 1;
  CuspList nbrs = neighbours(sigma, alphas);
  if (debug)
    cout<<"Neighbours of "<<sigma<<" are "<<nbrs<<endl;

  int ok = std::all_of(nbrs.begin(), nbrs.end(),
                     [sigma, nbrs](const RatQuad& a1)
                     {return std::any_of(nbrs.begin(), nbrs.end(),
                                         [sigma,a1](const RatQuad& a2) {return angle_under_pi(sigma,a1,a2);});});
  if (debug)
    {
      if (ok)
        cout<<" + "<<sigma<<" IS surrounded"<<endl;
      else
        cout<<" - "<<sigma<<" is NOT surrounded"<<endl;
    }
  return ok;
}

// test if all singular points (sigmas) are surrounded by alpha circles:
int are_sigmas_surrounded(const CuspList& sigmas, const CuspList& alphas, int debug)
{
  return std::all_of(sigmas.begin(), sigmas.end(),
                     [alphas,debug](const RatQuad& s) {return is_sigma_surrounded(s, alphas, debug);});
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

// reduce r mod s so that r/s is in the rectangle
Quad rectify(const Quad& r, const Quad& s)
{
  Quad q;
  RatQuad a = reduce_to_rectangle(RatQuad(r,s), q);
  return r - q*s;
}

// Return sorted list of (saturated covering) alphas (denom 1, denom 2, denom 3, larger denoms in pairs or fours)
// + pluspairs, minuspairs, fours
CuspList sort_alphas(const CuspList& A,
                     vector<vector<Quad>>& pluspairs, vector<vector<Quad>>& minuspairs, vector<vector<Quad>>& fours,
                     int verbose, int debug)
{
  CuspList alist, a1list, a2list, a3list;
  map<Quad, vector<Quad>, Quad_comparison> alist_by_denom; // list of numerators for each denominator
  for (const auto& a : A)
    {
      assert (a.in_rectangle());
      Quad r = a.num(), s = a.den();
      if (div(s,r)) // a is integral
        {
          a1list.push_back(a);
          continue;
        }
      if (div(s,2*r)) // 2a is integral
        {
          a2list.push_back(a);
          continue;
        }
      if (div(s,3*r)) // 3a is integral
        {
          a3list.push_back(a);
          continue;
        }
      // Now s (where a=r/s) is not 1,2,3; add s:r to alist_by_denom:
      if (alist_by_denom.find(s) == alist_by_denom.end())
        alist_by_denom[s] = {r};
      else
        alist_by_denom[s].push_back(r);
    }
  if (debug)
    {
      cout << " alphas with denominator 1: " <<a1list <<endl;
      cout << " alphas with denominator 2: " <<a2list <<endl;
      cout << " alphas with denominator 3: " <<a3list <<endl;
    }
  assert (a1list.size()==1 && a1list[0]==RatQuad(0));
  assert (compare_CuspLists_as_sets_mod_translation(a2list, denom_2_alphas()));
  if (!compare_CuspLists_as_sets_mod_translation(a3list, denom_3_alphas()))
    {
      cout <<"!!!\n";
      cout << "alphas with denom 3 (from full list):    " << a3list <<endl;
      cout << "alphas with denom 3 (from special list): " << denom_3_alphas() <<endl;
    }
  vector<vector<Quad>> rsfours; // {r,s} pairs in a foursome

  for (const auto& rlists : alist_by_denom)
    {
      Quad s = rlists.first;
      Quad sbar = s.conj();
      auto rlist = rlists.second;
      if (debug)
        cout << "s = " << s << ": rlist = " << rlist << endl;
      for (const auto& r : rlist)
        {
          assert (r==rectify(r,s));
          Quad mr = rectify(-r,s);
          vector<Quad> rs = {r,s}, mrs = {mr,s};
          if (debug) cout << "(r,s) = " << rs <<"; (-r,s) = "<< mrs << endl;
          if (std::find(pluspairs.begin(), pluspairs.end(), rs) != pluspairs.end())
            continue;
          if (std::find(pluspairs.begin(), pluspairs.end(), mrs) != pluspairs.end())
            continue;
          if (std::find(minuspairs.begin(), minuspairs.end(), rs) != minuspairs.end())
            continue;
          if (std::find(minuspairs.begin(), minuspairs.end(), mrs) != minuspairs.end())
            continue;
          if (std::find(rsfours.begin(), rsfours.end(), rs) != rsfours.end())
            continue;
          if (std::find(rsfours.begin(), rsfours.end(), mrs) != rsfours.end())
            continue;
          Quad r2=r*r;
          if (div(s, r2-1)) // r^2 = +1 (mod s) : "plus pair"
            {
              pluspairs.push_back(rs);
              if (verbose) cout << "+ pair "<<rs<<endl;
              continue;
            }
          if (div(s, r2+1)) // r^2 = -1 (mod s) : "minus pair"
            {
              minuspairs.push_back(rs);
              if (verbose) cout << "- pair "<<rs<<endl;
              continue;
            }
          // look for a four, i.e. r*r' = -1 (mod s) with r!=+-r'
          auto p = find_if(rlist.begin(), rlist.end(), [s,r](const Quad& rd) {return div(s, r*rd+1);});
          assert (p!=rlist.end());
          Quad rd = *p;
          assert (rd==rectify(rd,s));
          assert (div(s,r*rd+1));
          Quad mrd = rectify(-rd,s);
          vector<Quad> rds = {rd,s}, mrds = {mrd,s};
          if (debug) cout << "(r',s) = " << rds <<"; (-r',s) = "<< mrds << endl;
          if (std::find(rsfours.begin(), rsfours.end(), rds) == rsfours.end()
              &&
              std::find(rsfours.begin(), rsfours.end(), mrds) == rsfours.end()
              )
            {
              rsfours.push_back(rs);
              rsfours.push_back(rds);
              vector<Quad> sr1r2 = {s,r,rd};
              fours.push_back(sr1r2);
              if (verbose)
                cout << " foursome " << s << " " << r << " " << rd << endl;
            }
          else
            {
              if (verbose)
                cout<<" - skipping, seen before"<<endl;
            }
        }
    }

  // for homology edge relation computation include alphas with denom 1,2,3:
  Quad w = Quad::w;
  int d = Quad::d;
  // denom 1:
  minuspairs.push_back({ZERO,ONE});
  // denom 2:
  switch (d%4) {
  case 1:
    minuspairs.push_back({w,TWO});  // w^2=-1 (mod 2)
                                    // could also go into pluspairs
    break;
  case 2:
    minuspairs.push_back({w+1,TWO}); // (w+1)^2=-1 (mod 2)
                                     // could also go into pluspairs
    break;
  case 3:
    fours.push_back({TWO,w,w+1}); // w(w+1)=-1 (mod 2)
      };
  // denom 3:
  switch (d%12) {
  case 1: case 10:
    minuspairs.push_back({w, THREE});          // w^2=-1 (mod 3)
    fours.push_back({THREE, 1+w, 1-w});   // (1+w)(1-w)=-1 (mod 3)
    break;
  case 7:
    if (Quad::d>19)
      minuspairs.push_back({1+w, THREE});        // (1+w)^2=-1 (mod 3)
    if (Quad::d>31)
      fours.push_back({THREE, w, 1-w}); // w(1-w)=-1 (mod 3)
    break;
  case 2: case 5:
    if (d>5)
      pluspairs.push_back({w, THREE});           // w^2=+1 (mod 3)
    break;
  case 11:
    if (d>23)
      pluspairs.push_back({1+w, THREE});         // (1+w)^2=+1 (mod 3)
    break;
  case 3:
    if (d>15)
      fours.push_back({THREE, w, w-1});     // w(w-1)=-1 (mod 3)
    break;
  case 6: case 9:
    if (d>6)
      fours.push_back({THREE, w+1, w-1});   // (w+1)(w-1)=-1 (mod 3)
    break;
  }
  alist.insert(alist.end(), a1list.begin(), a1list.end());
  alist.insert(alist.end(), a2list.begin(), a2list.end());
  alist.insert(alist.end(), a3list.begin(), a3list.end());
  for (const auto& rs : pluspairs)
    {
      Quad r = rs[0], s=rs[1];
      RatQuad a(r,s);
      alist.push_back(a);
      alist.push_back(-a);
    }
  for (const auto& rs : minuspairs)
    {
      Quad r = rs[0], s=rs[1];
      RatQuad a(r,s);
      alist.push_back(a);
      alist.push_back(-a);
    }
  for (const auto& sr1r2 : fours)
    {
      Quad s = sr1r2[0], r1=sr1r2[1], r2=sr1r2[2];
      RatQuad a1(r1,s), a2(r2,s);
      alist.push_back(a1);
      alist.push_back(-a1);
      alist.push_back(a2);
      alist.push_back(-a2);
    }
  return alist;
}

// Output sorted list of alphas (denom > 3 in pairs or fours)

void output_alphas(vector<vector<Quad>>& pluspairs, vector<vector<Quad>>& minuspairs, vector<vector<Quad>>& fours,
                   int to_file, int to_screen)
{
  vector<Quad> small_denoms = {Quad::zero, Quad::one, 2*Quad::one, 3*Quad::one};
  ofstream geodata;
  stringstream ss;
  if (to_file)
    {
      ss << "geodata_" << Quad::d << ".dat";
      geodata.open(ss.str().c_str(), ios_base::app);
    }
  int nlines=0;
  for ( const auto& rs : pluspairs)
    {
      Quad s = rs[1];
      if (std::find(small_denoms.begin(), small_denoms.end(), s) != small_denoms.end())
        continue;
      Quad r = rectify(rs[0], s);
      if (!pos(r)) r = rectify(-r,s);
      nlines++;
      if (to_file) // output s, r, -r
        {
          geodata << Quad::d << " A ";
          geodata << s.re() << " " << s.im() << " ";
          geodata << r.re() << " " << r.im() << " ";
          geodata << -r.re() << " " << -r.im() << endl;
        }
      if (to_screen)
        {
          cout << Quad::d << " A ";
          cout << s.re() << " " << s.im() << " ";
          cout << r.re() << " " << r.im() << " ";
          cout << -r.re() << " " << -r.im() << endl;
        }
    }
  for ( const auto& rs : minuspairs)
    {
      Quad s = rs[1];
      if (std::find(small_denoms.begin(), small_denoms.end(), s) != small_denoms.end())
        continue;
      Quad r = rectify(rs[0], s);
      if (!pos(r)) r = rectify(-r,s);
      nlines++;
      if (to_file) // output s, r, r
        {
          geodata << Quad::d << " A ";
          geodata << s.re() << " " << s.im() << " ";
          geodata << r.re() << " " << r.im() << " ";
          geodata << r.re() << " " << r.im() << endl;
        }
      if (to_screen)
        {
          cout << Quad::d << " A ";
          cout << s.re() << " " << s.im() << " ";
          cout << r.re() << " " << r.im() << " ";
          cout << r.re() << " " << r.im() << endl;
        }
    }
  for ( const auto& sr1r2 : fours)
    {
      Quad s = sr1r2[0];
      if (std::find(small_denoms.begin(), small_denoms.end(), s) != small_denoms.end())
        continue;
      Quad r1 = rectify(sr1r2[1],s);
      Quad r2 = rectify(sr1r2[2],s);
      if (!pos(r1)) {r1 = rectify(-r1,s); r2 = rectify(-r2,s);}
      nlines++;
      if (to_file) // output s, r1, r2
        {
          geodata << Quad::d << " A ";
          geodata << s.re() << " " << s.im() << " ";
          geodata << r1.re() << " " << r1.im() << " ";
          geodata << r2.re() << " " << r2.im() << endl;
        }
      if (to_screen)
        {
          cout << Quad::d << " A ";
          cout << s.re() << " " << s.im() << " ";
          cout << r1.re() << " " << r1.im() << " ";
          cout << r2.re() << " " << r2.im() << endl;
        }
    }
}
