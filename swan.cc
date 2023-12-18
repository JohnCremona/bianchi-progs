// FILE SWAN.CC: implementation of Swan's algorithm

#include <iostream>

#include "swan.h"
#include "geometry.h"
#include "looper.h"

// Given an ideal I, return a list of singular points of class [I]
// (one representative for each orbit under integral translations).

CuspList singular_points_in_class(Qideal I, int verbose)
{
  if (I.is_principal())
    return {RatQuad(Quad::one, Quad::zero)};
  INT n = I.norm();
  Quad r, s(n);
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
      vector<Quad> rlist = Qideal(s).residues();
      if (verbose)
        cout<<"Residues modulo "<<s<<" are "<<rlist<<endl;
      for ( auto r : rlist)
        {
          r = r%s;
          if (I==Qideal({r,s}))
            {
              RatQuad sig = reduce_to_rectangle(RatQuad(r,s));
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

CuspList sort_singular_points(const CuspList S, int verbose)
{
  CuspList sorted_S;
  sorted_S.push_back(RatQuad(Quad::one, Quad::zero));
  if (Quad::class_number==1)
    return sorted_S;

  // denom 2 sigmas we construct directly
  Quad w = Quad::w;
  long d = Quad::d;

  switch (d%8) {
  case 1: case 5:
    sorted_S.push_back(RatQuad(1+w,TWO));
    break;
  case 2: case 6:
    sorted_S.push_back(RatQuad(w,TWO));
    break;
  case 7:
    sorted_S.push_back(RatQuad(w,TWO));
    sorted_S.push_back(RatQuad(1-w,TWO));
    break;
  default:
    ;
  }

  // denom 3 sigmas we construct directly
  switch (d%12) {
  case 2: case 5:
    if (d>5)
      {
        sorted_S.push_back(RatQuad(1+w,THREE));
        sorted_S.push_back(RatQuad(-1-w,THREE));
        sorted_S.push_back(RatQuad(1-w,THREE));
        sorted_S.push_back(RatQuad(-1+w,THREE));
      }
    break;
  case 3:
    if (d>15)
      {
        sorted_S.push_back(RatQuad(1+w,THREE));
        sorted_S.push_back(RatQuad(-1-w,THREE));
      }
    break;
  case 6: case 9:
    if (d>6)
      {
        sorted_S.push_back(RatQuad(w,THREE));
        sorted_S.push_back(RatQuad(-w,THREE));
      }
    break;
  case 11:
    if (d>23)
      {
        sorted_S.push_back(RatQuad(w,THREE));
        sorted_S.push_back(RatQuad(-w,THREE));
        sorted_S.push_back(RatQuad(1-w,THREE));
        sorted_S.push_back(RatQuad(-1+w,THREE));
      }
    break;
  default:
    ;
  }

  if (verbose)
    cout<<"Sigmas with small denominators: "<<sorted_S<<endl;

  // Now process the other sigmas (if any)
  RAT half(1,2);
  for ( auto s : S) // not const or reference as we may change it
    {
      if (verbose)
        cout <<"sigma = "<<s<<endl;

      // skip oo and denom 2 or 3 sigmas:
      if (s.is_infinity() or (TWO*s).is_integral() or (THREE*s).is_integral())
        {
          if (verbose)
            cout <<" - skipping (small denominator)"<<endl;
          continue;
        }

      // skip sigmas if we have seen its negative:
      Quad t;
      if (cusp_index_with_translation(-s, sorted_S, t) != -1)
        {
          if (verbose)
            cout <<" - skipping (negative seen already)";
          continue;
        }

      // Standardise if real or imaginary parts are +-1/2:
      RAT r(s.real()), i(s.imag());
      if (Quad::t==0)
        {
          if (r<0 && i==half)
            s -= w;
          else
            if (i<0 && r==half)
              s -= ONE;
        }
      else
        {
          if (i==half)
            {
              if (r>0)
                s -= w;
              else
                if (r<-half)
                  s += 1-w;
            }
          else
            if (2*r+i==1 && i<0)
              s -= ONE;
        }
      r = s.real();
      i = s.imag();
      if (verbose)
        cout <<"sigma = "<<s<<", -sigma = "<<-s<<endl;
      if (i>0)
        {
          sorted_S.push_back(s);
          sorted_S.push_back(-s);
        }
      else
        {
          sorted_S.push_back(-s);
          sorted_S.push_back(s);
        }
    }
  if (verbose)
    cout<<"All sigmas: "<<sorted_S<<endl;
  return sorted_S;
}


// Output sorted list of singular points (oo, denom 2, denom 3, larger denoms in +/- pairs)

void output_singular_points(const CuspList S, int to_file, int to_screen)
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
      if ((s.imag()>0) && (std::find(small_denoms.begin(), small_denoms.end(), sden) == small_denoms.end()))
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
    return 0;
  int t  = tau(a1,a2);
  return (strict? t==-2: t<0);
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
      Quad root_minus_d = (Quad::t ? 2*Quad::w-1 : Quad::w);
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
    ans.push_back(z/TWO + sqrtd2/(TWO*delta));
  return ans;
}


// return 1 iff a is [strictly] inside S_b
int is_inside(const RatQuad& a, const RatQuad& b, int strict)
{
  RAT t = radius_squared(b) - (a-b).norm();
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
  return principal_cusps_with_denominators(elements_of_norm(n));
}

// list of principal cusps with denominator norm up to given bound,
// omitting any whose circles are contained in an earlier one.

CuspList principal_cusps_of_dnorm_up_to(const INT& maxn)
{
  return principal_cusps_with_denominators(elements_of_norm_up_to(maxn));
}

// list of principal cusps with given denominator
CuspList principal_cusps_with_denominator(const Quad& s)
{
  CuspList alist;
  auto rlist = invertible_residues(s);
  for ( const auto& r : rlist)
    {
      RatQuad a = reduce_to_rectangle(RatQuad(r, s));
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

// Given principal cusps a1, a2, a such that the circles S_a1 and
// S_a2 intersect in distinct points, test whether S_a covers either
// or both these points.

// Returns 0 if neither, 2 if both, +1 or -1 if just one.  The signs
// are consistent so that if a returns +1 and a' returns -1 then each
// intersection point is covered by either S_a or S_a'.

//#define DEBUG_ARE_INTERSECTION_POINTS_COVERED_BY_ONE

int are_intersection_points_covered_by_one(const RatQuad& a1, const RatQuad& a2, const RatQuad& a)
{
#ifdef DEBUG_ARE_INTERSECTION_POINTS_COVERED_BY_ONE
  cout << "Testing if intersection points of "<<a1<<" and "<<a2<<" are covered by "<<a<<endl;
#endif

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
  RAT D2 = (D*D).real();          // negative rational
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
      RAT u = (D*(w-w.conj())).real();
      code = u.sign();
    }
#ifdef DEBUG_ARE_INTERSECTION_POINTS_COVERED_BY_ONE
  cout << " - test returns code "<<code<<endl;
#endif
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
  CuspList zlist = intersection_points_in_k(a0,a1);
  if (!zlist.empty())
    {
      for (auto z = zlist.begin(); z!=zlist.end(); ++z)
        {
          if (cusp_index(*z, sigmas)!=-1) // z is singular: OK
            continue;
          if (!is_inside_one(*z, alist))  // z is not covered: not OK
            return 0;
        }
      // we reach here if every z is either singular or covered: OK
      return 1;
    }

  // Now the intersection points are not in k. Check that either one
  // S_a covers both, or two cover one each:

  int t = 0; // will hold +1 or -1 if we have covered only one of the two
  for ( const auto& a2 : alist)
    {
      if (a2 == a1)
        continue;
      int t2 = are_intersection_points_covered_by_one(a0, a1, a2);
      if (t2==2) // both are covered by a2
        return 1;
      if (t2==0) // neither is covered by a2
        continue;
      // Now t2 is +1 or -1; we win if t is its negative
      if (t==-t2) // then they are (1,-1) or (-1,1) so we have covered both points
        return 1;
      t = t2;     // = +1 or -1: we have covered one of the points, so remember which
    }
  // If we reach here then none of the S_a covers both points
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
                        vector<CuspPair>& pairs_ok)
{
  // check if S_a0 is strictly entirely contained in one S_alpha:
  if (circle_inside_any_circle(a0, alist, 1))
    return 1;

  // extract the relevant alphas, if any, namely those for which
  // S_alpha and S_a0 properly intersect:

  CuspList a0list(alist.size()); // upper bound on size
  auto it1 = std::copy_if(alist.begin(), alist.end(), a0list.begin(),
                         [a0](RatQuad a){return circles_intersect(a0, a);});
  a0list.resize(std::distance(a0list.begin(),it1));  // shrink to new size

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
        pairs_ok.push_back(pr); // record that this pair is ok
      else
        all_ok = 0; // but continue checking the other a1s
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
                          const CuspList& slist, vector<CuspPair>& pairs_ok)
{
  // cout << "At start of are_alphas_surrounded()" << endl;
  // cout << "alist_ok = "<<alist_ok<<endl;
  // cout << "alist_open = "<<alist_open<<endl;

  // concatenate the two alists
  CuspList alist;
  alist.insert(alist.end(), alist_ok.begin(), alist_ok.end());
  alist.insert(alist.end(), alist_open.begin(), alist_open.end());

  // add translates
  vector<Quad> translates = { Quad::one, Quad::w, Quad::one+Quad::w, Quad::one-Quad::w };
  CuspList alistx = alist;
  for ( auto& a : alist)
    for (auto& t : translates)
      {
        alistx.push_back(a+t);
        alistx.push_back(a-t);
      }
  int i=0, n_open = alist_open.size();

  // We make a copy of alist_open to loop over, so we can delete elements from the original as we go
  CuspList new_alist_open = alist_open;
  for ( const auto& a : new_alist_open)
    {
      assert (a.in_quarter_rectangle());
      i++;
      cout <<"Testing alpha #"<<i<<"/"<<n_open<<" = "<<a<<"...";
      if (is_alpha_surrounded(a, alistx, slist, pairs_ok))
        {
          cout << " ok! surrounded\n";
          // add this alpha to the ok list end remove from the open list
          alist_ok.push_back(a);
          alist_open.erase(std::find(alist_open.begin(), alist_open.end(), a));
        }
      else
        {
          cout << " no, not surrounded" << endl;
          return 0;
        }
    }
  return 1;
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
      maxn = new_alphas.back().den().norm();

      if (verbose)
        cout << "----------------------------------------------------------\n"
             << "Considering "<<new_alphas.size()<<" extra principal cusps " << s << maxn
             // << ": "<<new_alphas
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
          if (circle_inside_any_circle(a, alist, 0)) // strict=0
              continue;
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

      if (are_alphas_surrounded(alphas_ok, alphas_open, sigmas, pairs_ok))
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
