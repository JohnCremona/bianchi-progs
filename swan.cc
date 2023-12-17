// FILE SWAN.CC: implementation of Swan's algorithm

#include <iostream>

#include "swan.h"
#include "geometry.h"
#include "looper.h"

// Given an ideal I, return a list of singular points of class [I]
// (one representative for each orbit under integral translations).

vector<RatQuad> singular_points_in_class(Qideal I, int verbose)
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
  vector<RatQuad> S;
  for (auto si=slist.begin(); si!=slist.end(); si++)
    {
      s = *si;
      if (verbose)
        cout<<" - using s = "<<s<<endl;
      vector<Quad> rlist = Qideal(s).residues();
      if (verbose)
        cout<<"Residues modulo "<<s<<" are "<<rlist<<endl;
      for (auto ri=rlist.begin(); ri!=rlist.end(); ri++)
        {
          r = (*ri)%s;
          Qideal J = Qideal({r,s});
          //cout<<" @ r = "<<r<<" gives (r,s) = "<<J<<endl;
          if (I==J)
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

vector<vector<RatQuad>> singular_points_by_class()
{
  vector<vector<RatQuad>> sigma_lists;
  for (auto I = Quad::class_group.begin(); I != Quad::class_group.end(); ++I)
    sigma_lists.push_back(singular_points_in_class(*I));
  return sigma_lists;
}

// Return one list of all singular points.

vector<RatQuad> singular_points()
{
  vector<RatQuad> sigma_list;
  for (auto I = Quad::class_group.begin(); I != Quad::class_group.end(); ++I)
    {
      vector<RatQuad> S = singular_points_in_class(*I);
      sigma_list.insert(sigma_list.end(), S.begin(), S.end());
    }
  return sigma_list;
}

// Return sorted list of singular points (oo, denom 2, denom 3, larger denoms in +/- pairs)

vector<RatQuad> sort_singular_points(const vector<RatQuad> S, int verbose)
{
  vector<RatQuad> sorted_S;
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
  for (auto si=S.begin(); si!=S.end(); ++si)
    {
      RatQuad s = *si;
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
      RAT r(s.real()), i(s.imag()), half(1,2);
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

void output_singular_points(const vector<RatQuad> S, int to_file, int to_screen)
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
  for (auto s = S.begin(); s!= S.end(); ++s)
    {
      Quad sden = s->den(), snum = s->num();
      if ((s->imag()>0) && (std::find(small_denoms.begin(), small_denoms.end(), sden) == small_denoms.end()))
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
vector<RatQuad> intersection_points_in_k(const RatQuad& a1, const RatQuad& a2)
{
  RAT r1sq = radius_squared(a1), r2sq = radius_squared(a2);
  RatQuad delta = a2-a1;
  RAT n = delta.norm();
  RAT d1 = n - (r1sq + r2sq);
  RAT d2 = d1*d1 - 4*r1sq*r2sq;

  vector<RatQuad> ans;
  if (d2 > 0)
    return ans;
  delta = delta.conj();
  RatQuad z = a1 + a2 + (r1sq-r2sq)/delta;
  // find sqrts of d2 in k, if any
  vector<RatQuad> d2sqrts = sqrts(d2);
  for (auto sqrtd2 = d2sqrts.begin(); sqrtd2!=d2sqrts.end(); ++sqrtd2)
    ans.push_back(z/TWO + (*sqrtd2)/(TWO*delta));
  return ans;
}


// return 1 iff a is [strictly] inside S_b
int is_inside(const RatQuad& a, const RatQuad& b, int strict)
{
  RAT t = radius_squared(b) - (a-b).norm();
  return (strict ? 0<t : 0<=t);
}

// return 1 iff a is [strictly] inside S_b for at least one b in blist
int is_inside_one(const RatQuad& a, const vector<RatQuad>& blist, int strict)
{
  for (auto b = blist.begin(); b!=blist.end(); ++b)
    if (is_inside(a, *b, strict))
      return 1;
  return 0;
}


// list of principal cusps with given denominator norm
vector<RatQuad> principal_cusps_of_norm(const INT& n)
{
  vector<RatQuad> alist;
  auto slist = elements_of_norm(n);
  // cout << "Elements of norm "<<n<<": "<<slist<<endl;
  for ( const auto& s : slist)
    {
      auto rlist = invertible_residues(s);
      // cout << "s="<<s<<": invertible residues: "<<rlist<<endl;
      for ( const auto& r : rlist)
        alist.push_back(reduce_to_rectangle(RatQuad(r, s)));
    }
  return alist;
}

// list of principal cusps with denominator norm up to given bound,
// omitting any whose circles are contained in an earlier one.
vector<RatQuad> principal_cusps_up_to(const INT& maxn)
{
  vector<RatQuad> alist;
  auto slist = elements_of_norm_up_to(maxn);
  // cout << "Elements of norm up to "<<maxn<<": "<<slist<<endl;
  for ( const auto& s : slist)
    {
      auto rlist = invertible_residues(s);
      // cout << "s="<<s<<": invertible residues: "<<rlist<<endl;
      for ( const auto& r : rlist)
        {
          RatQuad a = reduce_to_rectangle(RatQuad(r, s));
          if (!is_inside_one(a, alist))
            alist.push_back(a);
        }
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
    n = delta.norm();
  RatQuad z0 = ((a1+a2) + (r1sq-r2sq)/delta)/TWO;
  RAT
    d1 = n - (r1sq + r2sq),
    d2 = d1*d1 - 4*r1sq*r2sq;

  assert (d2 < 0);

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

int are_intersection_points_covered(const RatQuad& a0, const RatQuad& a1, const vector<RatQuad>& alist,
                                    const vector<RatQuad>& sigmas, int debug)
{
  vector<RatQuad> zlist = intersection_points_in_k(a0,a1);
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
  for (auto a2=alist.begin(); a2!=alist.end(); ++a2)
    {
      if (*a2 == a1)
        continue;
      int t2 = are_intersection_points_covered_by_one(a0, a1, *a2);
      if (t2==2)
        return 1;
      if (t2!=0) // then it is +1 or -1
        {
          if (t==-t2) // then they are (1,-1) or (-1,1) so we have covered both points
            return 1;
          t = t2;     // = +1 or -1: we have covered one of the points, so remember which
        }
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

// pairs_ok is a list of pairs (a1,a2) whose intersection points are
// known to be covered, so can be skipped.

// Returns 0 or 1, and pairs_ok is updated to a list of pairs
// whose intersections have now been shown to be covered.

int is_alpha_surrounded(const RatQuad& a0, const vector<RatQuad>& alist, const vector<RatQuad>& sigmas,
                        vector< pair<RatQuad,RatQuad> >& pairs_ok)
{
  // check if S_a0 is strictly entirely contained in one S_alpha:
  for (auto a = alist.begin(); a!=alist.end(); ++a)
    if (circle_inside_circle(a0, *a, 1))
      return 1;

  // extract the relevant alphas, if any, namely those for which
  // S_alpha and S_a0 properly intersect:

  vector<RatQuad> a0list;
  for (auto a = alist.begin(); a!=alist.end(); ++a)
    if (circles_intersect(a0, *a))
      a0list.push_back(*a);

  // extract the pairs in pairs_ok which contain a0:
  vector< pair<RatQuad,RatQuad> > a0_pairs_ok;
  for (auto pr=pairs_ok.begin(); pr!=pairs_ok.end(); ++pr)
    {
      if (a0==pr->second)
        a0_pairs_ok.push_back({a0, pr->first});
      if (a0==pr->first)
        a0_pairs_ok.push_back(*pr);
    }

  int all_ok = 1;
  for (auto a1=a0list.begin(); a1!=a0list.end(); ++a1)
    {
      pair<RatQuad,RatQuad> pr = {a0,*a1};
      if (std::find(a0_pairs_ok.begin(), a0_pairs_ok.end(), pr) != a0_pairs_ok.end())
        continue;
      if (are_intersection_points_covered(a0, *a1, alist, sigmas))
        pairs_ok.push_back(pr);
      else
        all_ok = 0;
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

// Returns 0/1

// NB All a in alist_open will be tested, i.e. we carry on after a
// failure.

int are_alphas_surrounded(vector<RatQuad>& alist_ok, vector<RatQuad>& alist_open,
                          const vector<RatQuad>& slist, vector< pair<RatQuad,RatQuad> >& pairs_ok)
{
  cout << "alist_ok = "<<alist_ok<<endl;
  cout << "alist_open = "<<alist_open<<endl;
  // concatenate the two alists
  vector<RatQuad> alist;
  alist.insert(alist.end(), alist_ok.begin(), alist_ok.end());
  alist.insert(alist.end(), alist_open.begin(), alist_open.end());

  // add translates
  RatQuad one(Quad::one), w(Quad::w);
  vector<RatQuad> translates = { one, w, one+w, one-w };
  cout << "translates = "<<translates<<endl;
  vector<RatQuad> alistx = alist;
  for ( auto& a : alist)
    for (auto& t : translates)
      {
        alistx.push_back(a+t);
        alistx.push_back(a-t);
      }
  cout << "alistx = "<<alistx<<endl;
  int ok, i=0, all_ok = 1;
  // We make a copy of this to loop over, so we can delete elements from the original as we go
  vector<RatQuad>& new_alist_open = alist_open;
  for ( const auto& a : new_alist_open)
    {
      i++;
      if (!all_ok) // then return False so don't check remaining alphas
        continue;
      if (a.in_quarter_rectangle())
        {
          cout <<"Testing alpha #"<<i<<"/"<<alist_open.size()<<" = "<<a<<"...";
          ok = is_alpha_surrounded(a, alistx, slist, pairs_ok);
          if (ok)
            cout << " ok! surrounded" << endl;
          else
            {
              all_ok = 0;
              cout << " no, not surrounded" << endl;
            }
        }
      else // alphas not in quarter rectangle are ignored by putting straight into the ok list
        {
          ok = 1;
        }
      if (ok) // add this alpha to the ok list end remove from the open list
        {
          alist_ok.push_back(a);
          alist_open.erase(std::find(alist_open.begin(), alist_open.end(), a));
        }
    }
  return all_ok;
}

// Returns a finite list of principal cusps a such that the S_{a+t}
// for all integral t cover CC apart from singular points.

// For n>=1 successively, we test as a candidate set all a=r/s with
// r,s coprime, r reduced mod s, N(s)<=n (omitting any for which S_a
// is contained in any earlier S_a') until we succeed.

// sigmas can be set to a list of singular points (up to translation),
// otherwise these will be computed.

// Other functions will then (1) saturate the set, (2) discard redundancies.

vector<RatQuad> covering_alphas(const vector<RatQuad>& sigmas, int verbose)
{
  INT maxn = 0;
  vector<RatQuad> alphas_ok, alphas_open, new_cusps, alist;
  vector< pair<RatQuad,RatQuad> > pairs_ok;
  int first = 1;
  while (1)
    {
      int nc = 0; // number of new alphas added to list
      if (first)
        {
          maxn = Quad::absdisc/4;
          new_cusps = principal_cusps_up_to(maxn);
          first = 0;
        }
      else
        {
          maxn += 1;
          new_cusps = principal_cusps_of_norm(maxn);
          while (new_cusps.empty())
            {
              maxn += 1;
              new_cusps = principal_cusps_of_norm(maxn);
            }
        }
      if (verbose)
        cout << "Added "<<new_cusps.size()<<" extra principal cusps alpha, norms up to "<<maxn<<endl;
      for (auto a : new_cusps)
        {
          cout << "Working on a="<<a<<endl;
          // test if a is redundant
          int red=0;
          for (auto b : alist)
            {
              red = circle_inside_circle(a, b, 0);
              if (red) break;
            }
          if (red) // skip a
            {
              cout << " - redundant, skipping"<<endl;
              continue;
            }
          alist.push_back(a);
          if (a.in_quarter_rectangle())
            {
              cout << " - in quarter rectangle, using"<<endl;
              alphas_open.push_back(a);
              nc += 1;
            }
          else
            {
              cout << " - not in quarter rectangle, not checking"<<endl;
              alphas_ok.push_back(a);
            }
        }
      if (verbose && nc)
        {
          cout << "Adding "<<nc<<" alphas of norm ";
          cout << (first? "up to " : "") << maxn;
          cout <<" (plus symmetrics); #alphas="<<alphas_ok.size()+alphas_open.size();
          cout <<" of which "<<alphas_ok.size()<<" are proved surrounded so far";
        }
      if (nc==0)
        continue;
      if (are_alphas_surrounded(alphas_ok, alphas_open, sigmas, pairs_ok))
        {
          if (verbose)
            cout << "Success using "<<alphas_ok.size()<<" alphas of with max norm "<<maxn<<"!";
          return alphas_ok;
        }
      else
        {
          if (verbose)
            cout << "Some alphas are not surrounded, continuing...";
        }
    }
}
