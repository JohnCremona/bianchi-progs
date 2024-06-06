// FILE SWAN_ALPHAS.CC: implementation of alpha-finding functions for Swan's algorithm

#include "swan_alphas.h"
#include "looper.h"

// direct lists of alphas and sigmas of denominator 2 or 3:
CuspList denom_2_alphas()
{
  if (Quad::is_Euclidean) return {};
  int d8 = (Quad::d)%8;
  if (d8==1 || d8==5) return {{0,1,2}}; // w/2
  if (d8==2 || d8==6) return {{1,1,2}};  // (1+w)/2
  if (d8==3) return {{0,1,2}, {-1,1,2}}; // w/2, (w-1)/2
  // now d8=7
  return {};
}

CuspList denom_3_alphas()
{
  if (Quad::is_Euclidean) return {};
  int d = Quad::d;
  if (d==5 || d==6 || d==15 || d==19 || d==23) return {};

  int d12 = (Quad::d)%12;
  // fill a3list up to sign
  switch (d12) {
  case 3:
    return {{0,1,3}, {0,-1,3}, {-1,1,3}, {1,-1,3}}; // {w, w-1}/3 and negs
  case 7:
    if (d>31)
      return {{1,1,3}, {-1,-1,3}, {0,1,3}, {0,-1,3}, {1,-1,3}, {-1,1,3}}; // {1+w, w, 1-w}/3 and negs
    else
      return {{1,1,3}, {-1,-1,3}}; // {1+w}/3 and negs
  case 11:
    return {{1,1,3}, {-1,-1,3}}; // {1+w}/3 and negs
  case 1: case 10:
    return {{0,1,3}, {0,-1,3}, {1,1,3}, {-1,-1,3}, {1,-1,3}, {-1,1,3}}; // {w, 1+w, 1-w}/3 and negs
  case 2: case 5:
    return {{0,1,3}, {0,-1,3}}; // {w}/3 and negs
  case 6: case 9:
    return {{1,1,3}, {-1,-1,3}, {-1,1,3}, {1,-1,3}}; // {w+1, w-1}/3 and negs
  default:
    return {};
  }
}

// Return sorted list of (saturated covering) alphas (denom 1, denom 2, denom 3, larger denoms in pairs or fours)
// + pluspairs, minuspairs, fours
CuspList sort_alphas(const CuspList& A,
                     vector<vector<Quad>>& pluspairs, vector<vector<Quad>>& minuspairs, vector<vector<Quad>>& fours,
                     int verbose, int debug)
{
  CuspList alist, a1list, a2list, a3list;
  map<Quad, vector<Quad>> alist_by_denom; // list of numerators for each denominator
  Quad temp;
  for (const auto& a0 : A)
    {
      RatQuad a = reduce_to_rectangle(a0, temp);
      if (!temp.is_zero() && debug)
        cout<<"sort_alphas replacing "<<a0<<" by "<<a<<", its translate by "<<temp<<endl;
      //assert (a.in_rectangle());
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
      alist_by_denom[s].push_back(r);
    }
  if (debug)
    {
      cout << " alphas with denominator 1: " <<a1list <<endl;
      cout << " alphas with denominator 2: " <<a2list <<endl;
      cout << " alphas with denominator 3: " <<a3list <<endl;
    }
  assert (a1list.size()==1 && a1list[0]==RatQuad(Quad(0)));
  auto d2s_expected = denom_2_alphas(), d3s_expected = denom_3_alphas();
  assert (compare_CuspLists_as_sets_mod_translation(a2list, d2s_expected));
  // cout << "a3list (actual):   " <<a3list << endl;
  // cout << "a3list (expected): " <<d3s_expected << endl;
  assert (compare_CuspLists_as_sets_mod_translation(a3list, d3s_expected));
  // We rely on the alphas of denom 1, 2, 3 begin in an exact order:
  if (a2list != d2s_expected)
    {
      if (verbose||debug) cout << "replacing denom 2 alphas by standard list "<<d2s_expected << endl;
      a2list = d2s_expected;
    }
  if (a3list != d3s_expected)
    {
      if (verbose||debug) cout << "replacing denom 3 alphas by standard list "<<d3s_expected << endl;
      a3list = d3s_expected;
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

  // Make a new sorted list of all alphas:

  alist.insert(alist.end(), a1list.begin(), a1list.end());
  alist.insert(alist.end(), a2list.begin(), a2list.end());
  alist.insert(alist.end(), a3list.begin(), a3list.end());
  if (debug)
    cout << "After inserting alphas of denom 1,2,3, alist = "<<alist<<endl;

  // NB the code in geometry.cc relies on these being exact a, -a pairs
  for (const auto& rs : pluspairs)
    {
      Quad r = rs[0], s=rs[1];
      alist.push_back(RatQuad(r,s));
      alist.push_back(RatQuad(-r,s));
    }
  if (debug)
    cout << "After inserting pluspairs, alist = "<<alist<<endl;
  for (const auto& rs : minuspairs)
    {
      Quad r = rs[0], s=rs[1];
      alist.push_back(RatQuad(r,s));
      alist.push_back(RatQuad(-r,s));
    }
  if (debug)
    cout << "After inserting minuspairs, alist = "<<alist<<endl;
  for (const auto& sr1r2 : fours)
    {
      Quad r = sr1r2[1], s=sr1r2[0];
      alist.push_back(RatQuad(r,s));
      alist.push_back(RatQuad(-r,s));
      r = sr1r2[2];
      alist.push_back(RatQuad(r,s));
      alist.push_back(RatQuad(-r,s));
    }
  if (debug)
    cout << "After inserting fours, alist = "<<alist<<endl;

  // include alphas with denom 1,2,3 in pluspairs, minuspairs and fours:
  Quad w = Quad::w, zero(0), one(1), two(2), three(3);
  int d = Quad::d;

  // denom 1:
  minuspairs.push_back({zero,one});

  // denom 2:
  if (!Quad::is_Euclidean)
    {
      switch (d%8) {
      case 1: case 5:
        minuspairs.push_back({w,two});  // w^2=-1 (mod 2)
        break;
      case 2: case 6:
        minuspairs.push_back({w+1,two}); // (w+1)^2=-1 (mod 2)
        break;
      case 3:
        fours.push_back({two,w,w-1}); // w(w-1)=-1 (mod 2)
      };
    }
  // denom 3:
  if (!Quad::is_Euclidean && !( (d==5 || d==6 || d==15 || d==19 || d==23) ))
    {
      switch (d%12) {
      case 1: case 10:
        minuspairs.push_back({w, three});          // w^2=-1 (mod 3)
        fours.push_back({three, 1+w, 1-w});   // (1+w)(1-w)=-1 (mod 3)
        break;
      case 7:
        minuspairs.push_back({1+w, three});        // (1+w)^2=-1 (mod 3)
        if (Quad::d>31)
          fours.push_back({three, w, 1-w}); // w(1-w)=-1 (mod 3)
        break;
      case 2: case 5:
        pluspairs.push_back({w, three});           // w^2=+1 (mod 3)
        break;
      case 11:
        pluspairs.push_back({1+w, three});         // (1+w)^2=+1 (mod 3)
        break;
      case 3:
        fours.push_back({three, w, w-1});     // w(w-1)=-1 (mod 3)
        break;
      case 6: case 9:
        fours.push_back({three, w+1, w-1});   // (w+1)(w-1)=-1 (mod 3)
        break;
      }
    }
  if (debug)
    cout << "sort_alphas() returns "<<alist<<endl;

  return alist;
}

// Output sorted list of alphas (denom > 3 in pairs or fours)

string make_A_line(const Quad& s, const Quad& r1, const Quad& r2)
{
  ostringstream ost;
  ost << Quad::d << " A ";
  ost << s.re() << " " << s.im() << " ";
  ost << r1.re() << " " << r1.im() << " ";
  ost << r2.re() << " " << r2.im();
  return ost.str();
}

CuspList alpha_orbits(const CuspList& alist, vector<vector<Quad>>& triples, int verbose, int debug)
{
  vector<vector<Quad>> pluspairs, minuspairs, fours;
  CuspList new_alist = sort_alphas(alist, pluspairs, minuspairs, fours, verbose, debug);
  triples.clear();
  triples.reserve(pluspairs.size() + minuspairs.size() + fours.size());

  // We don't just concatenate as we want denoms 1, 2, 3 first
  // denoms 1 and 2 are in minuspairs not pluspairs

  // Denominator 2 is special as r=-r (mod1); 0/1 was put into minuspairs
  const Quad one(1);
  for ( const auto& rs : minuspairs)
    if (rs[1]==one)
      {
        triples.push_back({rs[1], rs[0], rs[0]});
        break; // there's at most one
      }

  // Denominator 2 is special as r=-r (mod2); these were was put into minuspairs or fours
  const Quad two(2);
  for ( const auto& rs : minuspairs)
    if (rs[1]==two)
      {
        triples.push_back({rs[1], rs[0], rs[0]});
        break; // there's at most one
      }
  for ( const auto& sr1r2 : fours)
    if (sr1r2[0]==two)
      {
        triples.push_back(sr1r2);
        break; // there's at most one
      }

  // Denominator 3 does not need special treatment, but we do these
  // before general denominators just so that they come before them in
  // the list of triples
  const Quad three(3);
  for ( const auto& rs : pluspairs)
    if (rs[1]==three)
      triples.push_back({rs[1], rs[0], -rs[0]});
  for ( const auto& rs : minuspairs)
    if (rs[1]==three)
      triples.push_back({rs[1], rs[0], rs[0]});
  for ( const auto& sr1r2 : fours)
    if (sr1r2[0]==three)
      triples.push_back(sr1r2);

  for ( const auto& rs : pluspairs)
    {
      Quad s = rs[1];
      if (s!=one && s!=two && s!=three)
        triples.push_back({s, rs[0], -rs[0]});
    }
  for ( const auto& rs : minuspairs)
    {
      Quad s = rs[1];
      if (s!=one && s!=two && s!=three)
        triples.push_back({s, rs[0], rs[0]});
    }
  for ( const auto& sr1r2 : fours)
    {
      Quad s = sr1r2[0];
      if (s!=one && s!=two && s!=three)
        triples.push_back(sr1r2);
    }
  if (debug)
    cout << "alpha_orbits() returns "<<new_alist<<endl;
  return new_alist;
}

void output_alphas(const vector<vector<Quad>>& triples,
                   int to_file, int to_screen)
{
  vector<Quad> small_denoms = {Quad(0), Quad(1), Quad(2), Quad(3)};
  ofstream geodata;
  stringstream ss;
  if (to_file)
    {
      ss << "geodata_" << Quad::d << ".dat";
      geodata.open(ss.str().c_str()); //  , ios_base::app);
    }
  for ( const auto& sr1r2 : triples)
    {
      if (std::find(small_denoms.begin(), small_denoms.end(), sr1r2[0]) != small_denoms.end())
        continue;
      string st = make_A_line(sr1r2);
      if (to_file)
        geodata << st <<endl;
      if (to_screen)
        cout << st << endl;
    }
}

void output_alphas(const vector<vector<Quad>>& pluspairs,
                   const vector<vector<Quad>>& minuspairs,
                   const vector<vector<Quad>>& fours,
                   int to_file, int to_screen)
{
  vector<Quad> small_denoms = {Quad(0), Quad(1), Quad(2), Quad(3)};
  ofstream geodata;
  stringstream ss;
  if (to_file)
    {
      ss << "geodata_" << Quad::d << ".dat";
      geodata.open(ss.str().c_str()); //  , ios_base::app);
    }
  int nlines=0;

  auto output3 = [small_denoms, to_file, to_screen, &geodata](const Quad& s, const Quad& r1, const Quad& r2)
    {
      if (std::find(small_denoms.begin(), small_denoms.end(), s) != small_denoms.end())
        return 0;
      string st = make_A_line(s, r1, r2);
      if (to_file)
        geodata << st <<endl;
      if (to_screen)
        cout << st << endl;
      return 1;
    };
  for ( const auto& rs : pluspairs)
    nlines += output3(rs[1], rs[0], -rs[0]);
  for ( const auto& rs : minuspairs)
    nlines += output3(rs[1], rs[0], rs[0]);
  for ( const auto& sr1r2 : fours)
    nlines += output3(sr1r2[0], sr1r2[1], sr1r2[2]);
  // if (to_file || to_screen)
  //   cout << nlines << " A-lines output"<<endl;
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

  RatQuad z0 = ((a1+a2) + RatQuad(r1sq-r2sq)/delta)/2;
  RAT T = 2 * n * (rsq - (z0-a).norm()) + d2/2;
  int Tsign = T.sign();
  RAT T2 = T*T;
  RatQuad D = tri_det(a, a2, a1); // pure imaginary
  RAT D2 = (D*D).x_coord(1);      // negative rational
  RAT d2D2 = d2*D2;               // positive rational

  // the covering condition is \pm sqrt(d2)*D < T

  int code;
  switch ((T2-d2D2).sign()) {
  case 0:
  default:
    code = 0;
    break;
  case 1:
    code = 1+Tsign; // =0 or 2, never 1
    break;
  case -1:
    Quad w = Quad::w;
    code = ((D*(w-w.conj())).x_coord(1)).sign();
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
  CuspList a0list = intersecting_alphas(a0,alist);
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
      if (are_intersection_points_covered(a0, a1, alist, sigmas, debug))
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
  CuspList alistx = cusp_shifts(alist, Quad::shifts_by_one);

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

CuspList covering_alphas(const CuspList& slist, int verbose)
{
  CuspList alphas_ok;  // list that will be returned
  INT maxn = Quad::absdisc; ///4; // we first consider all alphas with dnorm up to this
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
                             principal_cusps_with_denominators(looper.values_with_norm_up_to(maxn, 1))
                             :
                             principal_cusps_with_denominators(looper.values_with_current_norm())
                             );
      // cout << "new_alphas = " << new_alphas << endl;
      maxn = new_alphas.back().den().norm();

      if (verbose)
        cout << "----------------------------------------------------------\n"
             << "Considering "<<new_alphas.size()<<(first?"":" extra") << " principal cusps " << s << maxn
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
              // cout<<"alpha = "<<a<<" is weakly inside one of "<<alist<<endl;
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

      if (are_alphas_surrounded(alphas_ok, alphas_open, slist, pairs_ok, verbose, verbose>1))
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

H3pointList triple_intersections(const CuspList& alphas, int debug)
{
  H3pointList old_points = old_triple_intersections(alphas, debug);
  if (0)
    {
      H3pointList new_points = new_triple_intersections(alphas, debug);
      cout << "Old list of corners (size "<<old_points.size()<<") : "<<old_points<<endl;
      cout << "New list of corners (size "<<new_points.size()<<") : "<<new_points<<endl;
      if (old_points.size()!=new_points.size())
        cout << "Numbers differ!" <<endl;
    }
  // return new_points;
  return old_points;
}

H3pointList old_triple_intersections(const CuspList& alphas, int debug)
{
    if (debug)
      cout << "Finding triple intersections for " <<alphas.size()<<" alphas (old code)..."<<endl;

    // Extract the alphas in F4.  NB the standard list of alphas has
    // the property that, after those of denom 1, 2, 3, they come in
    // pairs a, -a; the second of the pair is not always in the
    // rectangle (because of the rounding of 1/2) which causes a
    // problem here:

    // CuspList alphasF4(alphas.size());
    // auto it1 = std::copy_if(alphas.begin(), alphas.end(), alphasF4.begin(),
    //                         [](RatQuad a){return a.in_quarter_rectangle();});
    // alphasF4.resize(std::distance(alphasF4.begin(),it1));  // shrink to new size
    CuspList alphasF4;
    for ( auto a : alphas)
      {
        a = reduce_to_rectangle(a);
        if (a.in_quarter_rectangle())
          alphasF4.push_back(a);
      }
    if (debug)
      cout << alphasF4.size() <<" alphas are in the quarter rectangle: " << alphasF4 << endl;

    // Extend these by 8 translations:
    CuspList alphasF4X;
    for ( auto& z : alphasF4)
      {
        CuspList z_nbrs = F4nbrs(z);
        alphasF4X.insert(alphasF4X.end(), z_nbrs.begin(), z_nbrs.end());
      }
    if (debug)
      cout << alphasF4X.size() <<" neighbours of these: " << alphasF4X << endl;

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
        i = i_j_list.first;
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
        i=ijk[0]; int j=ijk[1], k=ijk[2];
        H3pointList points1 = tri_inter_points(alphasF4X[i], alphasF4X[j], alphasF4X[k]);
        if (points1.empty())
          continue;
        H3point P = *points1.begin();
        RatQuad z = P.z;
        RAT t2 = P.t2;
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

H3pointList new_triple_intersections(const CuspList& alphas, int debug)
{
    if (debug)
      cout << "Finding triple intersections for " <<alphas.size()<<" alphas (new code)..."<<endl;

    // Extract the alphas in F4.  NB the standard list of alphas has
    // the property that, after those of denom 1, 2, 3, they come in
    // pairs a, -a; the second of the pair is not always in the
    // rectangle (because of the rounding of 1/2).
    CuspList alphasF4;
    for ( auto a : alphas)
      {
        a = reduce_to_rectangle(a);
        if (a.in_quarter_rectangle())
          alphasF4.push_back(a);
      }
    if (debug)
      cout << alphasF4.size() <<" alphas are in the quarter rectangle: " << alphasF4 << endl;

    // Extend these by 8 translations:
    CuspList alphasF4X;
    for ( auto& z : alphasF4)
      {
        CuspList z_nbrs = F4nbrs(z);
        alphasF4X.insert(alphasF4X.end(), z_nbrs.begin(), z_nbrs.end());
      }
    if (debug)
      cout << alphasF4X.size() <<" neighbours of these: " << alphasF4X << endl;

    H3pointList points;

    // Loop through all a0 in alphasF4:
    for (const auto& a0 : alphasF4)
      {
        // find its neighbours
        CuspList blist = intersecting_alphas(a0, alphasF4X);
        // For each intersecting pair of these we may have a triple
        // intersection point:
        int n = blist.size();
        if (n<3)
          return points; // empty list
        for (auto bi = blist.begin(); bi!=blist.end(); ++bi)
          for (auto bj = bi+1; bj!=blist.end(); ++bj)
            {
              RatQuad b1 = *bi, b2 = *bj;
              if (!circles_intersect(b1,b2))
                continue;
              H3pointList points1 = tri_inter_points(a0, b1, b2);
              if (points1.empty())
                continue;
              H3point P = points1.front();
              RatQuad z = P.z;
              RAT t2 = P.t2;
              if (t2.sign()==0)
                continue;
              if (!z.in_quarter_rectangle())
                continue;
              if (std::find(points.begin(), points.end(), P) != points.end())
                continue;
              if (is_under_any(P, alphasF4X))
                continue;
              if (debug)
                cout << " found P = "<<P<<", tri-intersection of "<<a0<<", "<<b1<<", "<<b2<<"\n";
              points.push_back(P);
              // These corners are in F4, so we apply symmetries to get all those in F:
              RatQuad zbar = z.conj();
              for (const auto& z2 : {-z, zbar, -zbar})
                {
                  if (!z2.in_rectangle())
                    continue;
                  H3point P2 = {z2, t2};
                  if (std::find(points.begin(), points.end(), P2) != points.end())
                    continue;
                  if (is_under_any(P, alphasF4X))
                    continue;
                  if (debug)
                    cout << " adding P2 = "<<P2<<" from ("<<a0<<","<<b1<<","<<b2<<")" <<endl;
                  points.push_back(P2);
                }
            } // end of double loop over b's
      } // end of loop over alphasF4
    if (debug)
      cout << " returning "<<points.size() <<" corners" <<endl;
    return points;
}


// return sublist of a in alist which have t least 3 vertices in points
CuspList remove_redundants(const CuspList& alist, const H3pointList& points)
{
  CuspList new_alist;
  std::copy_if(alist.begin(), alist.end(), std::back_inserter(new_alist),
               [points](const RatQuad& a) {return nverts(a, points) >= 3;});
  return new_alist;
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

CuspList saturate_covering_alphas(const CuspList& alphas, const CuspList& slist, INT maxn, int debug, int verbose)
{
  if (verbose)
    {
      cout << "Saturating "<<alphas.size()<<" alphas with max dnorm "<< max_dnorm(alphas) <<endl;
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
  for ( const auto& s : slist)
    {
      H3point P = {s, RAT(0)};
      pointsx.push_back(P);
      for ( const auto& t : translates)
        {
          pointsx.push_back(translate(P,t));
          pointsx.push_back(translate(P,-t));
        }
    }
  CuspList new_alphas = remove_redundants(alphas, pointsx);
  if (verbose)
    cout << "After removing alphas which go through <3 vertices, we now have "
         <<new_alphas.size()<<" alphas with max norm "<< max_dnorm(new_alphas) <<endl;

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
                                  [maxn](H3point P) { return maxn*P.t2>=ONE;}),
                   points.end());
      if (debug)
        {
          cout << " -- of which "<<points.size()<<" are low enough to be properly covered by a new alpha"<<endl;
          // for (const auto& P : points)
          //   cout << "P = " << P << " is in first quadrant? "<< P.z.in_quarter_rectangle() << endl;
        }
      points.erase(std::remove_if(points.begin(), points.end(),
                                  [](H3point P) { return !P.z.in_quarter_rectangle();}),
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
                             [](const RatQuad& a) {return a.den().norm();});
              cout << " with denominator norms " << norms;
              cout << ";  heights above P (height "<<P.t2<<") are ";
              vector<RAT> hts;
              hts.resize(extras.size());
              std::transform(extras.begin(), extras.end(), hts.begin(),
                             [P](const RatQuad& a) {return height_above(a,P.z);});
              cout <<hts <<endl;
            }
          for ( auto& a : extras)
            {
              RatQuad ca = a.conj();
              CuspList blist = {a,-a,ca,-ca};
              for ( auto& b : blist)
                {
                  b = reduce_to_rectangle(b);
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
              cout << " alphas are saturated! "<<new_alphas.size()<<" alphas with max norm "
                   <<max_dnorm(new_alphas)<<endl;
            }
          else
            {
              cout << " alphas not saturated, "<<extra_alphas.size()<<" extras needed: "
                   <<extra_alphas<<" (with norms at most "<<max_dnorm(extra_alphas)<<")"<<endl;
            }
        }
      new_alphas.insert(new_alphas.end(), extra_alphas.begin(), extra_alphas.end());
    } // ends while(!sat)

  if (verbose)
    cout << "After saturation we now have "<<new_alphas.size()<<" alphas with max norm "<<max_dnorm(new_alphas)<<endl;

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
  for ( const auto& s : slist)
    {
      H3point P = {s, RAT(0)};
      pointsx.push_back(P);
      for ( const auto& t : translates)
        {
          pointsx.push_back(translate(P,t));
          pointsx.push_back(translate(P,-t));
        }
    }
  new_alphas = remove_redundants(new_alphas, pointsx);
  if (verbose)
    cout << "After removing alphas which go through <3 vertices, we now have "
         <<new_alphas.size()<<" alphas with max norm "<< max_dnorm(new_alphas) <<endl;
  std::sort(new_alphas.begin(), new_alphas.end());
  return new_alphas;
}

// return  a saturated irredundant list of alphas in the fundamental rectangle
CuspList find_alphas(const CuspList& slist, int debug, int verbose)
{
  auto alist = covering_alphas(slist, verbose);
  INT maxn = max_dnorm(alist);
  return saturate_covering_alphas(alist, slist, maxn, debug, verbose);
}

// return list of alphas (or translates) which pass through a finite cusp
CuspList neighbours(const RatQuad& a, const CuspList& alist)
{
  Quad r = a.num(), s=a.den();
  INT ns = s.norm();
  CuspList ans;
  const std::array<int,3> txy = {-1,0,1};
  for ( const auto& alpha : alist)
    {
      Quad c = alpha.num(), d=alpha.den();
      Quad f = r*d-s*c, g=s*d;
      // S_alpha goes through a iff N(r*d-s*c)=N(s)
      for (auto tx : txy)
        for (auto ty : txy)
          {
            Quad t(tx,ty);    // test alpha+t = RatQuad(c+dt,d)
            Quad b = f-g*t; // = rd-s(c+dt)
            if (b.norm()==ns)
              ans.push_back(RatQuad(c+d*t,d, 0)); // 0 means do not reduce
          }
    }
  return ans;
}

// return sorted list of alphas (or translates) which pass through a finite cusp,
// i.e. angle_under_pi(sigma, a[i-1], a[i]) for all 0<=i<n and angle_under_pi(sigma, a[n-1], a[0]).
CuspList sorted_neighbours(const RatQuad& sigma, const CuspList& alphas)
{
  return circular_sort(sigma, neighbours(sigma, alphas));
}

// test if one singular point (sigma) is surrounded by alpha circles:
int is_sigma_surrounded(const RatQuad& sigma, const CuspList& alphas, int debug)
{
  if (sigma.is_infinity())
    return 1;
  CuspList ans = neighbours(sigma, alphas);
  if (debug)
    cout<<"Neighbours of "<<sigma<<" are "<<ans<<endl;

  int ok = std::all_of(ans.begin(), ans.end(),
                     [sigma, ans](const RatQuad& a1)
                     {return std::any_of(ans.begin(), ans.end(),
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
int are_sigmas_surrounded(const CuspList& slist, const CuspList& alphas, int debug)
{
  return std::all_of(slist.begin(), slist.end(),
                     [alphas,debug](const RatQuad& s) {return is_sigma_surrounded(s, alphas, debug);});
}

