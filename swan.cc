// FILE SWAN.CC: implementation of Swan's algorithm

#include "swan.h"
#include "swan_sigmas.h"
#include "swan_alphas.h"

void SwanData::make_sigmas() {
  if (slist.empty())
    {
      string step = "SwanData::make_sigmas()";
      SwanTimer.start(step);
      slist = sort_singular_points(singular_points());
      slistx.clear();
      for (const auto& s : slist)
        {
          if (s.is_infinity())
            continue;
          CuspList s_sh = cusp_shifts(s, Quad::shifts_by_one);
          slistx.insert(slistx.end(), s_sh.begin(), s_sh.end());
        }
      SwanTimer.stop(step);
      if (showtimes) SwanTimer.show(1, step);
    }
}

void SwanData::make_alphas(int verbose) {
  if (alist.empty())
    {
      string step = "SwanData::make_alphas()";
      SwanTimer.start(step);
      make_sigmas();
      find_covering_alphas(verbose);
      saturate_alphas(verbose);
      SwanTimer.stop(step);
      if (showtimes) SwanTimer.show(1, step);
    }
}

// add one alpha
int SwanData::add_one_alpha(const RatQuad& a, int covered, int verbose)
{
  if (circle_inside_any_circle(a, alistx, 0)) // strict=0
    return 0;

  alist.push_back(a);
  auto a_sh = cusp_shifts(a, Quad::shifts_by_one);
  alistx.insert(alistx.end(), a_sh.begin(), a_sh.end());

  Quad u; RatQuad b0, at;
  if (verbose)
    cout<<"appending a="<<a<<" to alist and its translates to alistx"<<endl;
  nbrs[a] = intersecting_alphas(a,alistx);
  if (verbose)
    cout<<" - nbrs of "<<a<<" initially "<<nbrs[a]<<endl;
  if (a.in_quarter_rectangle())
    {
      alist_open.push_back(a);
      alistF4.push_back(a);
      if (verbose)
        cout<<"appending a="<<a<<" to alistF4"<<endl;
      if (covered)
        {
          nbrs_ok[a] = nbrs[a];
          nbrs_open[a] = CuspList();
        }
      else
        {
          nbrs_open[a] = nbrs[a];
          nbrs_ok[a] = CuspList();
        }
    }
  else
    {
      alist_ok.push_back(a);
    }
  // Update nbr lists of earlier alphas (in F4) if they intersect
  // one of the new translates:
  for (const auto& b : nbrs[a]) // NB these may have been translated
    {
      b0 = reduce_to_rectangle(b, u);
      // now b0+u intersects a, so a-u intersects b0
      if (b0.in_quarter_rectangle())
        {
          at = a-u;
          nbrs[b0].push_back(at);
          // Add at to nbrs_open[b0], unless nbrs_open[b0] is empty,
          // meaning that b0 is already surrounded, in which case we
          // add at to nbrs_ok[b0]:
          if (nbrs_open[b0].empty())
            {
              if (std::find(nbrs_ok[b0].begin(), nbrs_ok[b0].end(), at) == nbrs_ok[b0].end())
                {
                  nbrs_ok[b0].push_back(at);
                  if (verbose)
                    cout<<"added nbr "<<at<<" = "<<a<<"-"<<u<<" to "<<b0<<"'s ok list"<<endl;
                }
            }
          else
            {
              if (std::find(nbrs_open[b0].begin(), nbrs_open[b0].end(), at) == nbrs_open[b0].end())
                {
                  nbrs_open[b0].push_back(at);
                  if (verbose)
                    cout<<"added nbr "<<at<<" = "<<a<<"-"<<u<<" to "<<b0<<"'s open list"<<endl;
                }
            }
        }
    }
  return 1;
}

int SwanData::add_new_alphas(int verbose)
{
  INT first_maxn = 1+Quad::absdisc/4; // we first consider all alphas with dnorm up to this
  int first = maxn==0;
  string s = (first? " of norm up to " : " of norm ");

  // Get the next batch new_alphas, either of all dnorms up to maxn, or those with the next dnorm
  CuspList new_alphas = (first
                         ?
                         principal_cusps_with_denominators(denom_looper.values_with_norm_up_to(first_maxn, 1))
                         :
                         principal_cusps_with_denominators(denom_looper.values_with_current_norm())
                         );
  maxn = new_alphas.back().den().norm();

  if (verbose)
    {
      cout << "Considering "<<new_alphas.size()<<(first?"":" extra") << " principal cusps " << s << maxn
           << endl;
      // cout << new_alphas <<endl;
    }

  // Of these, check for redundancy; if not redundant, append it to
  // alist, and it and its 8 translates to alistx.

  // If it is in the quarter rectangle F4, append it to alistF4 and
  // alist_open; find its intersecting neighbours from alistx and add
  // this to the map nbrs_open.

  int nc = 0, nc4 = 0; // number of new alphas and of those in F4

  for ( const auto& a : new_alphas)
    {
      if (add_one_alpha(a, 0, verbose>1))
        {
          nc++;
          if (a.in_quarter_rectangle())
            nc4 ++;
        }
    }

  // report on which new alphas we will be testing, if any:

  if (verbose)
    {
      if (nc)
        {
          cout << "Adding "<<nc<<" alphas " << s << maxn;
          cout << ", of which " <<nc4 << " are in the quarter rectangle; #alphas="<<alist.size();
          cout << ", with "<<alist_open.size()<<" not yet shown to be surrounded"<<endl;
        }
      else
        {
          cout << "All in this batch are redundant\n";
        }
    }

  return nc;
}

// Sets alist to a list of principal cusps a such that the S_{a+t}
// for all integral t cover CC apart from singular points.

// For n>=1 successively, we test as a candidate set all a=r/s with
// r,s coprime, r reduced mod s, N(s)<=n (omitting any for which S_a
// is contained in any earlier S_a') until we succeed.

// Other methods then saturate the set and discard redundancies.

void SwanData::find_covering_alphas(int verbose)
{
  string step = "SwanData::find_covering_alphas()";
  SwanTimer.start(step);
  int ok=0;
  while (!ok)
    {
      // Get the next batch new_alphas, either of all dnorms up to maxn, or those with the next dnorm
      int nc = add_new_alphas(verbose);
      if (nc==0)
        continue;

      // Test whether all alist_open are now surrounded (by alistx).
      // As a side effect, some alphas will be moved from alist_open
      // to alist_ok, and the nbrs_open and nbrs_ok lists are updated.

      ok = are_alphas_surrounded(verbose);
      if (!ok &&verbose)
        cout << "Some alphas are not surrounded, continuing...\n";
    }

  alist = alist_ok;
  if (verbose)
    cout << "Success in covering using "<<alist.size()<<" alphas of with max norm "<<maxn<<"\n";
  SwanTimer.stop(step);
  if (showtimes) SwanTimer.show(1, step);
}

// test if a is singular by reducing to rectangle and comparing with
// slist (but there are two special cases)
int SwanData::is_singular(const RatQuad& a)
{
  Quad t;
  RatQuad a0 = reduce_to_rectangle(a,t);
  if (std::find(slist.begin(), slist.end(), a0) != slist.end())
    return 1;

  // annoying special cases for historical back-compatibility.
  long d = Quad::d;
  // When d%8==7, we use s=(1-w)/2 but (w-1)/2 is in the rectangle
  if (d%8==7 && a0==RatQuad(Quad(-1,1),Quad(2))) return 1;
  // When d%12==15 and d>15, we use s=(-1-w)/3 but (2-w)/3 is in the rectangle
  if (d%12==3 && d>15 && a0==RatQuad(Quad(2,-1),Quad(3))) return 1;
  return 0;
}

// Test if intersection points of S_a and S_b are covered by some S_c.
// Assume that a is in F4 so we know its list of intersecting
// neighbours, that b is in that list and any covering S_c is also in
// that list
int SwanData::are_intersection_points_covered(const RatQuad& a, const RatQuad& b, int verbose)
{
  int debug=verbose;
  if (debug)
    cout << "Testing if intersection points of "<<a<<" and "<<b<<" are covered"<<endl;

  const CuspList& a_nbrs = nbrs[a];

  // First see if the intersection points are in k, in which case we
  // test whether they are singular and otherwise have an easy test:
  CuspList zlist = intersection_points_in_k(a,b);
  if (!zlist.empty())
    {
      if (debug)
        cout << " intersection points are k-rational: "<<zlist<<endl;
      for ( const auto& z : zlist)
        {
          if (is_singular(z)) // OK
            {
              if (debug)
                cout << " ok, "<<z<<" is singular"<<endl;
              continue;
            }
          if (!is_inside_one(z, a_nbrs, 1 /*strict*/))  // z is not covered: not OK
            {
              if (debug)
                cout << " returning no, "<<z<<" is not covered"<<endl;
              return 0;
            }
        }
      // we reach here if both z values are either singular or covered: OK
      if (debug)
        cout << "+++returning yes, both are covered"<<endl;
      return 1;
    }

  // Now the intersection points are not in k. Check that either one
  // S_a covers both, or two cover one each:

  int t = 0; // will hold +1 or -1 if we have covered only one of the two
  for ( const auto& c : a_nbrs)
    {
      if (c == b)
        continue;
      if (!circles_intersect(b,c))
        continue;
      // call the global function: it returns 0 for neither, 2 for
      // both and +1,-1 consistently for just one
      int t2 = are_intersection_points_covered_by_one(a, b, c);
      if (t2==2) // both are covered by c
        {
          if (debug)
            cout << "+++returning yes, both are covered"<<endl;
          return 1;
        }
      if (t2==0) // neither is covered by c
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

int SwanData::is_alpha_surrounded(const RatQuad& a, int verbose)
{
  int debug = verbose>1;
  if (debug)
    cout<<"Testing if "<<a<<" is surrounded..."<<flush;

  CuspList& a_nbrs_open = nbrs_open[a]; // a reference, so changes are kept
  CuspList& a_nbrs_ok = nbrs_ok[a];     // a reference, so changes are kept
  if (nbrs_open[a].empty())
    {
      if (debug)
        cout<<"yes, already shown"<<endl;
      return 1;
    }

  int all_ok = 1;
  for ( auto it = a_nbrs_open.begin(); it!=a_nbrs_open.end(); )
    {
      RatQuad b = *it;
      if (are_intersection_points_covered(a, b, debug))
        {
          // record that this pair is ok
          a_nbrs_ok.push_back(b);
          it = a_nbrs_open.erase(it);
          if (debug)
            {
              cout<<" - intersection points of S_{"<<a<<"} and S_{"<<b<<"} are surrounded"<<endl;
              // cout<<" - alist_open["<<a<<"] = "<<a_nbrs_open<<endl;
              // cout<<" - alist_ok["<<a<<"] = "<<a_nbrs_ok<<endl;
            }
        }
      else
        {
          it++;
          if (debug)
            cout<<" - intersection points of S_{"<<a<<"} and S_{"<<b<<"} are NOT surrounded"<<endl;
          all_ok = 0; // but continue checking the other bs
        }
    }
  if (debug)
    {
      if (all_ok)
        {
          cout << " S_{"<<a<<"} IS surrounded"<<endl;
          cout<<" - nbrs_open["<<a<<"] = "<<nbrs_open[a]<<endl;
          cout<<" - nbrs_ok["<<a<<"] = "<<nbrs_ok[a]<<endl;
        }
      else
        cout << " S_{"<<a<<"} is NOT surrounded"<<endl;
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

// We assume that all a in alist_open are in the quarter rectangle

// Returns 0/1

// NB Not all a in alist_open will be tested as we quit and return 0
// after any failure

int SwanData::are_alphas_surrounded(int verbose)
{
  int debug = verbose>1;
  if (debug)
    {
      cout << "At start of are_alphas_surrounded()" << endl;
      // cout << "alist_ok = "<<alist_ok<<endl;
      // cout << "alist_open = "<<alist_open<<endl;
    }
  int i=0, n_open = alist_open.size();

  // We make a copy of alist_open to loop over, so we can delete elements from the original as we go
  CuspList new_alist_open = alist_open;
  for ( const auto& a : new_alist_open)
    {
      i++;
      if (verbose) cout <<"Testing alpha #"<<i<<"/"<<n_open<<" = "<<a<<"...";
      if (is_alpha_surrounded(a, verbose))
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
  int ok = are_sigmas_surrounded(debug);
  if (verbose)
    {
      if (ok)
        cout << " and singular points are surrounded!" <<endl;
      else
        cout << " but singular points are not yet surrounded" <<endl;
    }
  return ok;
}

// test if one singular point (sigma) is surrounded by alpha circles:
int SwanData::is_sigma_surrounded(const RatQuad& s, int verbose)
{
  if (s.is_infinity())
    return 1;
  CuspList ans = neighbours(s, alistx);
  if (verbose)
    cout<<"Neighbours of "<<s<<" are "<<ans<<endl;

  int ok = std::all_of(ans.begin(), ans.end(),
                       [s, ans](const RatQuad& a1)
                       {return std::any_of(ans.begin(), ans.end(),
                                           [s,a1](const RatQuad& a2) {return angle_under_pi(s,a1,a2);});});
  if (verbose)
    {
      if (ok)
        cout<<" + "<<s<<" IS surrounded"<<endl;
      else
        cout<<" - "<<s<<" is NOT surrounded"<<endl;
    }
  return ok;
}

int SwanData::are_sigmas_surrounded(int verbose)
{
  make_sigmas();
  return std::all_of(slist.begin(), slist.end(),
                     [verbose, this](const RatQuad& s) {return is_sigma_surrounded(s, verbose);});
}

// list of singular corners [s,0] on S_a (s in slist or a translate)
H3pointList SwanData::singular_corners(const RatQuad& a)
{
  H3pointList ans;
  for (const auto& s : slistx)
    {
      if (is_on(s, a))
        ans.push_back({s,RAT(0)});
    }
  return ans;
}

// Find corners from one alpha.  The new corners are not added to the
// class's points list, but are returned.  The parameter redundant is
// set to 1 if a has <3 corners (including singular ones).
H3pointList SwanData::find_corners_from_one(const RatQuad& a, int& redundant, int verbose)
{
  string step = "SwanData::find_corners_from_one()";
  SwanTimer.start(step);
  int debug = verbose;
  Quad t;
  H3pointList new_corners; // list of corners to be returned (and then added to class global corners)
  H3pointList a_corners = singular_corners(a); // list of all corners of a found

  // these two are different: new_corners includes symmetrics and are
  // all in the rectangle (suitable for merging onto the class's
  // corners list), while a_corners will be a list of all corners of
  // a, without symmetrics and not all necessarily in the
  // rectangle. We only need to know the size of a_corners as the
  // condition for a to be redundant is that it has <3 corners.

  // get its intersecting neighbours
  if (!nbrs_open[a].empty())
    {
      cout<<"Problem: "<<a<<" has open neighbours "<<nbrs_open[a]<<endl;
      cout<<"find_corners_from_one() is only to be used after finding a covering set of alphas"<<endl;
      return new_corners; // empty
    }
  CuspList a_nbrs = nbrs_ok[a]; // which = nbrs[a]
  int n = a_nbrs.size();

  if (debug)
    cout<<"-------------------------\n"
        <<"Finding corners of "<<a<<" from its "<<n<<" neighbours "<<a_nbrs
        <<"\n it has "<<a_corners.size()<<" singular corners: "<<a_corners<<endl;

  // For each intersecting pair of these we may have a triple intersection point:
  if (n<3)
    return new_corners;
  for (auto bi = a_nbrs.begin(); bi!=a_nbrs.end(); ++bi)
    for (auto bj = bi+1; bj!=a_nbrs.end(); ++bj)
      {
        RatQuad b1 = *bi, b2 = *bj;
        if (!circles_intersect(b1,b2))
          continue;
        H3pointList corners1 = tri_inter_points(a, b1, b2);
        if (corners1.empty())
          continue;
        H3point P = corners1.front();
        RatQuad z = P.z;
        RAT t2 = P.t2;
        if (t2.sign()==0)
          continue;
        if (debug)
          cout << " ----------found P = "<<P<<", tri-intersection of "<<a<<", "<<b1<<", "<<b2<<"\n";
        if (is_under_any(P, alistx))
          {
            if (debug)
              cout << " ignoring P as it is under some hemisphere"<<endl;
            continue;
          }
        if (debug)
          cout << " P is not under any hemisphere"<<endl;
        if (std::find(a_corners.begin(), a_corners.end(), P) == a_corners.end())
          {
            a_corners.push_back(P);
            if (debug)
              cout<<" adding P to list of corners of "<<a
                  <<" (which now contains "<<a_corners.size()<<" corners)"<<endl;
          }
        else
          if (debug)
                cout<<" P is already in list of corners of "<<a<<endl;
        // Now see if we want to add P (and its symmetrics) to the new_corners list.
        // Not if it is there already:
        if (std::find(new_corners.begin(), new_corners.end(), P) != new_corners.end())
          {
            if (debug)
              cout << " P is already in new_corners list"<<endl;
            continue;
          }
        if (debug)
          cout << " P is not in new_corners list"<<endl;

        // Store this if it is in F4 (then store up to 3 symmetrics):
        if (!z.in_quarter_rectangle())
          {
            if (debug)
              cout << " P is not in F4, not adding to new_corners list"<<endl;
            continue;
          }
        if (debug)
          cout << " P is in F4, adding to corners list"<<endl;
        new_corners.push_back(P);

        // These corners are in F4, so we apply symmetries to get all those in F:
        RatQuad zbar = z.conj();
        for (auto z2 : {-z, zbar, -zbar})
          {
            z2 = reduce_to_rectangle(z2, t);
            H3point P2 = {z2, t2};
            if (std::find(new_corners.begin(), new_corners.end(), P2) != new_corners.end())
              continue;
            if (is_under_any(P, alistx))
              continue;
            if (debug)
              cout << " adding symmetric "<<P2<<" from ("<<a<<","<<b1<<","<<b2<<") to new_corners list" <<endl;
            new_corners.push_back(P2);
          }
      } // end of double loop over b's
  n = a_corners.size(); // re-using variable n
  redundant = n<3;
  if (debug)
    {
      cout<<"  "<<a<<" has "<<n<<" corners (including singular corners) so ";
      cout<<(redundant?"IS":"is NOT")<<" redundant\n\n";
    }
  SwanTimer.stop(step);
  if (showtimes) SwanTimer.show(1, step);
  return new_corners;
}

// Find potential corners, store in class's corners list
void SwanData::new_find_corners(int verbose)
{
  string step = "SwanData::find_corners()";
  SwanTimer.start(step);
  corners = triple_intersections(alist, verbose);
  H3pointList cornersx;
  for ( const auto& P : corners)
    {
      auto Ps = H3point_shifts(P, Quad::shifts_by_one);
      cornersx.insert(cornersx.end(), Ps.begin(), Ps.end());
    }
  H3pointList s_corners(slistx.size());
  const RAT zero(0);
  std::transform(slistx.begin(), slistx.end(), s_corners.begin(),
                 [zero] (const RatQuad& s) {return H3point({s,zero});});
  cornersx.insert(cornersx.end(), s_corners.begin(), s_corners.end());

  CuspList alist0 = remove_redundants(alist, cornersx);
  if (alist.size() > alist0.size())
    {
      if (verbose)
        cout<<"...number of alphas reduced from "<<alist.size()<<" to "<<alist0.size()<<endl;
      alist = alist0;
    }
  SwanTimer.stop(step);
  if (showtimes) SwanTimer.show(1, step);
}

// Find potential corners, store in class's corners list, replacing
// alistF4 with sublist of alphas in F4 on >=3 corners
void SwanData::old_find_corners(int verbose)
{
  string step = "SwanData::old_find_corners()";
  SwanTimer.start(step);
  int debug = verbose>1;
  if (verbose)
    {
      cout << "In SwanData.find_corners() with " <<alist.size()<<" alphas..."<<endl;
      cout << alistF4.size() <<" alphas are in the quarter rectangle"<<endl;
    }

  CuspList not_redundants; // list of a in alistF4 which are on >=3 corners
  Quad t;

  // Loop through all a in alistF4 finding corners from each:
  for (const auto& a : alistF4)
    {
      int red;
      H3pointList a_corners = find_corners_from_one(a, red, debug);
      if (!red)
        not_redundants.push_back(a);
      for (const auto& P : a_corners)
        if (std::find(corners.begin(), corners.end(), P) == corners.end())
          corners.push_back(P);
    }

  if (verbose)
    {
      cout << " SwanData.find_corners() found "<<corners.size() <<" corners: " << corners<<endl;
      int nnotred = not_redundants.size();
      int nred =  alistF4.size() - nnotred;
      cout << " of " <<alistF4.size() << " alphas in quarter_rectangle, "
           <<nred<<" are redundant and " <<nnotred<<" are not"<<endl;
      cout<<"New alistF4: "<<not_redundants<<endl;
    }

  alistF4 = not_redundants;

  // re-expand to full alist excluding redundants
  // NB We will not need alistx again; if that changes, need to update it here
  auto nalist = alist.size();
  if (verbose>1) cout << "Old alist: "<<alist<<endl;
  alist = alistF4;

  for (const auto& a : alistF4)
    {
      RatQuad abar = a.conj();
      for (auto b : {-a, abar, -abar})
        {
          b = reduce_to_rectangle(b, t);
          if (std::find(alist.begin(), alist.end(), b) == alist.end())
            alist.push_back(b);
        }
    }
  if (verbose)
    {
      if (alist.size()<nalist)
        cout << "Number of alphas in alist reduced from "<<nalist<<" to "<<alist.size()<<endl;
      else
        {
          cout << "Number of alphas in alist remains "<<nalist<<endl;
          if (verbose>1)
            cout << "New alist: "<<alist<<endl;
        }
    }
  SwanTimer.stop(step);
  if (showtimes) SwanTimer.show(1, step);
}

// After an unsuccessful saturation loop which produces extra alphas a
// such that S_a properly covers an old corners (which was then
// deleted), use these to compute more corners.  These will (up to
// symmetry) come from a triple intersection involving at least one of
// the new alphas.  This function just merges the output of
// find_corners_from_one(a).
H3pointList SwanData::find_extra_corners(const CuspList& extra_alphas)
{
  H3pointList extra_corners;
  for (const auto& a : extra_alphas)
    {
      int red;
      H3pointList a_corners = find_corners_from_one(a, red);
      for (const auto& P : a_corners)
        if (std::find(extra_corners.begin(), extra_corners.end(), P) == extra_corners.end())
          extra_corners.push_back(P);
    }
  return extra_corners;
}

void SwanData::saturate_alphas(int verbose)
{
  string step = "SwanData::saturate_alphas()";
  SwanTimer.start(step);
  int debug = verbose>0;
  int already_saturated = 1; // will flip to 0 if any new alphas are added
  auto nalist = alist.size();
  INT m = max_dnorm(alist);
  Quad temp;
  if (verbose)
    cout << "Saturating "<<nalist<<" alphas with max dnorm "<< m <<endl;

  // Find triple intersections, reducing alist and alistF4 to exclude
  // any redundants (on <3 corners):

  find_corners(verbose);
  m = max_dnorm(alist);
  if (verbose)
    {
      cout << "Found "<<corners.size() << " corners (triple intersection points)" <<endl;
      if (alist.size()<nalist)
        cout << " (number of alphas reduced from "<<nalist<<" to "<<alist.size() << " with max dnorm " <<m<<")"<<endl;
      else
        cout << " (number of alphas remains "<<nalist<<")"<<endl;
    }

  int sat = 0, first_run=1;
  H3pointList corners_open = corners;
  CuspList extra_alphas; // will be filled with any extra alphas needed on each pass
  while (!sat)
    {
      if (!first_run)
        {
          corners_open = find_extra_corners(extra_alphas);
          if (verbose)
            {
              cout << "Found "<<corners_open.size()<<" new potential corners"<<endl;
              cout << " - number of corners was "<<corners_open.size()<<endl;
            }
          for (const auto& P : corners_open)
            if (std::find(corners.begin(), corners.end(), P) == corners.end())
              corners.push_back(P);
          if (verbose)
            cout << " - number of corners is now "<<corners_open.size()<<endl;
        }
      if (first_run)
        {
          if(debug)
            cout << "Extracting those of square height less than 1/"<<maxn<<endl;
          // Remove corners which cannot be better covered by an alpha with dnorm>maxn
          corners_open.erase(std::remove_if(corners_open.begin(), corners_open.end(),
                                            [this](H3point P) { return maxn*P.t2>=ONE;}),
                             corners_open.end());
          if (debug)
            cout << " -- of which "<<corners_open.size()
                 <<" are low enough to be properly covered by a new alpha"<<endl;
        }

      // only check corners in F4 (if we find new alphas we'll add
      // their symmetrics too):
      corners_open.erase(std::remove_if(corners_open.begin(), corners_open.end(),
                                        [](H3point P) { return !P.z.in_quarter_rectangle();}),
                         corners_open.end());
      if (debug)
        cout << " -- of which "<<corners_open.size()<<" lie in the first quadrant" << endl;

      first_run = 0;
      extra_alphas.clear();
      sat = 1;        // will be set to 0 if we find out that the alphas are not already saturated
      int iP = 0, nP = corners_open.size();
      for (const auto& P : corners_open)
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
              continue; // on to the next P
            }
          sat = 0; // we now know the alphas are not saturated
          already_saturated = 0;
          if (debug)
            {
              auto& a = extras.front();
              cout << "   - found " <<extras.size()<<" maximally properly covering hemispheres"
                   << " with denominator norm " << a.den().norm()
                   << ", height " << height_above(a,P.z) << " above P (height "<<P.t2<<")" << endl;
            }
          // for each a in extras, add it and its 4 symmetrics to extra_alphas
          for ( auto& a : extras)
            {
              RatQuad ca = a.conj();
              CuspList blist = {a,-a,ca,-ca};
              for ( auto& b : blist)
                {
                  b = reduce_to_rectangle(b, temp);
                  if (debug)
                    cout << " - testing potential new alpha " << b << endl;
                  if (std::find(extra_alphas.begin(), extra_alphas.end(), b) == extra_alphas.end())
                    {
                      if (debug)
                        cout << " - adding alpha " << b << endl;
                      extra_alphas.push_back(b);
                    }
                  else
                    if (debug)
                      cout << " - not adding it, it's a repeat" <<endl;
                } // loop over 4 flips of a
            } // loop over a in extras

          //  Delete this P (and symmetrics) from corners:
          if (verbose)
            cout<<"Deleting "<<P<<" from corners list"<<endl;

          Quad t;
          RAT t2 = P.t2;
          RatQuad z = P.z;
          RatQuad zbar = z.conj();
          for (auto z2 : {z, -z, zbar, -zbar})
            {
              H3point Q = {reduce_to_rectangle(z2, t), t2};
              auto it = std::find(corners.begin(), corners.end(), Q);
              if ( it != corners.end())
                corners.erase(it);
            }

        } // checked all P in corners_open

      // Now we have finished checking all potential corners; any
      // which could be properly covered have been deleted and some
      // new alphas collected in extra_alphas; if there were any such,
      // sat has been set to 0. ALSO some of the original alphas not
      // thought to be redundant may now be so, since we have deleted
      // some corners. Deal with this after finishing saturation.
      if (verbose)
        {
          if (sat)
            {
              m = max_dnorm(alist);
              cout << " alphas are saturated! "<<alist.size()<<" alphas with max norm "<<m<<endl;
            }
          else
            {
              m = max_dnorm(extra_alphas);
              cout << " alphas not saturated, "<<extra_alphas.size()<<" extras needed: "
                   <<extra_alphas<<" with norms at most "<<m<<endl;
            }
        }
      if (sat)
        break; // out of the while(!sat) loop
      // add the extra alphas to alist and alist4, updating the nbrs
      // lists:
      for (const auto& a : extra_alphas)
        add_one_alpha(a, 1, verbose>1);
      if (verbose)
        cout<<"alist now has size "<<alist.size()<<endl;
    } // ends while(!sat)

  if (already_saturated) // the the original alist was saturated
    {
      if (verbose)
        cout<<"Already saturated: we have "<< alist.size()<<" alphas with max norm "<<m<<endl;
    }
  else
    {
      if (verbose)
        cout<<"Final check for alphas now redundant..."<<endl;
      // recheck that all remaining alphas are not redundant, which
      // they might have become on deleting some potential corners:
      H3pointList cornersx;
      for ( const auto& P : corners)
        {
          auto Ps = H3point_shifts(P, Quad::shifts_by_one);
          cornersx.insert(cornersx.end(), Ps.begin(), Ps.end());
        }
      H3pointList s_corners(slistx.size());
      RAT zero(0);
      std::transform(slistx.begin(), slistx.end(), s_corners.begin(),
                     [zero] (const RatQuad& s) {return H3point({s,zero});});
      cornersx.insert(cornersx.end(), s_corners.begin(), s_corners.end());

      CuspList alist0 = remove_redundants(alist, cornersx);
      if (alist.size() > alist0.size())
        {
          if (verbose)
            cout<<"...number of alphas reduced from "<<alist.size()<<" to "<<alist0.size()<<endl;
          alist = alist0;
        }
      else
        if (verbose)
          cout<<"...number of alphas reduced from "<<alist.size()<<" to "<<alist0.size()<<endl;
    }

  // (re)sort alist before returning:
  std::sort(alist.begin(), alist.end(), Cusp_cmp);
  SwanTimer.stop(step);
  if (showtimes) SwanTimer.show(1, step);
}
