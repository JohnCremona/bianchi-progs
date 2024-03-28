// FILE SWAN.CC: implementation of Swan's algorithm

#include "swan.h"
#include "swan_sigmas.h"
#include "swan_alphas.h"

SwanData::SwanData() {
  maxn = 0;
  for (auto x : {-1,0,1}) for (auto y : {-1,0,1}) shifts.push_back(Quad(x,y));
}

void SwanData::make_sigmas() {
  if (slist.empty())
    {
      slist = sort_singular_points(singular_points());
      slistx.clear();
      for (const auto& s : slist)
        {
          if (s.is_finite())
            {
              for (const auto& t : shifts)
                slistx.push_back(s+t);
            }
        }
    }
}

void SwanData::make_alphas(int verbose) {
  if (alist.empty())
    {
      make_sigmas();
      find_covering_alphas(verbose);
      saturate_alphas(verbose);
    }
}

// add one alpha
int SwanData::add_one_alpha(const RatQuad& a, int covered, int verbose)
{
  if (circle_inside_any_circle(a, alistx, 0)) // strict=0
    return 0;

  alist.push_back(a);
  for (const auto& t : shifts)
    alistx.push_back(a+t);

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
  int nc, ok=0;
  while (!ok)
    {
      // Get the next batch new_alphas, either of all dnorms up to maxn, or those with the next dnorm
      nc = add_new_alphas(verbose);
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
  if (d%8==7 && a0==RatQuad(Quad(-1,1),TWO)) return 1;
  // When d%12==15 and d>15, we use s=(-1-w)/3 but (2-w)/3 is in the rectangle
  if (d%12==3 && d>15 && a0==RatQuad(Quad(2,-1),THREE)) return 1;
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
      //if (::are_intersection_points_covered(a, b, alistx, slist, debug))
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
        ans.push_back({s,0});
    }
  return ans;
}

// Find and fill corners list, replacing alistF4 with sublist of alphas in F4 on >=3 corners
void SwanData::find_corners(int verbose)
{
  int debug = verbose>1;
  if (verbose)
    cout << "In SwanData.find_corners() with " <<alist.size()<<" alphas..."<<endl;

  // Extract the alphas in F4.  NB the standard list of alphas has
  // the property that, after those of denom 1, 2, 3, they come in
  // pairs a, -a; the second of the pair is not always in the
  // rectangle (because of the rounding of 1/2).
  Quad t;
  if (verbose)
    cout << alistF4.size() <<" alphas are in the quarter rectangle"<<endl;

  CuspList not_redundants; // list of a in alistF4 which are on >=3 corners

  // Loop through all a0 in alistF4:
  for (const auto& a0 : alistF4)
    {
      // get its intersecting neighbours
      if (!nbrs_open[a0].empty())
        {
          cout<<"Problem: "<<a0<<" has open neighbours "<<nbrs_open[a0]<<endl;
          exit(1);
        }
      CuspList a_nbrs = nbrs_ok[a0];
      H3pointList a0_corners = singular_corners(a0); // list of all corners of a0 found

      if (debug)
        cout<<"-------------------------\n"
            <<"Finding corners of "<<a0<<" from its "<<a_nbrs.size()<<" neighbours "<<a_nbrs
            <<"\n it has "<<a0_corners.size()<<" singular corners: "<<a0_corners<<endl;

      // For each intersecting pair of these we may have a triple
      // intersection point:
      int n = a_nbrs.size();
      if (n<3)
        return; // empty list
      for (auto bi = a_nbrs.begin(); bi!=a_nbrs.end(); ++bi)
        for (auto bj = bi+1; bj!=a_nbrs.end(); ++bj)
          {
            RatQuad b1 = *bi, b2 = *bj;
            if (!circles_intersect(b1,b2))
              continue;
            H3pointList corners1 = tri_inter_points(a0, b1, b2);
            if (corners1.empty())
              continue;
            H3point P = corners1.front();
            RatQuad z = P.z;
            RAT t2 = P.t2;
            if (t2.sign()==0)
              continue;
            if (debug)
              cout << " ----------found P = "<<P<<", tri-intersection of "<<a0<<", "<<b1<<", "<<b2<<"\n";
            if (is_under_any(P, alistx))
              {
                if (debug)
                  cout << " ignoring P as it is under some hemisphere"<<endl;
                continue;
              }
            if (debug)
              cout << " P is not under any hemisphere"<<endl;
            if (std::find(a0_corners.begin(), a0_corners.end(), P) == a0_corners.end())
              {
                a0_corners.push_back(P);
                if (debug)
                  cout<<" adding P to list of corners of "<<a0
                      <<" (which now contains "<<a0_corners.size()<<" corners)"<<endl;
              }
            else
              if (debug)
                cout<<" P is already in list of corners of "<<a0<<endl;
            // Now see if we want to add P (and its symmetrics) to the main corners list.
            // Not if it is there already:
            if (std::find(corners.begin(), corners.end(), P) != corners.end())
              {
                if (debug)
                  cout << " P is already in corners list"<<endl;
                continue;
              }
            if (debug)
              cout << " P is not in corners list"<<endl;

            // Store this if it is in F4 (then store up to 3 symmetrics):
            if (!z.in_quarter_rectangle())
              {
                if (debug)
                  cout << " P is not in F4, not adding to corners list"<<endl;
                continue;
              }
            if (debug)
              cout << " P is in F4, adding to corners list"<<endl;
            corners.push_back(P);
            // These corners are in F4, so we apply symmetries to get all those in F:
            RatQuad zbar = z.conj();
            for (auto z2 : {-z, zbar, -zbar})
              {
                z2 = reduce_to_rectangle(z2, t);
                H3point P2 = {z2, t2};
                if (std::find(corners.begin(), corners.end(), P2) != corners.end())
                  continue;
                if (is_under_any(P, alistx))
                  continue;
                if (debug)
                  cout << " adding symmetric "<<P2<<" from ("<<a0<<","<<b1<<","<<b2<<") to corners list" <<endl;
                corners.push_back(P2);
              }
          } // end of double loop over b's
      if (debug)
        cout<<"  "<<a0<<" has "<<a0_corners.size()<<" corners\n\n";
      if (a0_corners.size()>2)
        not_redundants.push_back(a0);
    } // end of loop over alistF4
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
  if (verbose>1)
    {
      cout << "Old alist: "<<alist<<endl;
    }
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
}

void SwanData::saturate_alphas(int verbose)
{
  int debug = verbose>0;
  int already_saturated = 1; // will flip to 0 if any new alphas are added
  auto nalist = alist.size();
  INT m = max_dnorm(alist);
  Quad temp;
  if (verbose)
    {
      cout << "Saturating "<<nalist<<" alphas with max dnorm "<< m <<endl;
    }

  // Find triple intersections, reducing alist and alistF4 to exclude
  // any redundants (on <3 corners):

  find_corners(verbose);
  m = max_dnorm(alist);
  if (verbose)
    {
      cout << "Found "<<corners.size() << " triple intersection points" <<endl;
      if (alist.size()<nalist)
        cout << " (number of alphas reduced from "<<nalist<<" to "<<alist.size() << " with max dnorm " <<m<<")"<<endl;
      else
        cout << " (number of alphas remains "<<nalist<<")"<<endl;
    }

  // add translates of these corners and singular points:

  // cornersx.clear();
  // for ( const auto& P : corners)
  //   for ( const auto& t : shifts)
  //     cornersx.push_back(translate(P,t));
  // for ( const auto& s : slist)
  //   {
  //     H3point P = {s, ZERO};
  //     for ( const auto& t : shifts)
  //       cornersx.push_back(translate(P,t));
  //   }

  int sat = 0, first_run=1;
  H3pointList corners_open = corners, corners_ok;
  while (!sat)
    {
      if (!first_run)
        {
          find_corners(debug);
          corners_open = corners;
          if (verbose)
            cout << "Found "<<corners_open.size()<<" potential vertices"<<endl;
        }
      first_run = 0;
      if(debug)
        cout << "Extracting those of square height less than 1/"<<maxn<<endl;
      // Remove corners which cannot be better covered by an alpha with dnorm>maxn
      corners_open.erase(std::remove_if(corners_open.begin(), corners_open.end(),
                                        [this](H3point P) { return maxn*P.t2>=ONE;}),
                         corners_open.end());
      if (debug)
        cout << " -- of which "<<corners_open.size()<<" are low enough to be properly covered by a new alpha"<<endl;

      corners_open.erase(std::remove_if(corners_open.begin(), corners_open.end(),
                                        [](H3point P) { return !P.z.in_quarter_rectangle();}),
                         corners_open.end());
      if (debug)
        cout << " -- of which "<<corners_open.size()<<" lie in the first quadrant" << endl;

      corners_open.erase(std::remove_if(corners_open.begin(), corners_open.end(),
                                        [corners_ok](H3point P)
                                        {return std::find(corners_ok.begin(), corners_ok.end(), P) != corners_ok.end();}),
                         corners_open.end());
      if (debug)
        cout << " -- of which "<<corners_open.size()<<" have not already been checked" << endl;

      sat = 1;        // will be set to 0 if we find out that the alphas are not already saturated
      CuspList extra_alphas; // will be filled with any extra alphas needed on this pass
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
              corners_ok.push_back(P);
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
          for ( auto& a : extras)
            {
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
          // Now delete this P from corners and corners_open (in the
          // next loop we'll run find_corners but this adds to the
          // existing list
          if (verbose)
            cout<<"Deleting "<<P<<" from corners list"<<endl;
          corners.erase(std::find(corners.begin(), corners.end(), P));
        } // checked all P in corners_open
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
      for (const auto& a : extra_alphas)
        add_one_alpha(a, 1, verbose>1);
      if (verbose)
        cout<<"alist now has size "<<alist.size()<<endl;
    } // ends while(!sat)

  if (already_saturated)
    {
      if (verbose)
        cout<<"Already saturated"<<endl;
      std::sort(alist.begin(), alist.end(), Cusp_cmp);
      return;
    }

  if (verbose)
    cout << "After saturation we now have "<< alist.size()<<" alphas with max norm "<<m<<endl;

  // Now again delete any alphas with <3 vertices, allowing for translates
  // cornersx.clear();
  // for ( const auto& P : corners)
  //   for ( const auto& t : shifts)
  //     cornersx.push_back(translate(P,t));
  // for ( const auto& s : slist)
  //   {
  //     H3point P = {s, ZERO};
  //     for ( const auto& t : shifts)
  //       cornersx.push_back(translate(P,t));
  //   }
  // alist = remove_redundants(alist, cornersx);
  // m = max_dnorm(alist);
  // if (verbose)
  //   cout << "After removing alphas which go through <3 vertices, we now have "
  //        <<alist.size()<<" alphas with max norm "<< m <<endl;

  std::sort(alist.begin(), alist.end(), Cusp_cmp);
}
