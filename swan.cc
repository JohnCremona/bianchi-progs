// FILE SWAN.CC: implementation of Swan's algorithm

#include "swan.h"
#include "swan_sigmas.h"
#include "swan_alphas.h"
#include "swan_tess.h"
#include "swan_hom.h" // for homology_invariants() only
#include "pari_snf.h" // for homology_invariants_via_pari() only

// clears all alphas-sigma data, polyhedra, faces
void SwanData::clear()
{
  alist.clear(); alistx.clear(); alistF4.clear(); alist_ok.clear(); alist_open.clear(); maxn=0;
  nbrs.clear(); nbrs_ok.clear(); nbrs_open.clear();
  a_denoms.clear(); a_ind.clear(); a_inv.clear(); a_flip.clear(); Mlist.clear();
  slist.clear(); slistx.clear(); s_ind.clear(); s_flip.clear();
  edge_pairs_plus.clear();  edge_pairs_minus.clear();  edge_fours.clear();
  alpha_sets.clear();
  corners.clear();
  singular_polyhedra.clear(); principal_polyhedra.clear(); all_polyhedra.clear();
  aaa.clear();  aas.clear(); sqs.clear(); hexs.clear();
  T_faces.clear(); U_faces.clear(); Q_faces.clear(); H_faces.clear();
  all_faces.clear();
  M32.clear();
}

// clears all and recomputes alphas-sigma data, polyhedra, faces
void SwanData::create(int verbose)
{
  clear();
  make_sigmas();
  make_alphas(verbose);        // calls make_sigmas()
  make_all_polyhedra(verbose); // calls make_alphas()
  make_all_faces(verbose);     // calls make_all_polyhedra()
  encode_all_faces(1, verbose);  // (check=1) calls make_all_faces()
}

// read from geodata file, or create from scratch and store if not successful (file absent)
void SwanData::read_or_create(string subdir, int verbose)
{ if (verbose)
    cout<<"In SwanData::read_or_create()" << endl;
  if (!read(subdir, verbose))
    { if (verbose)
        cout << "Creating SwanData from scratch"<<endl;
      create(verbose);
      if (verbose)
        cout << "Storing SwanData in geodata file"<<endl;
      output_geodata(subdir);
      read(subdir, verbose);
    }
  else
    {
      if (verbose)
        cout << "SwanData read data OK" << endl;
    }
}

void SwanData::make_sigmas() {
  if (slist.empty())
    {
      string step = "SwanData::make_sigmas()";
      SwanTimer.start(step);
      slist = sort_singular_points(singular_points());
      s_flip.reserve(slist.size());
      for (int i=0; i< (int)slist.size(); i++)
        s_flip.push_back(cusp_index_upto_translation(-slist[i], slist));
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

void SwanData::make_alphas(int verbose)
{
  if (alist.empty())
    {
      string step = "SwanData::make_alphas()";
      SwanTimer.start(step);
      make_sigmas();
      find_covering_alphas(verbose);
      saturate_alphas(verbose);
      make_alpha_orbits();
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
  if (a.in_quarter_rectangle() && !is_alpha_redundant(a))
    {
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
      alist_open.insert(a);
      alistF4.push_back(a);
      if (verbose)
        cout<<"appending a="<<a<<" to alistF4"<<endl;
    }
  else
    {
      alist_ok.insert(a);
    }
  // Update nbr lists of earlier alphas (in F4) if they intersect
  // one of the new translates:
  auto a_nbrs = nbrs.at(a); // make a copy as the loop may change the nbrs map
  for (const auto& b : a_nbrs ) // NB these may have been translated
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

int SwanData::add_new_alphas(Quadlooper& denom_looper, int verbose)
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
      if (verbose>1)  cout << new_alphas <<endl;
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

// Return i such that alist[i]=a, else -1
int SwanData::alpha_index(const RatQuad& a) const
{
  auto s = a_ind.find(a.coords());
  if (s==a_ind.end())
    return -1;
  return s->second;
}

// Return i and set t such that alist[i]+t=a, else -1
int SwanData::alpha_index_with_translation(const RatQuad& a, Quad& t) const
{
  auto s = a_ind.find(a.coords());
  if (s==a_ind.end())
    s = a_ind.find(reduce_to_rectangle(a,t).coords());
  if (s==a_ind.end())
    return -1;
  int i = s->second;
  (a-alist[i]).is_integral(t);
  assert (t+alist[i]==a);
  return i;
}

// Return i such that alist[i]+t=a, else -1
// (same as above for when t is not needed)
int SwanData::alpha_index_upto_translation(const RatQuad& a) const
{
  auto s = a_ind.find(a.coords());
  if (s==a_ind.end())
    s = a_ind.find(reduce_to_rectangle(a).coords());
  if (s==a_ind.end())
    return -1;
  return s->second;
}

// Return i such that slist[i]=s, else -1
int SwanData::sigma_index(const RatQuad& s) const
{
  auto i = s_ind.find(s.coords());
  if (i==s_ind.end())
    return -1;
  return i->second;
}

// Return i and set t such that slist[i]+shift=a/b, else -1
int SwanData::sigma_index_with_translation(const Quad& a, const Quad& b, Quad& shift) const
{
  shift = 0;
  if (b.is_zero()) {return 0;}
  int t = 0;
  for ( const auto& sigma : slist )
    {
      Quad r=sigma.num(), s=sigma.den(); // sigma = r/s
      if (s.is_zero())
        {
          t++;
          continue;
        }
      shift = mms(a,s,b,r); //a*s-b*r
      if (div(b*s,shift,shift))  // success! NB shift is divided by b*s before returning
        return t;
      t++;
    } // end of loop over singular points sigma
  return -1;
}

// Return i such that slist[i]+shift=s, else -1
// (same as above for when shift is not needed)
int SwanData::sigma_index_upto_translation(const RatQuad& s) const
{
  auto is = s_ind.find(s.coords());
  if (is==s_ind.end())
    is = s_ind.find(reduce_to_rectangle(s).coords());
  if (is==s_ind.end())
    return -1;
  return is->second;
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
  Quadlooper denom_looper; // default init: norms from 1 to oo, both conjugates
  int ok=0;
  while (!ok)
    {
      // Get the next batch new_alphas, either of all dnorms up to maxn, or those with the next dnorm
      int nc = add_new_alphas(denom_looper, verbose);
      if (nc==0)
        continue;

      // Test whether all alist_open are now surrounded (by alistx).
      // As a side effect, some alphas will be moved from alist_open
      // to alist_ok, and the nbrs_open and nbrs_ok lists are updated.

      ok = are_alphas_surrounded(verbose);
      if (!ok &&verbose)
        {
          cout << "Some alphas are not surrounded, continuing...\n";
          if (verbose>1)
            {
              cout << "alist_ok: " <<alist_ok<<endl;
              cout << "alist_open: " <<alist_open<<endl;
            }
        }
    }

  alist.clear();
  std::copy(alist_ok.begin(), alist_ok.end(), std::back_inserter(alist));
  if (verbose)
    cout << "Success in covering using "<<alist.size()<<" alphas of with max norm "<<maxn<<endl;
  if (verbose>1)
    cout <<alist<<endl;
  SwanTimer.stop(step);
  if (showtimes) SwanTimer.show(1, step);
}

// test if a is singular by reducing to rectangle and comparing with
// slist (but there are two special cases)
int SwanData::is_singular(const RatQuad& a)
{
  RatQuad a0 = reduce_to_rectangle(a);
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

int SwanData::is_alpha_redundant(const RatQuad& a, int verbose)
{
  int debug = verbose>1;
  if (debug)
    cout<<"Testing if "<<a<<" is redundant, with S_a covered by 2 neighbouring S_b..."<<flush;
  const CuspList& a_nbrs = nbrs[a];
  for (auto it = a_nbrs.begin(); it!=a_nbrs.end(); ++it)
    {
      RatQuad b = *it;
      for (auto jt = it+1; jt!=a_nbrs.end(); ++jt)
        {
          RatQuad c = *jt;
          if ( (are_intersection_points_covered_by_one(a, b, c)==2)
               &&
               (are_intersection_points_covered_by_one(a, c, b)==2)
               )
            {
              if (debug)
                cout << "Yes, " << a << " is redundant as covered by S_"<<b<<" and S_"<<c<<endl;
              return 1;
            }
        }
    }
  if (debug)
    cout << "No,  " << a << " is not proved to be redundant"<<endl;
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
  std::set<RatQuad> new_alist_open = alist_open;
  int alphas_redundant = 0;
  int alphas_surrounded = 0;
  int alphas_not_surrounded = 0;
  for ( const auto& a : new_alist_open)
    {
      i++;
      if (debug) cout <<"Testing alpha #"<<i<<"/"<<n_open<<" = "<<a<<"...";
      if (is_alpha_redundant(a, verbose))
        {
          alphas_redundant++;
          if (debug) cout << " ok! redundant" <<endl;
          // add this alpha to the ok list end remove from the open list and full list
          alist_ok.insert(a);
          alist_open.erase(a);
          alist.erase(std::find(alist.begin(), alist.end(), a));
        }
      else
        {
          if (is_alpha_surrounded(a, verbose))
            {
              alphas_surrounded++;
              if (debug) cout << " ok! surrounded" << endl;
              // add this alpha to the ok list end remove from the open list
              alist_ok.insert(a);
              alist_open.erase(a);
            }
          else
            {
              if (debug) cout << " no, not surrounded" << endl;
              alphas_not_surrounded++;
              //return 0;
            }
        }
    }

  if (alphas_not_surrounded)
    {
      if (verbose)
        {
          cout<<alphas_redundant<<" redundant alphas"<<endl;
          cout<<alphas_surrounded<<" surrounded alphas"<<endl;
          cout<<alphas_not_surrounded<<" alphas not yet surrounded"<<endl;
        }
      return 0;
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
  CuspList ans = neighbours(s, alist); // includes translates of a in alist
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
            z2 = reduce_to_rectangle(z2);
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

// Populate the corners list, using preexisting alist and alistF4

void SwanData::find_corners(int debug)
{
  string step = "SwanData::find_corners()";
  SwanTimer.start(step);

  if (debug)
    cout << "Finding corners (triple intersections) for " <<alist.size()<<" alphas..."<<endl;

  // Extend alistF4 by 8 translations:
  CuspList alistF4X;
  for ( auto& z : alistF4)
    {
      CuspList z_nbrs = F4nbrs(z);
      alistF4X.insert(alistF4X.end(), z_nbrs.begin(), z_nbrs.end());
    }
  int n=alistF4X.size();
  if (debug)
    {
      cout << n <<" in first quadrant + neighbours of these" <<endl;
      if (debug>1) cout << alistF4X << endl;
    }

  // For each i get a list of j>i for which S_ai and S_aj intersect properly
  vector<std::set<int> > i2j;
  int i=-1;
  for ( const auto& a : alistF4X)
    {
      i++;
      std::set<int> i2j_i;
      for (int j = i+1; j<n; j++)
        if (circles_intersect(a, alistF4X[j]))
          i2j_i.insert(i2j_i.end(), j);
      i2j.push_back(i2j_i);
    }
  if (debug)
    cout << " finished making i2j" <<endl;

  // Hence make a list of triples (i,j,k) with i<j<k with pairwise
  // proper intersections; test if they make a triple intersection
  // (corner) P=[z,t2] with t2>0 and z in the quarter-rectangle:

  std::set<H3point> cornersF4;
  i=-1;
  for (const auto& ij_list : i2j)
    {
      i++;
      for (const auto& j : ij_list) // so (i,j) intersect
        {
          std::set<int> jk_list = i2j[j], ik_list = ij_list, k_list;
          // so (i,k) and (j,k) also intersect
          std::set_intersection(ik_list.begin(), ik_list.end(),
                                jk_list.begin(), jk_list.end(),
                                std::inserter(k_list, k_list.begin()));
          for (const auto& k : k_list) // so all pairs intersect
            { // Now i<j<k and they intersect pairwise, not necessarily triply
              // cout << "\n(i,j,k) = ("<<i<<","<<j<<","<<k<<")" <<endl;
                H3pointList points1 = tri_inter_points(alistF4X[i], alistF4X[j], alistF4X[k]);
                if (points1.empty()) // it has size 0 or 1
                  continue;
                H3point P = points1.front();
                if (P.t2.sign()==0)
                  continue;
                // cout << " testing P = "<<P<<"..."<<flush;
                if (!P.z.in_quarter_rectangle())
                  continue;
                // cout << " z(P) in F4..."<<flush;
                if (is_under_any(P, alistF4X))
                  continue;
                // cout << " P not under any hemisphere..."<<flush;
                if (debug>2)
                  cout << " found P = "<<P<<"\n";
                auto res = cornersF4.insert(P);
                if (debug && res.second)
                  cout << " adding P = "<<P<<" from (i,j,k) = ("<<i<<","<<j<<","<<k<<")" <<endl;
              }
        }
    }

  if (debug)
    cout << " found "<<cornersF4.size() <<" corners in first quadrant" << endl;
  if (debug>1)
    cout<<cornersF4<<endl;

  // These corners are in F4, so we apply symmetries to get all those in F:
  corners.reserve(4*cornersF4.size());
  static const RAT half(1,2);
  for (const auto& P : cornersF4)
    {
      corners.push_back(P);
      RatQuad z = P.z;
      RatQuad zbar = z.conj();
      vector<RAT> xy = z.coords(1); // so z = x+y*sqrt(-d)
      RAT x=xy[0], y=xy[1];
      if (Quad::t) y *= 2;
      int x_on_edge = (x==half);
      int y_on_edge = (y==half);
      int xy0 = (x==0)||(y==0);

      if (!(x_on_edge || y_on_edge)) // then -P is in the rectangle too
        corners.push_back({-z, P.t2});
      if (!(y_on_edge || xy0)) // then Pbar is in the rectangle too and !=P, !=-P
        corners.push_back({zbar, P.t2});
      if (!(x_on_edge || xy0)) // then -Pbar is in the rectangle too and !=P, !=-P
        corners.push_back({-zbar, P.t2});

    }

  if (debug)
    {
      cout << " found "<<corners.size() <<" corners" <<endl;
      if (debug>2) cout << " Corners (before sorting):\n" << corners << endl;
    }
  std::sort(corners.begin(), corners.end());
  if (debug>2)
    cout << " Corners (after  sorting):\n" << corners << endl;

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
      if (debug)
        cout<<"...number of alphas reduced from "<<alist.size()<<" to "<<alist0.size()<<endl;
      // Test that the singular points are still surrounded:
      int ok = are_sigmas_surrounded();
      if (!ok)
        cout << " but singular points are no longer surrounded" <<endl;

      alist = alist0;
    }

  SwanTimer.stop(step);
  if (showtimes) SwanTimer.show(1, step);
}

#if(0)
// Find potential corners, store in class's corners list, replacing
// alistF4 with sublist of alphas in F4 on >=3 corners
void SwanData::old_find_corners(int verbose)
{
  string step = "SwanData::old_find_corners()";
  SwanTimer.start(step);
  int debug = verbose>1;
  if (verbose)
    {
      cout << "In SwanData.old_find_corners() with " <<alist.size()<<" alphas..."<<endl;
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
      cout << " SwanData.old_find_corners() found "<<corners.size() <<" corners: " << corners<<endl;
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
          b = reduce_to_rectangle(b);
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
#endif // old_find_corners()

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
                  b = reduce_to_rectangle(b);
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
              H3point Q = {reduce_to_rectangle(z2), t2};
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

  // (re)sort alist and corners before returning:
  if (verbose)
    {
      cout << " resorting " << corners.size() << " corners "<<endl;
      if (verbose>1)
        cout << " Corners (before sorting):\n" << corners << endl;
    }
  std::sort(corners.begin(), corners.end());
  if (verbose>1)
    cout << " Corners (after  sorting):\n" << corners << endl;
  std::sort(alist.begin(), alist.end());
  SwanTimer.stop(step);
  if (showtimes) SwanTimer.show(1, step);
}

void SwanData::output_sigmas(int include_small_denoms, string subdir)
{
  make_sigmas();
  vector<Quad> small_denoms;
  if (!include_small_denoms) small_denoms  = {Quad(0), Quad(1), Quad(2), Quad(3)};
  ofstream geodata;
  stringstream ss;
  if (!subdir.empty()) ss << subdir << "/";
  ss << "geodata_" << Quad::d << ".dat";
  geodata.open(ss.str().c_str(), ios_base::app);

  int nlines=0;
  CuspList sigmas_output; // record when s is output so -s is not also output
  for ( const auto& s : slist)
    {
      if (std::find(small_denoms.begin(), small_denoms.end(), s.den()) != small_denoms.end())
        continue;
      // do not output s if s or -s already output
      if (std::find(sigmas_output.begin(), sigmas_output.end(), s) != sigmas_output.end())
        continue;
      nlines++;
      sigmas_output.push_back(s);
      sigmas_output.push_back(-s);
      geodata << make_S_line(s) << endl;
    }
  geodata.close();
}

void SwanData::make_alpha_orbits()
{
  string step = "SwanData::make_alpha_orbits()";
  SwanTimer.start(step);
  CuspList new_alist = alpha_orbits(alist, alpha_sets); // sorts alist as well as finding the orbits
  // cout<<"alpha_sets = "<<alpha_sets<<endl;
  // cout<<"sorted alphas: "<<new_alist<<endl;
  alist.clear();
  for (const auto& sr1r2 : alpha_sets)
    process_alpha_orbit(sr1r2);
  assert (alist==new_alist);
  if (alist!=new_alist)
    {
      cout<<"after processing alpha orbits, alist = "<<alist<<endl;
      exit(1);
    }
  SwanTimer.stop(step);
  if (showtimes) SwanTimer.show(1, step);
}

void SwanData::output_alphas(int include_small_denoms, string subdir)
{
  make_alphas();
  vector<Quad> small_denoms;
  if (!include_small_denoms) small_denoms  = {Quad(0), Quad(1), Quad(2), Quad(3)};

  ofstream geodata;
  stringstream ss;
  if (!subdir.empty()) ss << subdir << "/";
  ss << "geodata_" << Quad::d << ".dat";
  geodata.open(ss.str().c_str()); //  , ios_base::app);

  for ( const auto& sr1r2 : alpha_sets)
    {
      if (std::find(small_denoms.begin(), small_denoms.end(), sr1r2[0]) == small_denoms.end())
        geodata << make_A_line(sr1r2) <<endl;
    }
  geodata.close();
}

void SwanData::process_sigma_orbit(const Quad& r, const Quad& s)
{
  int ns = slist.size();
  RatQuad sigma(r,s), msigma(-r,s);
  slist.push_back(sigma);
  if (!s.is_zero()) // so oo is ignored
    {
      s_ind[sigma.coords()] = ns;
      s_ind[reduce_to_rectangle(sigma).coords()] = ns;
    }
  static const Quad two(2);
  if (s.is_zero() || s==two) // don't also include -sigma
    {
      s_flip.push_back(ns);     // identity
    }
  else
    {
      slist.push_back(msigma);
      s_ind[msigma.coords()] = ns+1;
      s_ind[reduce_to_rectangle(msigma).coords()] = ns+1;
      s_flip.push_back(ns+1);   // transposition with next
      s_flip.push_back(ns);     // transposition with previous
    }
  // cout<<"After process_sigma_orbit("<<r<<","<<s<<"):"<<endl;
  // cout<<"slist = "<<slist<<endl;
  // cout<<"s_flip = "<<s_flip<<endl;
}

void SwanData::process_alpha_orbit(const Quad& s, const Quad& r1, const Quad& r2, int verbose)
{
  if (verbose)
    cout<<"Processing alpha orbit ("<<s<<","<<r1<<","<<r2<<")\n";

  auto add_alpha = [this](const Quad& a, const Quad& b, const Quad& c, const Quad& d)
  {
    RatQuad alpha(-d,c);
    int n = alist.size(); // = index of this alpha once we have added it
    a_ind[alpha.coords()] = n;
    alist.push_back(alpha);
    mat22 M(a,b,c,d);  // maps alpha = -d/c to oo
    assert (M.is_unimodular());
    assert (M(alpha).is_infinity());
    Mlist.push_back(M);
    a_denoms.insert(c);
    a_ind[alpha.coords()] = n;
    alpha = reduce_to_rectangle(alpha);
    a_ind[alpha.coords()] = n;
  };

  Quad t = -(r1*r2+Quad::one), two(2);
  assert(div(s,t));
  t /= s;
  int s_divides_2 = div(s,two);
  int n = alist.size(); // = index of this alpha once we have added it

  if (r1==r2) // "-" pair, r1^2=-1 (mod s)
    {
      // General case alpha[n] = r/s, alpha[n+1] = -r/s are each self-inverse
      // Special case alpha[n] = r/s = -r/s (if s|2) is self-inverse
      edge_pairs_minus.push_back(n);
      if (s_divides_2) // inv and flip are identity (special case: only one alpha)
        {
          a_inv.push_back(n);
          a_flip.push_back(n);
        }
      else // inv is identity, flip swaps n,n+1 (general case: two alphas)
        {
          // n,n+1 map to n,n+1
          a_inv.push_back(n);
          a_inv.push_back(n+1);
          // n,n+1 map to n+1,n
          a_flip.push_back(n+1);
          a_flip.push_back(n);
        }
      add_alpha( r1, t, s, -r1); // alpha =  r1/s
      if (!s_divides_2)
        add_alpha(-r1, t, s,  r1); // alpha = -r1/s
      if (verbose>1)
        cout<<"(-pair): alist is now "<<alist<<", #alphas="<<alist.size()
            <<", edge_pairs_minus="<<edge_pairs_minus<<endl;
      return;
    }

  if (r1==-r2) // "+" pair, r1^2=+1 (mod s) -- here we do not have s|2
    {
      // General case alpha[n] = r/s, alpha[n+1] = -r/s are inverses of each other
      // Special case alpha[n] = r/s = -r/s (if s|2) does not occur
      edge_pairs_plus.push_back(n);
      // n,n+1 map to n+1,n
      a_inv.push_back(n+1);
      a_inv.push_back(n);
      // n,n+1 map to n+1,n
      a_flip.push_back(n+1);
      a_flip.push_back(n);
      add_alpha(-r1, t, s, -r1); // alpha =  r1/s
      add_alpha( r1, t, s,  r1); // alpha = -r1/s
      if (verbose>1)
        cout<<"(+pair): alist is now "<<alist<<", #alphas="<<alist.size()<<", edge_pairs_plus="<<edge_pairs_plus<<endl;
      return;
    }

  // Now we have in general four distinct alphas r1/s, -r1/s, r2/s, -r2/s with r1*r2=-1 (s).
  // Special case (s|2): just two alphas, r1/s, r2/s with r1*r2=-1=+1 (mod s) and r1!=r2.
  edge_fours.push_back(n);
  if (s_divides_2) // inv swaps n,n+1; flip is identity
    {
      // n,n+1 map to n+1,n
      a_inv.push_back(n+1);
      a_inv.push_back(n);
      // n,n+1 map to n,n+1
      a_flip.push_back(n);
      a_flip.push_back(n+1);
    }
  else // inv swaps n,n+2 and n+1,n+3; flip swaps n,n+1 and n+2,n+3
    {
      // n,n+1,n+2,n+3 map to n+2,n+3,n,n+1
      a_inv.push_back(n+2);
      a_inv.push_back(n+3);
      a_inv.push_back(n);
      a_inv.push_back(n+1);
      // n,n+1,n+2,n+3 map to n+1,n,n+3,n+2
      a_flip.push_back(n+1);
      a_flip.push_back(n);
      a_flip.push_back(n+3);
      a_flip.push_back(n+2);
    }
  // cout<<"(four): before adding alphas, #alphas="<<n<<endl;
  add_alpha( r2, t, s, -r1); // alpha =  r1/s
  if (!s_divides_2)
    add_alpha(-r2, t, s,  r1); // alpha = -r1/s
  add_alpha( r1, t, s, -r2); // alpha =  r2/s
  if (!s_divides_2)
    add_alpha(-r1, t, s,  r2); // alpha = -r2/s
  if (verbose>1)
    cout<<"(four): alist is now "<<alist<<", #alphas="<<alist.size()<<", edge_fours="<<edge_fours<<endl;
}

int SwanData::read_alphas_and_sigmas(int include_small_denoms, string subdir, int verbose)
{
  alist.clear(); a_denoms.clear(); a_ind.clear(); a_inv.clear(); a_flip.clear(); Mlist.clear();
  slist.clear(); s_ind.clear(); s_flip.clear();
  edge_pairs_plus.clear();  edge_pairs_minus.clear();  edge_fours.clear();

  Quad w = Quad::w, zero(0), one(1), two(2), three(3);

  // Start alist with 0/1 amd slist with sigma = oo:

  process_alpha_orbit(one, zero, zero, verbose);
  slist = {RatQuad::infinity()};
  s_flip = {0};

  if (Quad::is_Euclidean) return 1;

  ifstream geodata;
  stringstream ss;
  if (!subdir.empty()) ss << subdir << "/";
  ss << "geodata_" << Quad::d << ".dat";
  geodata.open(ss.str().c_str());
  if (!geodata.is_open())
    {
      if (verbose) cout << "No geodata file!" <<endl;
      clear();
      return 0;
    }
  else
    {
      if (verbose)
        cout << "Reading A- and S- lines from geodata file" <<endl;
    }

  // Add in alphas and sigmas with small denoms "manually" if they are not on file:

  if (!include_small_denoms) // means they are *not* in the geodata file
    {
      if (verbose)
        cout << "Constructing alphas and sigmas of denominator 1,2,3..."<<endl;
      int d = Quad::d;
      int d8 = d%8, d12 = d%12;

      // sigmas of denom 2 and 3

      if (Quad::class_number > 1)
        {
          // sigmas of denom 2 (if any)

          switch (d8) {
          case 1: case 5:
            process_sigma_orbit(1+w,two);
            break;
          case 2: case 6:
            process_sigma_orbit(w,two);
            break;
          case 7:
            process_sigma_orbit(w,two);
            process_sigma_orbit(1-w,two); // NB not in rectangle: (w-1)/2 is
            break;
          default:
            ;
          } // d8 switch

          // sigmas of denom 3 (if any)

          if (!(d==5 || d==6 || d==23 || d==15))
            {
              switch (d12) {
              case 2: case 5:
                process_sigma_orbit(1+w,three);
                process_sigma_orbit(1-w,three);
                break;
              case 3:
                process_sigma_orbit(1+w,three);
                // NB (-1-w)/3 is not in rectangle: (2-w)/3 is
                break;
              case 6: case 9:
                process_sigma_orbit(w,three);
                break;
              case 11:
                process_sigma_orbit(w,three);
                process_sigma_orbit(-1+w,three);
                break;
              default:
                ;
              } // d12 switch
            } // excluded d=5,6,15,23
        }  // class number>1

      // alphas of denom 2 (if any)

      switch (d8) {
      case 1: case 5:
        process_alpha_orbit(two, w, w, verbose);
        break;
      case 2: case 6:
        process_alpha_orbit(two, 1+w, 1+w, verbose);
        break;
      case 3:
        process_alpha_orbit(two, w, w-1, verbose);
        break;
      default:
        ;
      } // d8 switch

      // alphas of denom 3 (if any)

      if (!( (d==5 || d==6 || d==15 || d==19 || d==23) ))
        {
          switch (d%12) {
          case 1: case 10:
            process_alpha_orbit(three, w, w, verbose);       // w^2=-1 (mod 3)
            process_alpha_orbit(three, 1+w, 1-w, verbose);   // (1+w)(1-w)=-1 (mod 3)
            break;
          case 7:
            if (d>31)
              process_alpha_orbit(three, w, 1-w, verbose);    // w(1-w)=-1 (mod 3)
            process_alpha_orbit(three, 1+w, 1+w, verbose);    // (1+w)^2=-1 (mod 3)
            break;
          case 2: case 5:
            process_alpha_orbit(three, w, -w, verbose);         // w^2=+1 (mod 3)
            break;
          case 11:
            process_alpha_orbit(three, 1+w, -1-w, verbose);     // (1+w)^2=+1 (mod 3)
            break;
          case 3:
            process_alpha_orbit(three, w, w-1, verbose);     // w(w-1)=-1 (mod 3)
            break;
          case 6: case 9:
            process_alpha_orbit(three, w+1, w-1, verbose);   // (w+1)(w-1)=-1 (mod 3)
            break;
          } // d12 switch
        } // small fields
    } // small denoms

  // Now read A- and S- lines from geodata file for alphas and sigmas
  // of other denominators:

  string line;
  int nA=0, nS=0;
  getline(geodata, line);
  while (!geodata.eof())
    {
      int file_d;
      char G;
      POLYGON poly =  parse_geodata_line(line, file_d, G, verbose>1);
      auto data = poly.shifts;
      switch(G) {
      case 'A': // alpha orbit
        {
          process_alpha_orbit(data[0], data[1], data[2], verbose>1);
          nA++;
          break;
        }
      case 'S': // sigma orbit
        {
          process_sigma_orbit(data[0], data[1]);
          nS++;
          break;
        }
      default:
        {
          break;
        }
      }
      if (G=='X')
        break; // don't read any more lines
      else
        getline(geodata, line);
    }
  geodata.close();
  if (verbose)
    {
      cout << alist.size() << " alphas in " << nA << " orbits" <<endl;
      cout << slist.size() << " sigmas in " << nS << " orbits" << endl;
      if (verbose>1)
        {
          cout << "alphas: " << alist << endl;
          cout << "sigmas: " << slist << endl;
          cout << "alpha_inv: " << a_inv << endl;
          cout << "edge_pairs_plus: " << edge_pairs_plus << endl;
          cout << "edge_pairs_minus: " << edge_pairs_minus << endl;
          cout << "edge_fours: " << edge_fours << endl;
        }
    }
  return 1;
}

void SwanData::make_singular_polyhedra(int verbose)
{
  string step = "SwanData::make_singular_polyhedra()";
  SwanTimer.start(step);
  singular_polyhedra = ::singular_polyhedra(slist, alist, verbose);
  SwanTimer.stop(step);
  if (showtimes) SwanTimer.show(1, step);
}

POLYHEDRON SwanData::make_principal_polyhedron(const H3point& P, std::set<int>& orbit, int verbose)
{
  if (verbose)
    cout << " - constructing principal polyhedron around corner "<<P<<"...\n";
  POLYHEDRON poly;
  RatQuad infty = RatQuad::infinity();

  // Find all a with P on S_a; these & oo are the vertices of the polyhedron
  CuspList alistP = covering_hemispheres(P);
  poly.vertices = alistP;
  poly.vertices.insert(poly.vertices.begin(), infty);

  int nv = poly.vertices.size();
  if (verbose)
    cout << " - polyhedron has "<< nv << " vertices ("<<alistP.size()
         <<" S_a go through P: "<<alistP<<")"<<endl;

  // local function to test for being fundamental or oo:
  Quad x;
  auto is_fund = [this, &x](const RatQuad& a) {
    return a.is_infinity() || cusp_index_with_translation(a, alist, x)>=0;
  };

  // check off flags j and k where P=corners[j], u*P=corners[k] (mod translation)
  int i, j, k;
  j = point_index_with_translation(P, corners, x);
  orbit.insert(j);
  if (verbose>1)
    cout<<"Checking off corner #"<<j<<" ("<<corners[j]<<")\n";
  H3point Q = negate(P);
  if (verbose>1)
    cout<<" Looking for corner "<<Q<<endl;
  k = point_index_with_translation(Q, corners, x);
  orbit.insert(k);
  if (verbose>1)
    cout<<" Q = corner #"<<k<<" ("<<corners[k]<<")"<<endl;

  for ( const auto& a : alistP)
    {
      // First add edges {oo,a} and {a,oo} for a fundamental:
      if (is_fund(a))
        {
          poly.edges.push_back({infty,a});
          poly.edges.push_back({a,infty});
        }
      // Then add edges {a,b} for finite a,b, when M_a(b) is fundamental.
      // NB this is equivalent to M_b(a) fundamental, so {b,a} will also be added.
      mat22 M = Malpha(a, P, corners, i); // M(a)=oo, M(P)=corners[i]
      for ( const auto& b : alistP)
        {
          if ((a!=b) && is_fund(M(b)))
            poly.edges.push_back({a,b});
        }
      // add flag i to orbit
      orbit.insert(i);
      if (verbose>1)
        cout<<"Checking off corner #"<<i<<" ("<<corners[i]<<")\n";

      // check off flag k where u*corners[i]=corners[k] (mod translation)
      Q = negate(P);
      if (verbose>1)
        cout<<" Looking for corner "<<Q<<endl;
      k = point_index_with_translation(Q, corners, x);
      orbit.insert(k);
      if (verbose>1)
        cout<<" Q = corner #"<<k<<" ("<<corners[k]<<")"<<endl;
    }
  if (verbose)
    {
      cout << " - orbit " << orbit << " of size " << orbit.size() << endl;
      if (verbose>1)
        {
          int ne = poly.edges.size()/2;
          int nf = 2+ne-nv; // Euler's formula!
          cout << " - polyhedron has (V,E,F)=("<<nv<<","<<ne<<","<<nf<<"):\n"; //<<poly << endl;
          cout << " - now filling in face data..."<<endl;
        }
    }
  fill_faces(poly, verbose);
  if (verbose)
    {
      cout << "After filling in faces, polyhedron is a " << poly_name(poly) << " with "<<poly.faces.size()<<" faces:\n"<<poly.faces << endl;
      // cout << "VEF(poly) = " << VEF(poly) <<endl;
      cout << "VEFx(poly) = " << VEFx(poly) <<endl;
    }
  return poly;
}

void SwanData::make_principal_polyhedra(int verbose)
{
  string step = "SwanData::make_principal_polyhedra()";
  SwanTimer.start(step);
  if (verbose)
    cout<<"Finding principal polyhedra from "<<alist.size()<<" alphas with "
        << corners.size()<<" corners..."<<endl; //<<"Corners:\n"<<corners<<endl;

  vector<int>flags(corners.size(), 0);
  int j=-1;
  for ( const auto& P : corners)
    {
      j++;
      if (flags[j])
        continue;
      if (verbose)
        cout<<"Using corner #"<<j<<"..."<<flush;
      std::set<int> orbit;
      principal_polyhedra.push_back(make_principal_polyhedron(P, orbit, verbose));
      if (verbose)
        cout<<"Corner orbit "<<orbit<<"..."<<endl;
      for ( int i : orbit) flags[i] = 1;
    }
  if (verbose)
    cout<<principal_polyhedra.size()<<" principal polyhedra constructed"<<endl;
  SwanTimer.stop(step);
  if (showtimes) SwanTimer.show(1, step);
}

void SwanData::make_all_polyhedra(int verbose)
{
  if (!all_polyhedra.empty()) return;
  make_alphas(verbose); // makes sure all sigma and alpha data is present

  string step = "SwanData::make_all_polyhedra()";
  SwanTimer.start(step);

  make_singular_polyhedra(verbose);
  all_polyhedra = singular_polyhedra;
  make_principal_polyhedra(verbose);
  all_polyhedra.insert(all_polyhedra.end(), principal_polyhedra.begin(), principal_polyhedra.end());

  SwanTimer.stop(step);
  if (showtimes) SwanTimer.show(1, step);
}

void SwanData::make_all_faces(int verbose)
{ verbose=1;
  if (!all_faces.empty()) return;
  make_all_polyhedra(verbose); // makes polyhedra (and all sigma and alpha data if necessary)

  string step = "SwanData::make_all_faces()";
  SwanTimer.start(step);

  vector<int> redundant_faces;
  all_faces = get_faces(all_polyhedra, alist, slist, M32, redundant_faces, verbose);

  // Split up faces into 4 types:
  int i=0;
  for (const auto& face: all_faces)
    {
      int sing = is_face_singular(face, slist);
      int red = std::find(redundant_faces.begin(), redundant_faces.end(), i)!=redundant_faces.end();
      int n = face.size();
      string s;
      switch (n) {
      case 4:
        {
          if (!red) sqs.push_back(face);
          s = "square";
          break;
        }
      case 6:
        {
          if (!red) hexs.push_back(face);
          s = "hexagon";
          break;
        }
      case 3: default:
        {
          if (sing)
            {
              if (!red) aas.push_back(face);
              s = "aas triangle";
            }
          else
            {
              if (!red) aaa.push_back(face);
              s = "aaa triangle";
            }
        }
      }
      if (verbose)
        {
          cout<<i<<" ("<<s;
          if (red) cout << " (redundant)";
          cout<<"): "<<face<<endl;
        }
      i++;
    }

  SwanTimer.stop(step);
  if (showtimes) SwanTimer.show(1, step);
}

// Encode all faces found as POLYGONs in T_faces, U_faces, H_faces,
// Q_faces; report if verbose; check their encodings/decodings for
// consistency if check.

// NB T_faces will not include the universal triangle {0,oo,1} which
// is handled separately in face_relations code (for historic
// reasons).
int SwanData::encode_all_faces(int check, int verbose)
{
  make_all_faces(verbose); // will first make polyhedra if necessary

  int sing, ok, all_ok = 1;
  POLYGON P;

  if (verbose)
    cout<<aaa.size()<<" aaa-triangles\n";
  T_faces.clear();
  for ( const auto& face : aaa)
    {
      if (verbose) cout <<"T "<<face << " --> ";
      P = make_polygon(face, alist, slist, sing);
      if (verbose) cout <<"["<<P.indices<<","<<P.shifts<<"] : ";
      ok = (!check) || check_aaa_triangle(P, verbose);
      if (ok)
        {
          if (is_universal(face, alist, slist))
            {
              if (verbose)
                cout << " (universal triangle)" << endl;
              continue;
            }
          // if (is_standard(face, alist, slist))
          //   {
          //     if (verbose)
          //       cout << " (standard triangle)" << endl;
          //     continue;
          //   }
          T_faces.push_back(P);
        }
      else
        {
          cout<<"aaa-triangle fails encoding check"<<endl;
          all_ok = 0;
        }
    }

  if (verbose)
    cout<<aas.size()<<" aas-triangles\n";
  U_faces.clear();
  for ( const auto& face : aas)
    {
      if (verbose) cout <<"U "<<face << " --> ";
      P = make_polygon(face, alist, slist, sing);
      if (verbose) cout <<"["<<P.indices<<","<<P.shifts<<"] : ";
      ok = (!check) || check_aas_triangle(P, verbose);
      if (ok)
        U_faces.push_back(P);
      else
        {
          cout<<"aas-triangle fails encoding check"<<endl;
          all_ok = 0;
        }
    }

  if (verbose)
    cout<<sqs.size()<<" squares\n";
  Q_faces.clear();
  for ( const auto& face : sqs)
    {
      if (verbose) cout <<"Q "<<face << " --> ";
      P = make_polygon(face, alist, slist, sing);
      if (verbose) cout <<"["<<P.indices<<","<<P.shifts<<"] : ";
      ok = (!check) || check_square(P, verbose);
      if (ok)
        Q_faces.push_back(P);
      else
        {
          cout<<"square fails encoding check"<<endl;
          all_ok = 0;
        }
    }

  if (verbose)
    cout<<hexs.size()<<" hexagons\n";
  H_faces.clear();
  for ( const auto& face : hexs)
    {
      if (verbose) cout <<"H "<<face << " --> ";
      P = make_polygon(face, alist, slist, sing);
      if (verbose) cout <<"["<<P.indices<<","<<P.shifts<<"] : ";
      ok = (!check) || check_hexagon(P, verbose);
      if (ok)
        H_faces.push_back(P);
      else
        {
          cout<<"hexagon fails encoding check"<<endl;
          all_ok = 0;
        }
    }
  return all_ok;
}

// For use after reading face data from a geodata file.  From T_faces,
// U_faces, Q_faces, H_faces reconstruct all_faces.  Must include the
// the universal triangle {0,oo,1} not included in geodata files.
void SwanData::decode_all_faces()
{
  all_faces.clear();

  // Universal triangle
  static const CuspList tri0 = {{0,0,1}, {1,0,0}, {1,0,1}}; // {0,oo,1}
  all_faces.push_back(tri0);

  std::transform(T_faces.begin(), T_faces.end(), std::back_inserter(all_faces),
                 [this] (const POLYGON& P) {return remake_triangle(P, alist, slist, 0);});

  std::transform(U_faces.begin(), U_faces.end(), std::back_inserter(all_faces),
                 [this] (const POLYGON& P) {return remake_triangle(P, alist, slist, 1);});

  std::transform(Q_faces.begin(), Q_faces.end(), std::back_inserter(all_faces),
                 [this] (const POLYGON& P) {return remake_quadrilateral(P, alist);});

  std::transform(H_faces.begin(), H_faces.end(), std::back_inserter(all_faces),
                 [this] (const POLYGON& P) {return remake_hexagon(P, alist);});
}

void SwanData::output_face_data(string subdir, int verbose)
{
  ofstream geodata;
  stringstream ss;
  if (!subdir.empty()) ss << subdir << "/";
  ss << "geodata_" << Quad::d << ".dat";
  geodata.open(ss.str().c_str(), ios_base::app); // append
  int nlines=0;

  for (const POLYGON& T : T_faces)
    {
      string s = polygon_string(T, 0);
      nlines++;
      geodata << s << endl;
      if (verbose) cout << s << endl;
    }
  for (const POLYGON& U : U_faces)
    {
      string s = polygon_string(U, 1);
      nlines++;
      geodata << s << endl;
      if (verbose) cout << s << endl;
    }
  for (const POLYGON& Q : Q_faces)
    {
      string s = polygon_string(Q, 0);
      nlines++;
      geodata << s << endl;
      if (verbose) cout << s << endl;
    }
  for (const POLYGON& H : H_faces)
    {
      string s = polygon_string(H, 0);
      nlines++;
      geodata << s << endl;
      if (verbose) cout << s << endl;
    }

  if (verbose)
    cout << nlines << " lines output" <<endl;

  geodata.close();
}

// Read from subdir/geodata_d.dat all TUQH lines and fill
// aaa, aas, sqs, hexs (not all_faces)
int SwanData::read_face_data(string subdir, int verbose)
{
  T_faces.clear(); U_faces.clear(); Q_faces.clear(); H_faces.clear();

  ifstream geodata;
  stringstream ss;
  if (!subdir.empty()) ss << subdir << "/";
  ss << "geodata_" << Quad::d << ".dat";
  string infile = ss.str();
  geodata.open(infile.c_str());
  if (!geodata.is_open())
    {
      cout << "No geodata file " << infile << " exists!" <<endl;
      return 0;
    }
  if (verbose)
    cout << "reading from geodata file " << infile <<endl;

  string line;
  int nT=0, nU=0, nQ=0, nH=0;
  getline(geodata, line);
  while (!geodata.eof())
    {
      int file_d;
      char G;
      POLYGON poly = parse_geodata_line(line, file_d, G, verbose);
      switch(G) {
      case 'X':
        {
          if (verbose>1)
            cout<<"Skipping rest of file"<<endl;
          break;
        }
      case '%':
      case 'A': // alpha orbit
      case 'S': // sigma orbit
      default:
        {
          if (verbose>1)
            cout<<"Skipping line: "<<line<<endl;
          break;
        }
      case 'T': // aaa-triangle
        {
          T_faces.push_back(poly);
          nT++;
          break;
        }
      case 'U': // aas-triangle
        {
          U_faces.push_back(poly);
          nU++;
          break;
        }
      case 'Q': // square
        {
          Q_faces.push_back(poly);
          nQ++;
          break;
        }
      case 'H': // hexagon
        {
          H_faces.push_back(poly);
          nH++;
          break;
        }
      } // end of switch
      if (G=='X')
        break; // don't read any more lines
      else
        getline(geodata, line);
    }
  geodata.close();
  if (verbose)
    cout<<"Read "<<nT<<" T triangles, "<<nU<<" U-triangles, "<<nQ<<" squares and "<<nH<<" hexagons"<<endl;
  if (check_faces(verbose>1))
    {
      if(verbose)
        cout<<"All polygons on file check OK"<<endl;
      return nT || nU || nQ || nH;
    }
  else
    {
      cout<<"Some polygons on file fail to check OK"<<endl;
      return 0;
    }
}

// Return the index of an edge {a,b} in the range 0..#alphas+#sigmas-2
// For e={a,b} with a,b principal return i (>=0), where [oo,alphas[i]]
// is SL(2,Ok)-equivalent to e.

// If a is principal and b singular, return i+len(alphas)-1 (>=1),
// where [oo,sigmas[i]] is SL(2,Ok)-equivalent to e.

// If a is singular and b principal, return -(i+len(alphas)-1) (<0),
// where [sigmas[i],oo] is SL(2,Ok)-equivalent to e.

// Raise an error if both a and b are singular.

// Note that we assume sigmas[0]=oo which is not a singular point;
// the number of singular points is len(sigmas)-1.

int SwanData::edge_index(const EDGE& e)
{
  RatQuad a = e.alpha(), b = e.beta();
  Quad temp;
  int i;

  auto is_finite_singular = [this, &temp](const RatQuad& c)
  {
    int j = cusp_index_with_translation(c, slist, temp);
    return (j>0? j : 0);
  };

  if (is_finite_singular(b))
    {
      assert (!is_finite_singular(a));
      mat22 M = Malpha(a);  // M(a)=oo
      i = cusp_index_with_translation(M(b), slist, temp);
      assert (i>0);
      return alist.size() + i-1;
    }

  // now b is principal
  if (is_finite_singular(a))
    return - edge_index(EDGE(b,a));

  // now a and b are principal
  mat22 M = Malpha(a);  // M(a)=oo
  i = cusp_index_with_translation(M(b), alist, temp);
  assert (i>=0);
  return i;
}

// Return the edge boundary matrix M10 (matrix of delta: 1-chains -> 0-chains).
// The edge basis consists of first the [a,oo] for a in alist (whose
// boundary is trivial), then the [s,oo] for s in slist[1:].  The
// cusp basis is indexed by ideal classes. Same for SL2 and GL2.
vector<vector<int>> SwanData::edge_boundary_matrix()
{
  int na = n_alph(), ns = n_sig(1); // omit oo
  unsigned int nrows = na+ns;
  unsigned int h = Quad::class_number;
  // Each row has h columns, and there are na+ns rows
  vector<vector<int>> M;
  M.reserve(na+ns);
  for (int i=0; i<na; i++) // first na rows are 0 (boundary of principal edge)
    {
      vector<int> row(h, 0);
      M.push_back(row);
    }
  for (int i=0; i<ns; i++) // next ns rows are boundaries of edge {oo,s}
    {
      vector<int> row(h, 0);
      row[0] = -1;
      row[slist[i+1].ideal_class()] = +1;
      M.push_back(row);
    }
  assert (M.size()==nrows);
  return M;
}

// Return the image under delta of the face, as a vector of length #alist+#slist-1
vector<int> SwanData::face_boundary_vector(const CuspList& face)
{
  vector<int> v(n_alph()+n_sig(1),0);
  int n=face.size();
  for (int i=0; i<n; i++)
    {
      int j = edge_index(EDGE(face[i], face[(i+1)%n]));
      if (j<0) v[-j]-=1; else v[j]+=1;
    }
  return v;
}

// Use alist, slist, edge_pairs_plus, edge_pairs_minus, fours to
// return a matrix with one row per pair of glued oriented edges:
vector<vector<int>> SwanData::edge_pairings(int GL2, int debug)
{
  if(debug) cout<<"In edge_pairings("<<GL2<<")"<<endl;
  int nplus = edge_pairs_plus.size(), nminus = edge_pairs_minus.size(), nfours = edge_fours.size();
  unsigned int ncols = n_alph() + n_sig(1); // size of edge basis
  unsigned int nrows = nplus + nminus + nfours + (GL2? ncols : nfours);
  vector<vector<int>> M;
  M.reserve(nrows);

  int i;
  int j;
  Quad temp, s, r, r1, r2;

  // edge identifications
  // (0) GL2 only: {a,oo}={-a,oo} for a in alist and {s,oo}={-s,oo} for s in slist
  // (1) {a,oo}+{-a,oo}=0     for a=r/s, (r,s) in pluspairs
  // (2) {a,oo}=0            for a=r/s, (r,s) in minuspairs
  // (3) {a1,oo}+{a2,oo}=0    for a1=r1/s, a2=r2/s, (s,r1,r2) in fours
  //   and if not GL2: {-a1,oo}+{-a2,oo}=0 for a1=r1/s, a2=r2/s, (s,r1,r2) in fours

  // NB When d=3 (mod 8), d>11 both r=w/2 and r'=(w-1)/2 are alphas,
  // these form a "four" which only has 2 distinct elements
  long d = Quad::d;
  int special = (d%8==3 && d>11);
  if (debug&&special)
    cout<<"Special case for d="<<d<<" as there's a four of size only 2"<<endl;
  if (GL2) // type (0)
    {
      int nalphas = n_alph();
      for (i=0; i<nalphas; i++)
        {
          vector<int> row(ncols, 0);
          j = cusp_index_with_translation(-alist[i], alist, temp);
          assert ((j>=0)&&(j<nalphas));
          row[i] +=1;
          row[j] -=1;
          if (debug)
            cout<<"i="<<i<<": alphas[i]="<<alist[i]<<"; j="<<j<<": alphas[j]="<<alist[j]<<" so row is "<<row<<endl;
          M.push_back(row);
        }
      int nsigmas = n_sig(1);
      for (i=0; i<nsigmas; i++)
        {
          vector<int> row(ncols, 0);
          j = cusp_index_with_translation(-slist[i+1], slist, temp) -1;
          assert ((j>=0)&&(j<nsigmas));
          row[nalphas+i] +=1;
          row[nalphas+j] -=1;
          M.push_back(row);
        }
    }
  for (i=0; i<nplus; i++)
    {
      vector<int> row(ncols, 0);
      j = edge_pairs_plus[i];
      row[j] +=1;
      row[j+1] +=1;
      if (debug) cout<<"edge_pairs_plus["<<i<<"]="<<j<<" so row is "<<row<<endl;
      M.push_back(row);
    }
  // type (2)
  for (i=0; i<nminus; i++)
    {
      vector<int> row(ncols, 0);
      j = edge_pairs_minus[i];
      // row[j] +=2; // no: if g(e)=-e then e=0, not 2e=0
      row[j] +=1;
      if (debug) cout<<"edge_pairs_minus["<<i<<"]="<<j<<" so row is "<<row<<endl;
      M.push_back(row);
    }
  // types (3) and (3')
  for (i=0; i<nfours; i++)
    {
      vector<int> row(ncols, 0);
      j = edge_fours[i];
      row[j] +=1;
      if (i==0 && special)
        row[j+1] +=1;
      else
        row[j+2] +=1;
      M.push_back(row);
      if (debug) cout<<"edge_fours["<<i<<"]="<<j<<" so row is "<<row<<endl;
      if (!GL2)
        {
          vector<int> row2(ncols, 0);
          if (!(i==0 && special))
            {
              row2[j+1] +=1;
              row2[j+3] +=1;
            }
          M.push_back(row2);
        }
    }
  assert (M.size()==nrows);
  if (debug) cout<<"edge_pairings() returns a "<<nrows<<" x "<<ncols<<" matrix \n"<<M<<endl;
  return M;
}

// Return a matrix with one row per face in all_faces, giving its
// boundary as a Z-linear combination of oriented edges
vector<vector<int>> SwanData::face_boundaries(int GL2, int debug)
{
  if (debug) cout<<"In face_boundaries("<<GL2<<")"<<endl;
  unsigned int nrows = all_faces.size() * (GL2? 1 : 2);
  if (debug) cout<<"nrows = "<<nrows<<endl;
  vector<vector<int>> M;
  M.reserve(nrows);
  for (const auto& face : all_faces)
    {
      M.push_back(face_boundary_vector(face));
      if (!GL2)
        M.push_back(face_boundary_vector(negate_polygon(face)));
    }
  if (debug) cout << "Matrix has size "<<nrows<<" x "<< M[0].size() << ":\n" << M << endl;
  assert (M.size()==nrows);
  return M;
}

// Row concatenation of previous 2, giving the matrix M21 of the
// boundary map from 2-cells to 1-cells
vector<vector<int>> SwanData::face_boundary_matrix(int GL2, int debug)
{
  vector<vector<int>> M = edge_pairings(GL2, debug);
  if (debug)
    {
      cout << "edge pairing matrix has size " << M.size() << " x " << M[0].size() << endl;
      if (debug>1)
        cout << M << endl;
    }
  vector<vector<int>> M2 = face_boundaries(GL2, debug);
  if (debug)
    {
      cout << "face boundaries matrix has size " << M2.size() << " x " << M2[0].size() << endl;
      if (debug>1)
        cout << M2 << endl;
    }
  M.insert(M.end(), M2.begin(), M2.end());
  return M;
}

// return the invariants of H_1 as a Z-module for either GL2
// (group=1) or SL2 (group=2) or both (group=3)
vector<vector<INT>> SwanData::old_integral_homology(int group, int debug)
{
  if (all_faces.empty()) decode_all_faces();

  string step = "SwanData::integral_homology()";
  SwanTimer.start(step);

  vector<vector<int>> M10 = edge_boundary_matrix();
  if (debug)
    {
      cout << "edge boundary matrix M10 has size " << M10.size() << " x " << M10[0].size() << endl;
      if (debug>1) cout << "M10 = \n" << M10 << endl;
    }

  vector<vector<INT>> invs;
  int
    GL2 = group&1, // i.e. 1 or 3
    SL2 = group&2; // i.e. 2 or 3
  if (GL2)
    {
      vector<vector<int>> M21 = face_boundary_matrix(1, debug>1);
      if (debug)
        cout << "GL2 face boundary matrix M21 has size " << M21.size() << " x " << M21[0].size() << endl;
      invs.push_back(homology_invariants(M10, M21, debug));
    }
  if (SL2)
    {
      vector<vector<int>> M21 = face_boundary_matrix(0, debug>1);
      if (debug)
        cout << "SL2 face boundary matrix M21 has size " << M21.size() << " x " << M21[0].size() << endl;
      invs.push_back(homology_invariants(M10, M21, debug));
    }
  SwanTimer.stop(step);
  if (showtimes)
    SwanTimer.show(1, step);
  return invs;
}

// return the invariants of H_1 as a Z-module for either GL2
// (group=1) or SL2 (group=2) or both (group=3)
vector<vector<INT>> SwanData::integral_homology(int group, int debug)
{
  if (all_faces.empty()) decode_all_faces();

  string step = "SwanData::new_integral_homology()";
  SwanTimer.start(step);

  // find 'redundant' edges
  vector<int> class_flags(Quad::class_number, 0);
  int n_edges = n_alph()+n_sig(1);
  vector<int> edge_flags(n_edges, 0);
  int i = n_alph(), j;
  for (auto s: slist)
    {
      if (s.is_finite()) // ignore s = oo in trivial class
        {
          j = s.ideal_class();
          if (!class_flags[j]) // this is a new ideal class
            {
              class_flags[j] = 1; // check that this class has been flagged
              edge_flags[i] = j;  // mark the i'th edge as the first for the j'th ideal class
            }
          i++;
        }
    }
  // Now edge_flags[i] is set to 1 iff the i'th edge is {sigma,oo} where
  // sigma is the first singular point in one ideal class.

  vector<vector<INT>> invs;
  int
    GL2 = group&1, // i.e. 1 or 3
    SL2 = group&2; // i.e. 2 or 3
  if (GL2)
    {
      vector<vector<int>> M21 = face_boundary_matrix(1, debug>1);
      if (debug)
        cout << "GL2 face boundary matrix M21 has size " << M21.size() << " x " << M21[0].size() << endl;
      // Set the entries to 0 in positions j where edge_flags[j]>0
      i = 0;
      for (auto Mi : M21)
        {
          vector<int> Minew;
          Minew.reserve(n_edges);
          for (j=0; j<n_edges; j++)
            {
              if (!edge_flags[j])
                Minew.push_back(Mi[j]);
            }
          M21[i] = Minew;
          i++;
        }
      invs.push_back(invariants(M21));
    }
  if (SL2)
    {
      vector<vector<int>> M21 = face_boundary_matrix(0, debug>1);
      if (debug)
        cout << "SL2 face boundary matrix M21 has size " << M21.size() << " x " << M21[0].size() << endl;
      // Set the entries to 0 in positions j where edge_flags[j]>0
      i = 0;
      for (auto Mi : M21)
        {
          vector<int> Minew;
          Minew.reserve(n_edges);
          for (j=0; j<n_edges; j++)
            {
              if (!edge_flags[j])
                Minew.push_back(Mi[j]);
            }
          M21[i] = Minew;
          i++;
        }
      invs.push_back(invariants(M21));
    }
  SwanTimer.stop(step);
  if (showtimes)
    SwanTimer.show(1, step);
  return invs;
}

// This assumes that a is reduced mod b and that N(a)>N(b) so alpha=0 will not work
int SwanData::find_best_alpha(const Quad& a, const Quad& b, Quad& shift, int lucky) const
{
  Quad a0 = a, b0=b, b1, s, best_shift(0);
  INT best_bnorm = b.norm(), b0norm=b0.norm();
  int t=0, best_t=-1;
  // Look for a suitable alpha, trying all in turn, returning the type
  // of the one which gives best reduction
  for (const auto& M: Mlist)
    {
      if (t==0) // skipping alpha=0
        {
          t++;
          continue;
        }
      // First see if we can reduce without a further shift
      shift = 0;
      b1 = M.apply_left_den(a0,b0);
      if (b1.norm() >= b0norm) // not successful yet
        {
          // Find the shift taking a/b closest to alpha
          s = M.entry(1,0);
          shift = b1/(b0*s);      // closest integer to (a1/b1)-(r/s)
          if (!shift.is_zero())  // do the extra translation by q
            b1 = M.apply_left_den(a0-shift*b0,b0);
        }
      if (b1.norm() < best_bnorm)
        {
          best_bnorm = b1.norm();
          best_shift = shift;
          best_t     = t;
          if (lucky)
            break;
        }
      t++;
    }
  shift = best_shift;
  return best_t;
}

void SwanData::pseudo_euclidean_step(Quad& a, Quad& b, int& t,  Quad& c1, Quad& d1, Quad& c2, Quad& d2) const
{
  // We update c1,d1 unless they are both 0 and similarly c2,d2.
  // We record the type of the transformation in t unless it is initialised to -1.
  // For simple gcd, we need none of these;  for bezout (extended EA) we need c1,d1 and c2,d2
  // For convergents we need c1,d1 and t
#ifdef DEBUG_PSEA
  cout<<"Entering pseudo_euclidean_step with a="<<a<<", N(a)="<<a.norm()<<", b="<<b<<", N(b)="<<b.norm()<<endl;
#endif
  t = 0;
  if (b.is_zero())
    return;

  int compute_c1d1 = !(c1.is_zero() && d1.is_zero());
  int compute_c2d2 = !(c2.is_zero() && d2.is_zero());
  Quad u, shift = a/b;  // rounded, so N(a/b - shift) is minimal

  a.subprod(shift,b);
  if (compute_c1d1) d1.addprod(shift,c1);
  if (compute_c2d2) d2.addprod(shift,c2);
#ifdef DEBUG_PSEA
  cout<<" - translation = "<<shift<<endl;
#endif

  if (a.norm()<b.norm()) // standard Euclidean M works
    {
      u = a; a=-b; b=u;
      if (compute_c1d1) {u = -d1; d1=c1; c1=u;}
      if (compute_c2d2) {u = -d2; d2=c2; c2=u;}
#ifdef DEBUG_PSEA
      cout<<" - after inverting by S, returning (a,b) = ("<<a<<","<<b<<") ";
      if (compute_c1d1) cout << "(c1,d1)=("<<c1<<","<<d1<<") ";
      if (compute_c2d2) cout << "(c2,d2)=("<<c2<<","<<d2<<") ";
      cout <<" type=0" << endl;
#endif
      t = 0;
      return;
    }

  // If a/b is a translate of a singular point sigma, we apply an
  // extra translation if necessary and set t=-i where the translate
  // of a/b is the i'th singular point.  This is quick to check so we
  // do it first; if this fails then a/b can be reduced using at least
  // one alpha.

  t = sigma_index_with_translation(a,b, shift);
  if (t!=-1)
    {
      a.subprod(shift,b);
      if (compute_c1d1) d1.addprod(shift,c1);
      if (compute_c2d2) d2.addprod(shift,c2);
#ifdef DEBUG_PSEA
      cout<<" - a/b is singular, a translate of sigma["<<t<<"]"<<endl;
#endif
      t = -t;
      return;
    }

  // Now look for a suitable alpha (skipping 0/1 which we already tested)

  t = find_best_alpha(a, b, shift, 1); // 1 means use the first one found, 0 find the best
  if (t==-1)
    {
#ifdef DEBUG_PSEA
      cout<<" - all alphas failed thoough a/b is not singular"<<endl;
#endif
      // We should never arrive here, as it means that all alphas have
      // failed and a/b is not singular.
      cerr<<"Pseudo-Euclidean step fails for ("<<a<<", "<<b<<")"<<endl;
      exit(1);
    }

  mat22 M = Mlist[t];
  a.subprod(shift,b);
  M.apply_left(a,b);
  if (compute_c1d1)
    {
      d1.addprod(shift,c1);
      M.apply_right_inverse(c1,d1);
    }
  if (compute_c2d2)
    {
      d2.addprod(shift,c2);
      M.apply_right_inverse(c2,d2);
    }
#ifdef DEBUG_PSEA
  cout<<" - reduction success (with shift="<<shift<<" and alpha="<<alphas[t]
      <<"), returning (a,b) = ("<<a<<","<<b<<"), type "<<t<<endl;
#endif
  return;
}

// Each relation is a signed sum of edges (M)_alpha = {M(alpha},
// M(oo)} for M in the list mats and alpha=alphas[t] (when t>=0) or
// sigmas[-t] (when t<0), for t in the list types.  Here we check that
// such a relation holds identically in H_3 (not just modulo the
// congruence subgroup!)

// Special case: all signs +1
int SwanData::check_rel(const vector<mat22>& mats, const vector<int>& types) const
{
  vector<int> signs(mats.size(), 1);
  return check_rel(mats, types, signs);
}

// General case:

//#define DEBUG_FACE_RELATION

int SwanData::check_rel(const vector<mat22>& mats, const vector<int>& types, const vector<int>& signs) const
{
#ifdef DEBUG_FACE_RELATION
  int n = mats.size();
  cout<<"  - Checking "<<(n==2? "edge": "face")<<" relation...\n";
  cout<<"    mats: "<<mats<<endl;
  cout<<"    types: "<<types<<endl;
  cout<<"    signs: "<<signs<<endl;
#endif
  auto mi = mats.begin();
  auto ti = types.begin(), si = signs.begin();
  vector<RatQuad> as, bs;
  while (ti!=types.end())
    {
      mat22 M = *mi++;
      RatQuad c = base_point(*ti++);
      int s = *si++;
#ifdef DEBUG_FACE_RELATION
      cout<<"    M = "<<M<<" maps {"<< c <<",oo} to ";
#endif
      RatQuad a = M(c), b = M.image_oo();
#ifdef DEBUG_FACE_RELATION
      cout<<"{"<<a<<","<<b<<"}"<<endl;
#endif
      if (s>0) // use {a,b} when sign is +1
        {
          as.push_back(a);
          bs.push_back(b);
        }
      else // use {b,a} when sign is -1
        {
          as.push_back(b);
          bs.push_back(a);
        }
    }

  auto ai = as.begin()+1, bi = bs.begin();
  int ok=1;
  while ( bi!=bs.end() &&ok)
    {
      RatQuad next_alpha = (ai==as.end()? as[0]: *ai++);
      ok = ok && (*bi++==next_alpha);
    }
  if (!ok)
    {
      int nsides = mats.size();
      cout<<"\n*************Bad "<< (nsides==2? "edge": "face") << " relation!\n";
      cout<<"alphas: "<<as<<endl;
      cout<<"betas:  "<<bs<<endl;
      assert(ok);
    }
#ifdef DEBUG_FACE_RELATION
  else
    {
      cout<<"  - Good "<< (n==2? "edge": "face") << " relation!\n";
    }
#endif
  return ok;
}

// Generalization of extended Euclidean algorithm.

// Given a,b, returns M=[d1,-d2;-c1,c2] (the inverse of [c2,d2;c1,d1])
// such that g/h = M(a/b), i.e. g=d1*a-d2*b, h = -c1*a+c2*b, where
//
// (1) if (a,b) is principal then (a,b)=(g) and h=0, and s=0;
// (2) otherwise (a,b)=(g,h) and g/h is  the s'th singular point (s>=1).
//
// So in case (1) we get essentially the same information as
// quadbezout_psea(aa, bb, xx, yy) with xx,yy the top row of M.
//
// Note that in case 2, g/h is equal to the singular point sigma=g0/h0
// as an element of the field, but not as a fraction, since the ideal
// (g,h)=(a,b) is in the same (non-principal) ideal class as
// (g0,h0). but is not the same ideal.  In fact, (g,h)=lambda*(g0,h0)
// with lambda = g/g0 = h/h0 (since g*h0=h*g0), but in general lambda
// is not integral.

mat22 SwanData::generalised_extended_euclid(const Quad& aa, const Quad& bb, int& s) const
{
  Quad a(aa), b(bb), c1(Quad::zero), d1(Quad::one), c2(Quad::one), d2(Quad::zero);
  int t=0;
  while (b.norm()>0 && t>=0)
    {
      pseudo_euclidean_step(a, b, t, c1, d1, c2, d2);
      assert ((c2*d1-c1*d2).is_one());
      assert (c1*a+d1*b==bb);
      assert (c2*a+d2*b==aa);
    }

  // Now (1) c2*d1-c1*d2 = 1;
  //     (2) c2/c1 = aa/bb (as a reduced fraction), since b = c2*bb-c1*aa = 0;
  //     (3) a = gcd(aa,bb), since (aa,bb)=(a,b)=(a);
  //     (4) aa=a*c2, bb=a*c1;
  //     (5) aa*d1-bb*d2 = a.

  // Note the matrix inversion involved here:
  // [c2,d2;c1,d1]*[a,b] = [aa,bb], so
  // [d1,-d2;-c1,c2]*[aa,bb] = [a,b]

  s = (b.norm().is_zero()? 0 : -t);

  mat22 M(d1,-d2,-c1,c2); // maps aa/bb to a/b
  assert (d1*aa-d2*bb == a);
  assert (-c1*aa+c2*bb == b);
  return M;
}

int SwanData::check_aaa_triangle(const POLYGON& T, int verbose) const
{
  const vector<int>& tri = T.indices;
  const Quad& u = T.shifts[0];
  if (verbose)
    cout<<"Checking aaa-triangle ("<<tri<<","<<u<<")"<<"..."<<flush;
  mat22 Mi=Mlist[tri[0]], Mj=Mlist[tri[1]], Mk=Mlist[tri[2]];
  RatQuad x = (Mi(Mj.preimage_oo()+u) - Mk.preimage_oo());
  int ok = (x.is_integral());
  if (ok && verbose) cout << "OK"<<endl;
  if (!ok) cout << "BAD (x=" <<x<<")" << endl;
  return ok;
}

int SwanData::check_aas_triangle(const POLYGON& T, int verbose) const
{
  const vector<int>& tri = T.indices;
  const Quad& u = T.shifts[0];
  if (verbose)
    cout<<"Checking aas-triangle ("<<tri<<","<<u<<")"<<"..."<<flush;
  int i=tri[0], j=tri[1], k=tri[2];
  RatQuad x = Mlist[i](slist[j]+u) - slist[k];
  int ok = (x.is_integral());
  if (ok && verbose) cout << "OK"<<endl;
  if (!ok) cout << "BAD (x=" <<x<<")" << endl;
  return ok;
}

int SwanData::check_triangles(int verbose) const
{
  return
    std::all_of(T_faces.begin(), T_faces.end(),
                [this, verbose](const POLYGON& T) {return check_aaa_triangle(T, verbose);})
    &&
    std::all_of(U_faces.begin(), U_faces.end(),
                [this, verbose](const POLYGON& U) {return check_aas_triangle(U, verbose);});
}

int SwanData::check_square(const POLYGON& squ, int verbose) const
{
  // Check:  the square {{i,j,k,l},{x,y,z}} has vertices {alpha_i, oo, alpha[j']+z, beta}
  // where beta = z + M_j(x+alpha[k']) = M_i'(y+alpha_l),
  // so that M_i(T^z(M_j(x+alpha[k']))) = y+alpha_l.

  // Edges:

  // {alpha_i, oo} = (I)_i
  // {oo, alpha_j'+z} = (T^z*M_j)_j
  // {alpha_j'+z, beta} = (T^z*M_j*T^x*M_k)_k
  // {beta, alpha_i} = (M_i'*T^y)_l

  const vector<int>& ijkl = squ.indices;  // int i=ijkl[0], j=ijkl[1], k=ijkl[2], l=ijkl[3];
  const vector<Quad>& xyz = squ.shifts;
  if (verbose)
    cout<<"Checking square "<<ijkl<<", "<<xyz<<"..."<<flush;
  Quad x = xyz[0], y=xyz[1], z=xyz[2];
  mat22 Mi=Mlist[ijkl[0]], Mj=Mlist[ijkl[1]], Mk=Mlist[ijkl[2]], Ml=Mlist[ijkl[3]];
  RatQuad alpha1 = x + RatQuad(Mk.entry(0,0),Mk.entry(1,0));  // = x+alpha_k'
  RatQuad alpha2 = y + RatQuad(-Ml.entry(1,1),Ml.entry(1,0)); // = y+alpha_l
  mat22 M = Mi*mat22::Tmat(z)*Mj;
  RatQuad test = (M.entry(0,0)*alpha1+M.entry(0,1))/(M.entry(1,0)*alpha1+M.entry(1,1));
  int ok = (test == alpha2);
  if (ok && verbose) cout << "OK"<<endl;
  if (!ok) cout << "BAD (test="<<test<<", alpha2="<<alpha2<<")" << endl;
  return ok;
}

int SwanData::check_squares(int verbose) const
{
  return std::all_of(Q_faces.begin(), Q_faces.end(),
                     [this, verbose](const POLYGON& Q) {return check_square(Q, verbose);});
}

int SwanData::check_hexagon(const POLYGON& hex, int verbose) const
{
  // Check: the hexagon {{i,j,k,l,m,n},{u,x1,y1,x2,y2}} has vertices
  // {beta_1, alpha_i, oo, u+alpha[j], beta_2, gamma} where

  // beta1 = M_i'(x1+alpha[k]),
  // beta2 = M_j'(x2+alpha[l]),
  // gamma = M_i'*T^x1*M_k'*T^y1(alpha[m]) = T^u*M_j'*T^x2*M_l'*T^y2(alpha[n]).

  // Edges:

  // +{alpha_i, oo} = (I)_i
  // +{beta1, alpha_i} = (M_i'*T^x1)_k
  // +{gamma, beta1} = (M_i'*T^x1*M_k'*T^y1)_m
  // -{u+alpha_j, oo} = - (T^u)_j
  // -{beta2, alpha_j} = - (T^u*M_j'*T^x2)_l
  // -{gamma, beta2} = - (T^u*M_j'*T^x2*M_l'*T^y2)_n

  const vector<int>& ijklmn = hex.indices;
  const vector<Quad>& ux1y1x2y2 = hex.shifts;
  if (verbose)
    cout<<"Checking hexagon "<<ijklmn<<", "<<ux1y1x2y2<<"..."<<flush;
  Quad u = ux1y1x2y2[0], x1 = ux1y1x2y2[1], y1 = ux1y1x2y2[2], x2 = ux1y1x2y2[3], y2 = ux1y1x2y2[4];
  int i=ijklmn[0], j=ijklmn[1], k=ijklmn[2], l=ijklmn[3], m=ijklmn[4], n=ijklmn[5];
  RatQuad gamma1 = (Mlist[a_inv[i]]*mat22::Tmat(x1)*Mlist[a_inv[k]])(y1+alist[m]);
  RatQuad gamma2 = (mat22::Tmat(u)*Mlist[a_inv[j]]*mat22::Tmat(x2)*Mlist[a_inv[l]])(y2+alist[n]);
  int ok = gamma1==gamma2;
  if (ok && verbose) cout << "OK"<<endl;
  if (!ok) cout << "BAD (gamma1="<<gamma1<<", gamma2="<<gamma2<<")" << endl;
  return ok;
}

int SwanData::check_hexagons(int verbose) const
{
  return std::all_of(H_faces.begin(), H_faces.end(),
                     [this, verbose](const POLYGON& H) {return check_hexagon(H, verbose);});
}


