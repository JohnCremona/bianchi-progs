// FILE SWAN.CC: implementation of Swan's algorithm

#include <iostream>

#include "swan.h"
#include "geometry.h"

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
      ss << "geodata/geodata_" << Quad::d << ".dat";
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
