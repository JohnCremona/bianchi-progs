// FILE SWAN.CC: implementation of Swan's algorithm

#include <iostream>

#include "swan.h"
#include "geometry.h"

// Given an ideal I, return a list of singular points of class [I]
// (one representative for each orbit under integral translations).

vector<RatQuad> singular_points_in_class(Qideal I)
{
  if (I.is_principal())
    return {RatQuad::oo};
  vector<RatQuad> sigmas;
  vector<Quad> slist;
  INT Inorm = I.norm();
  Qideal Ibar = I.conj();
  Quad r, s(Inorm);
  slist.push_back(s);
  if (I!=Ibar)
    {
      Qideal I2 = I*I;
      if (I2.is_principal(s))
        {
          slist.push_back(s);
        }
    }
  for (auto si=slist.begin(); si!=slist.end(); si++)
    {
      s = *si;
      vector<Quad> rlist = Qideal(s).invertible_residues();
      for (auto ri=rlist.begin(); ri!=rlist.end(); ri++)
        {
          r = (*ri)%s;
          RatQuad sig(r,s);
          sig.normalise();
          sigmas.push_back(sig);
        }
    }
  return sigmas;
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
      vector<RatQuad> sigmas = singular_points_in_class(*I);
      sigma_list.insert(sigma_list.end(), sigmas.begin(), sigmas.end());
    }
  return sigma_list;
}

// Return sorted list of singular points (oo, denom 2, denom 3, larger denoms in +/- pairs)

vector<RatQuad> sort_singular_points(const vector<RatQuad> S)
{
  vector<RatQuad> sorted_S;
  sorted_S.push_back(RatQuad::oo);
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

  // Now process the other sigmas (if any)
  for (auto si=sigmas.begin(); si!=sigmas.end(); ++si)
    {
      RatQuad s = *si;
      // skip oo and denom 2 or 3 sigmas:
      if (s.is_infinity() or (TWO*s).is_integral() or (THREE*s).is_integral())
        continue;

      // skip sigmas if we have seen its negative:
      RatQuad ms = -s;
      if (std::find(sorted_S.begin(), sorted_S.end(), ms) != sorted_S.end())
        continue;

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
      ms = -s;
      if (i>0)
        {
          sorted_S.push_back(s);
          sorted_S.push_back(ms);
        }
      else
        {
          sorted_S.push_back(ms);
          sorted_S.push_back(s);
        }
    }
  return sorted_S;
}

