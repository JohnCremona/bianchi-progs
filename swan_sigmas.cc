// FILE SWAN_SIGMAS.CC: implementation of singular point functions for Swan's algorithm

#include "qideal.h"
#include "swan_sigmas.h"

CuspList denom_2_sigmas()
{
  if (Quad::class_number == 1) return {};
  int d8 = (Quad::d)%8;
  if (d8==1 || d8==5) return {{1,1,2}}; // (1+w)/2
  if (d8==2 || d8==6) return {{0,1,2}}; // w/2
  if (d8==7) return {{0,1,2}, {1,-1,2}}; // w/2, (1-w)/2
  return {};
}

CuspList denom_3_sigmas()
{
  if (Quad::class_number == 1) return {};
  int d = Quad::d;
  if (d==5 || d==6 || d==23 || d==15) return {};
  int d12 = (Quad::d)%12;
  switch (d12) {
  case 2: case 5:
    return {{1,1,3}, {-1,-1,3}, {1,-1,3}, {-1,1,3}}; // {1+w, 1-w}/3 and negs
  case 3:
    return {{1,1,3}, {-1,-1,3}}; // {w+1}/3 and negs
  case 6: case 9:
    return {{0,1,3}, {0,-1,3}}; // {w}/3 and negs
  case 11:
    return {{0,1,3}, {0,-1,3}, {-1,1,3}, {1,-1,3}}; // {w, w-1}/3 and negs
  default:
    return {};
  }
}

// Given an ideal I, return a list of singular points of class [I]
// (one representative for each orbit under integral translations).

CuspList singular_points_in_class(Qideal I, int verbose)
{
  if (I.is_principal())
    return {RatQuad::infinity()};
  INT n = I.norm();
  Quad temp, s1(n);
  vector<Quad> slist = {s1};
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
  vector<CuspList> sigma_lists(Quad::class_group.size());
  std::transform(Quad::class_group.begin(), Quad::class_group.end(), sigma_lists.begin(),
                 [] (const Qideal& I) {return singular_points_in_class(I);});
  return sigma_lists;
}

// Return one list of all singular points (excluding oo).

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
  Quad w = Quad::w, two(2), three(3);
  long d = Quad::d;

  switch (d%8) {
  case 1: case 5:
    sorted_slist.push_back(RatQuad(1+w,two));
    break;
  case 2: case 6:
    sorted_slist.push_back(RatQuad(w,two));
    break;
  case 7:
    sorted_slist.push_back(RatQuad(w,two));
    sorted_slist.push_back(RatQuad(1-w,two)); // NB not in rectangle: (w-1)/2 is
    break;
  default:
    ;
  }

  // denom 3 sigmas we construct directly
  switch (d%12) {
  case 2: case 5:
    if (d>5)
      {
        sorted_slist.push_back(RatQuad(1+w,three));
        sorted_slist.push_back(RatQuad(-1-w,three));
        sorted_slist.push_back(RatQuad(1-w,three));
        sorted_slist.push_back(RatQuad(-1+w,three));
      }
    break;
  case 3:
    if (d>15)
      {
        sorted_slist.push_back(RatQuad(1+w,three));
        sorted_slist.push_back(RatQuad(-1-w,three)); // NB not in rectangle: (2-w)/3 is
      }
    break;
  case 6: case 9:
    if (d>6)
      {
        sorted_slist.push_back(RatQuad(w,three));
        sorted_slist.push_back(RatQuad(-w,three));
      }
    break;
  case 11:
    if (d>23)
      {
        sorted_slist.push_back(RatQuad(w,three));
        sorted_slist.push_back(RatQuad(-w,three));
        sorted_slist.push_back(RatQuad(-1+w,three));
        sorted_slist.push_back(RatQuad(1-w,three));
      }
    break;
  default:
    ;
  }

  if (verbose)
    cout<<"Sigmas with small denominators: "<<sorted_slist<<endl;

  // Now process the other sigmas (if any)
  Quad temp;
  for ( auto s0 : slist) // not const or reference as we may change it
    {
      RatQuad s = reduce_to_rectangle(s0, temp);
      if (!temp.is_zero() && verbose)
        cout<<"sort_singular_points replacing "<<s0<<" by "<<s<<", its translate by "<<temp<<endl;

      if (verbose)
        cout <<"sigma = "<<s<<endl;
      assert (s.in_rectangle());

      // skip oo and denom 2 or 3 sigmas:
      if (s.is_infinity() or (two*s).is_integral() or (three*s).is_integral())
        {
          if (verbose)
            cout <<" - skipping (small denominator)"<<endl;
          continue;
        }

      RatQuad ms = -s; // do not reduce to rectangle!

      // skip sigmas if we have seen its negative:
      if (cusp_index_with_translation(ms, sorted_slist, temp)>=0)
        {
          if (verbose) cout <<" - skipping (negative seen already)";
          continue;
        }
      sorted_slist.push_back(s);
      sorted_slist.push_back(ms);
    }
  if (verbose)
    cout<<"All sigmas: "<<sorted_slist<<endl;

  return sorted_slist;
}


// Output sorted list of singular points (oo, denom 2, denom 3, larger denoms in +/- pairs)

string make_S_line(const RatQuad& s)
{
  ostringstream ost;
  ost << Quad::d << " S ";
  ost << s.num().re() << " " << s.num().im() << " ";
  ost << s.den().re() << " " << s.den().im();
  return ost.str();
}


void output_singular_points(const CuspList& S, int to_file, int to_screen)
{
  vector<Quad> small_denoms = {Quad(0), Quad(1), Quad(2), Quad(3)};
  ofstream geodata;
  stringstream ss;
  if (to_file)
    {
      ss << "geodata_" << Quad::d << ".dat";
      geodata.open(ss.str().c_str(), ios_base::app);
    }
  int nlines=0;
  CuspList sigmas_output; // record when s is output so -s is not also output
  for ( const auto& s : S)
    {
      // do not output oo or sigmas with denom <=3
      if (std::find(small_denoms.begin(), small_denoms.end(), s.den()) != small_denoms.end())
        continue;
      // do not output s if s or -s already output
      if (std::find(sigmas_output.begin(), sigmas_output.end(), s) != sigmas_output.end())
        continue;
      nlines++;
      sigmas_output.push_back(s);
      sigmas_output.push_back(-s);
      string st = make_S_line(s);
      if (to_file)
        geodata << st << endl;
      if (to_screen)
        cout << st << endl;
    }
  if (to_file)
    geodata.close();
  if (to_screen)
    cout << nlines << " S lines output" <<endl;
}


vector<RatQuad> test_singular_points(int output_level)
{
  if (output_level>=3)
    for (auto I : Quad::class_group)
      {
        cout<<"Ideal class ["<<I<<"]: ";
        cout<<"singular points "<<singular_points_in_class(I,(output_level>3))<<endl;
      }
  auto sigs = singular_points();
  if (output_level>=1)
    cout << "Number of singular points, including oo: "<<sigmas.size()<<endl;
  if (output_level>=2)
    cout << "Unsorted singular points: "<<sigmas<<endl;
  sigs = sort_singular_points(sigs);
  if (output_level>=1)
    cout << "Sorted singular points: "<<sigmas<<endl;
  int to_file=0; //(output_level>=1);
  int to_screen=0; //(output_level>=2);
  output_singular_points(sigs, to_file, to_screen);
  return sigs;
}
