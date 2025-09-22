// FILE nfd.cc: implementation of class nfd for newforms of any dimension
//////////////////////////////////////////////////////////////////////////
//
// Adapted from the similar class (over Q) in eclib
//
//////////////////////////////////////////////////////////////////////////

#include <NTL/mat_ZZ.h>
#include <NTL/mat_poly_ZZ.h>
#include <NTL/ZZXFactoring.h>
#include "eclib/subspace.h"
#include "matprocs.h"
#include "nfd.h"

#define OUTPUT_PARI_STYLE

nfd::nfd(homspace* h1, int verb)
  : H1(h1), verbose(verb)
{
  N = H1->N;
  dimH = H1->h1dim();
  cdimH = H1->h1cuspdim();
  dH = H1->h1cdenom();
  hmod = H1->hmod;

  if(verbose && dH>1)
    cout<<"H has dimension "<<dimH<<", cuspidal dimension "<<cdimH<<", denominator "<<dH<<endl;
  Hscales.resize(dimH+1);
  Hscales[0]=1;
  for(int i=1; i<=dimH; i++) Hscales[i]=Hscales[i-1]*dH;
  if (verbose>1)
    cout << "Hscales = "<<Hscales<<endl;
}

// Compute T, either one T_P or a linear combination of T_P, and its
// char poly and the irreducible factors of multiplicity 1:
void nfd::make_T()
{
  Quadprime P;
  int one_p;
  cout << "Use just one prime (1) or a linear combination (0)? ";
  cin >> one_p;
  if(one_p) // Compute just one Tp:
    {
      cout<<"Enter a prime P (label or generator): ";
      cin>>P;
      // QuadprimeLooper L(N); // loop over primes not dividing N
      // P = L; // first one
      if (verbose)
        cout << "Computing T_P for P = " << P << "..." << flush;
      T = heckeop(P);
      if (verbose)
        cout<<"done."<<endl;
    }
  else // a linear combination:
    {
      int nP;
      scalar cP;
      cout<<"Enter a linear combination of I and T_P for one or more primes P.\n";
      cout<<"First enter the coefficient of the identity: ";
      cin>>cP;
      T = mat::scalar_matrix(dimH, cP);
      cout<<"Now enter the number of P: "; cin>>nP;
      for (int iP=0; iP<nP; iP++)
        {
          cout<<"Enter a prime P (label or generator): ";
          cin>>P;
          cout<<"Enter the coefficient of "<<opname(P,N)<<": ";
          cin>>cP;
          if(verbose)
            cout << "Computing "<<opname(P,N)<<" for P = " << ideal_label(P) << "..." << flush;
          mat TP = heckeop(P);
          if(verbose)
            cout<<"done."<<endl;
          T += cP*TP;
        }
    }
  if (verbose)
    cout<<"Computing charpoly(T)... to_ZZ(dH)="<<to_ZZ(dH)<<endl;
  // Compute scaled char poly of T ( = char poly of T/dH, monic in ZZ[X])
  ZZX cpT = scaled_charpoly(mat_to_mat_ZZ(T), to_ZZ(dH), hmod);
  if (verbose)
    cout<<"(scaled) char poly = "<<cpT<<endl;

  // factor the charpoly, just the factors of multiplicity 1:
  vec_pair_ZZX_long factors_with_multiplicities;
  SquareFreeDecomp(factors_with_multiplicities,cpT);
  if(verbose>1)
    cout<<"NTL char poly square-free factors = "<<factors_with_multiplicities<<endl;
  if(factors_with_multiplicities[0].b>1)
    {
      cout<<"No factors of multiplicity 1"<<endl;
      return;
    }
  factors.clear();
  vec_pair_ZZX_long NTL_factors;
  ZZ cont;
  cout<<"Irreducible factors of multiplicity 1 are:"<<endl;
  factor(cont,NTL_factors,factors_with_multiplicities[0].a);
  ::sort(NTL_factors.begin(), NTL_factors.end(), fact_cmp);
  nfactors = NTL_factors.length();
  for(int i=0; i<nfactors; i++)
    {
      ZZX fi = NTL_factors[i].a;
      cout<<(i+1)<<":\t"<<fi<<"\t(degree "<<deg(fi)<<")"<<endl;
      factors.push_back(fi);
    }
  return;
}

// Select one factor f(X), set S to be ker(f(T)) and A the restriction
// of T to S, and compute the eigenvalue basis.  Return 0 if no factor
// is chosen.
int nfd::make_S()
{
  // Select one factor to determine the subspace S.
  // If the chosen factor is f then S = ker(f(A)).

  int j = 0;
  while((j<1)||(j>nfactors))
    {
      cout<<"Enter factor number, between 1 and "<<nfactors<<" (or 0 to stop): ";
      cin>>j;
      if (j==0)
        return 0;
    }
  f = factors[j-1];
  int d = deg(f);
  if (verbose)
    cout << "Factor "<<j<<" is f = "<<f<<" of degree "<<d<<endl;

  // Compute f(T); since T is scaled by dH and f(X) is not, we
  // evaluate dH^d*f(X/dH) at T; that is, we scale the coefficient of
  // X^i by dH^(d-i):
  mat fT = evaluate(scale_poly_up(f, to_ZZ(dH)), T);
  if (verbose)
    cout << "Computed f(T), finding its kernel..."<<flush;
  S = kernel(fT);
  dimS=dim(S);
  dS=denom(S);
  dHS=dH*dS;
  if (verbose)
    cout << "done" << endl;
  if(dimS!=d)
    {
      cout<<"Problem: eigenspace has wrong dimension "<<dimS<<endl;
      exit(1);
    }
  Sscales.resize(dimS+1);
  Sscales[0]=1;
  for(int i=1; i<=dimS; i++)
    Sscales[i]=Sscales[i-1]*dS;
  if (verbose>1)
    cout << "Sscales = "<<Sscales<<endl;

  if (verbose)
    {
      cout<<"Finished constructing S of dimension "<<dimS<<endl;
      cout<<"Computing A, the restriction of T to S..." <<flush;
    }

  A = transpose(restrict_mat(T,S)); // matrix of T on chosen irreducible subspace of dual space

  if(verbose)
    {
      cout<<"done."<<endl;
    }

  // Check that (scaled) charpoly(A) = fT

  ZZX cpA = scaled_charpoly(mat_to_mat_ZZ(A), to_ZZ(dHS), hmod);
  if (cpA!=f)
    {
      cout<<"Error: f(X) =            "<<f<<endl;
      cout<<"but scaled_charpoly(A) = "<<cpA<<endl;
    }

  if(verbose)
    {
      cout<<"S has denom "<<dS<<", cumulative denom = "<<dHS<<endl;
      if(dHS>1) cout<<dHS<<" * ";
      cout<<"A (the matrix of T restricted to S) = ";
      A.output_pari(cout);
      cout<<endl;
      cout<<"f(X) is the min poly of alpha, the eigenvalue of A"<<endl;
    }

  // compute projcoord, precomputed projections the basis of S

  mat projcoord = transpose(H1->FR.get_coord() * basis(S));
  if (verbose>1)
    cout<<"Before removing contents, projcoord = "<<projcoord<<endl;
  Scontents.clear();
  for (int i=1; i<=projcoord.nrows(); i++)
    {
      scalar ci = projcoord.row_content(i);
      Scontents.push_back(ci);
      if (verbose>1)
        cout << "Column "<<i<<" of projcoord has content "<<ci<<endl;
      projcoord.divrow(i,ci);
    }
  if (verbose>1)
    cout<<"After removing contents "<<Scontents<<",\nprojcoord = "<<projcoord<<endl;
  H1-> projcoord = transpose(projcoord);

  // Compute change of basis matrix, expressing the basis on which we
  // will express eigenvalues in terms of the power basis on the roots
  // of f(X):

  W.init(dimS,dimS);
  Winv.init(dimS,dimS);
  vec v(dimS);  v[1]=1; // so v=[1,0,...,0]
  W.setcol(1,v);
  for(int i=2; i<=dimS; i++)
    {
      v = A*v;
      W.setcol(i,v);
    }
  Wdetnum = inverse(W,Winv); // so W*Winv = Wdetnum*identity
  if(verbose>1)
    {
      cout<<"W     = ";
      W.output_pari(cout);
      cout<<endl;
      cout<<"Winv  = ";
      Winv.output_pari(cout);
      cout<<endl;
      cout<<"W^(-1)= (1/"<<Wdetnum<<") * Winv"<<endl;
    }

  Winv_scaled=Winv;
  for(int i=0; i<dimS; i++)
    {
      Winv_scaled.multrow(i+1,Hscales[i] * Sscales[i]); // * Scontents[i]);
    }
  // scale *columns* by Scontents
  for(int i=1; i<=dimS; i++)
    for(int j=1; j<=dimS; j++)
      Winv_scaled(i,j ) *= Scontents[j-1];
  if(verbose>1)
    {
      cout << "Winv_scaled = ";
      Winv_scaled.output_pari(cout);
      cout<<endl;
    }
  scalar c = Winv_scaled.content();
  Winv_scaled /= c;
  Wdetdenom = c;
  Wdetnum*=dHS;
  scalar g = gcd(Wdetnum, Wdetdenom);
  Wdetnum /= g;
  Wdetdenom /= g;
  if(verbose>1)
    cout<<"Wdetdenom = "<<Wdetdenom<<", Wdetnum = "<<Wdetnum<<endl;
  cout<<"Basis for Hecke eigenvalues, in terms of powers of alpha:"<<endl;
  for(int i=1; i<=dimS; i++)
    {
      cout << "("<<Wdetdenom<<"/"<<Wdetnum<<")*";
      cout << Winv_scaled.col(i)<<endl;
    }
  return 1;
}

// ap_vec has length dim(S)
vec nfd::ap(Quadprime& P)
{
  return H1->applyop(HeckePOp(P,N), H1->freemods[pivots(S)[1] -1], 1); // 1: proj to S
}

mat nfd::heckeop(Quadprime& P)
{
  return H1->calcop(HeckePOp(P, N), 0, 1, 0); // 1 cuspidal, 1 transpose, 0 display
}

mat nfd::heckeop_S(Quadprime& P)
{
  return H1->calcop_restricted(HeckePOp(P, N), S, 1, 0); // 1 transpose, 0 display
}

map<Qideal,homspace*> H1_dict;
map<pair<Qideal,Quadprime>, ZZX> full_poly_dict;
map<pair<Qideal,Quadprime>, ZZX> new_poly_dict;

int is_class_number_one(long d)
{
  return std::find(class_number_one_fields.begin(), class_number_one_fields.end(), d)
    != class_number_one_fields.end();
}

homspace* get_homspace(const Qideal& N, scalar mod)
{
  auto res = H1_dict.find(N);
  if (res==H1_dict.end())
    {
      homspace* H = new homspace(N, mod, 1); // cuspidal=1
      H1_dict[N] = H;
      return H;
    }
  else
    return H1_dict[N];
}

ZZX get_full_poly(const Qideal& N,  Quadprime& P, const scalar& mod)
{
  pair<Qideal,Quadprime> NP = {N,P};
  auto res = full_poly_dict.find(NP);
  if (res==full_poly_dict.end())
    {
      homspace* H = get_homspace(N, mod);
      ZZX full_poly = H->charpoly(HeckePOp(P, N), 1); // 1 for cuspidal
      full_poly_dict[NP] = full_poly;
      return full_poly;
    }
  else
    return full_poly_dict[NP];
}

ZZX get_new_poly(Qideal& N, Quadprime& P, const scalar& mod)
{
  pair<Qideal,Quadprime> NP = {N,P};
  auto res = new_poly_dict.find(NP);
  if (res==new_poly_dict.end())
    {
      ZZX new_poly = get_full_poly(N, P, mod);
      vector<Qideal> DD = alldivs(N);
      for( auto D : DD)
        {
          if (D==N)
            continue;
          ZZX new_poly_D = get_new_poly(D, P, mod);
          Qideal M = N/D;
          int mult = alldivs(M).size();
          for (int i=0; i<mult; i++)
            new_poly /= new_poly_D;
        }
      new_poly_dict[NP] = new_poly;
      return new_poly;
    }
  else
    return new_poly_dict[NP];
}
