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
    cout<<"Computing charpoly(T)..."<<flush;
  // Compute scaled char poly of T ( = char poly of T/dH, monic in ZZ[X])
  ZZX cpT = scaled_charpoly(mat_to_mat_ZZ(T), to_ZZ(dH), hmod);
  if (verbose)
    {
      cout << "done.";
      if (verbose>1)
        cout<<" scaled char poly = "<<cpT;
      cout<<endl;
    }

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
  if (verbose)
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

// For each factor f(X), set S to be ker(f(T)) and A the restriction
// of T to S, and compute the eigenvalue basis.
void nfd::make_irreducible_subspaces()
{
  // Loop over the factors f, setting S = ker(f(T)) and A=T|S.
  for (int j = 0; j<nfactors; j++)
    {
      ZZX fj = factors[j];
      int dj = deg(fj);
      if (verbose)
        cout << "Factor "<<j+1<<" is f = "<<fj<<" of degree "<<dj<<endl;

      // Compute f(T); since T is scaled by dH and f(X) is not, we
      // evaluate dH^d*f(X/dH) at T; that is, we scale the coefficient of
      // X^i by dH^(d-i):
      mat fT = evaluate(scale_poly_up(fj, to_ZZ(dH)), T);
      if (verbose)
        cout << "Computed f(T), finding its kernel..."<<flush;
      subspace Sj = kernel(fT);
      int dimSj=dim(Sj);
      scalar dSj=denom(Sj);
      scalar dHSj=dH*dSj;
      S.push_back(Sj);
      dS.push_back(dSj);
      dHS.push_back(dHSj);
      dimS.push_back(dimSj);
      if (verbose)
        cout << "done" << endl;
      if(dimSj!=dj)
        {
          cout<<"Problem: eigenspace has wrong dimension "<<dimSj<<endl;
          exit(1);
        }
      vector<scalar> Sscalesj(dimSj+1);
      Sscalesj[0]=1;
      for(int i=1; i<=dimSj; i++)
        Sscalesj[i]=Sscalesj[i-1]*dSj;
      if (verbose>1)
        cout << "Sscales = "<<Sscalesj<<endl;
      Sscales.push_back(Sscalesj);
      if (verbose)
        {
          cout<<"Finished constructing S of dimension "<<dimSj<<endl;
          cout<<"Computing A, the restriction of T to S..." <<flush;
        }

      mat Aj = transpose(restrict_mat(T,Sj)); // matrix of T on chosen irreducible subspace of dual space

      if(verbose)
        cout<<"done."<<endl;
      A.push_back(Aj);

      // Check that (scaled) charpoly(A) = fT

      ZZX cpA = scaled_charpoly(mat_to_mat_ZZ(Aj), to_ZZ(dHSj), hmod);
      if (cpA!=fj)
        {
          cout<<"Error: f(X) =            "<<fj<<endl;
          cout<<"but scaled_charpoly(A) = "<<cpA<<endl;
        }

      if(verbose)
        {
          cout<<"S has denom "<<dSj<<", cumulative denom = "<<dHSj<<endl;
          if(dHSj>1) cout<<dHSj<<" * ";
          cout<<"A (the matrix of T restricted to S) = ";
          Aj.output_pari(cout);
          cout<<endl;
          cout<<"f(X) is the min poly of alpha, the eigenvalue of A"<<endl;
        }

      // compute projcoord, precomputed projections the basis of S

      mat projcoordj = transpose(H1->FR.get_coord() * basis(Sj));
      vector<scalar> Scontentsj;
      for (int i=1; i<=projcoordj.nrows(); i++)
        {
          scalar ci = projcoordj.row_content(i);
          Scontentsj.push_back(ci);
          projcoordj.divrow(i,ci);
        }
      Scontents.push_back(Scontentsj);
      projcoord.push_back(transpose(projcoordj));

      // Compute change of basis matrix, expressing the basis on which we
      // will express eigenvalues in terms of the power basis on the roots
      // of f(X):

      mat Wj(dimSj,dimSj);
      mat Winvj(dimSj,dimSj);
      vec v(dimSj);  v[1]=1; // so v=[1,0,...,0]
      Wj.setcol(1,v);
      for(int i=2; i<=dimSj; i++)
        {
          v = Aj*v;
          Wj.setcol(i,v);
        }
      scalar Wdetnumj = inverse(Wj,Winvj); // so W*Winv = Wdetnum*identity
      if(verbose>1)
        {
          cout<<"W     = ";
          Wj.output_pari(cout);
          cout<<endl;
          cout<<"Winv  = ";
          Winvj.output_pari(cout);
          cout<<endl;
          cout<<"W^(-1)= (1/"<<Wdetnum<<") * Winv"<<endl;
        }
      W.push_back(Wj);
      Winv.push_back(Winvj);

      mat Winvj_scaled = Winvj;
      for(int i=0; i<dimSj; i++)
        {
          Winvj_scaled.multrow(i+1,Hscales[i] * Sscales[j][i]);
        }
      // scale *columns* by Scontents
      for(int i=1; i<=dimSj; i++)
        for(int k=1; k<=dimSj; k++)
          Winvj_scaled(i,k) *= Scontents[j][k-1];
      scalar c = Winvj_scaled.content();
      Winvj_scaled /= c;
      if(verbose>1)
        {
          cout << "Winv_scaled = ";
          Winvj_scaled.output_pari(cout);
          cout<<endl;
        }
      Winv_scaled.push_back(Winvj_scaled);
      scalar Wdetdenomj = c;
      Wdetnumj*=dHSj;
      scalar g = gcd(Wdetnumj, Wdetdenomj);
      Wdetnumj /= g;
      Wdetdenomj /= g;
      if(verbose>1)
        cout<<"Wdetdenom = "<<Wdetdenom<<", Wdetnum = "<<Wdetnum<<endl;
      Wdetnum.push_back(Wdetnumj);
      Wdetdenom.push_back(Wdetdenomj);
      if (verbose)
        {
          cout<<"Basis for Hecke eigenvalues, in terms of powers of alpha:"<<endl;
          for(int i=1; i<=dimSj; i++)
            {
              cout << "("<<Wdetdenomj<<"/"<<Wdetnumj<<")*";
              cout << Winvj_scaled.col(i)<<endl;
            }
        }
    } // end of loop over factors, index j
}

// ap_vec has length dim(S)
vector<vec> nfd::ap(Quadprime& P)
{
  vector<vec> ans;
  for (int j=0; j<nfactors; j++)
    {
      H1-> projcoord = projcoord[j];
      ans.push_back(H1->applyop(HeckePOp(P,N), H1->freemods[pivots(S[j])[1] -1], 1)); // 1: proj to S
    }
  return ans;
}

void nfd::display_basis(int j) const // output basis info for subspace j (1<=j<=nfactors)
{
  cout << "Factor "<<j<<":\t";
  int deg = dimS[j-1];
  ZZX f = factors[j-1];
  if (deg==1)
    {
      cout << "Hecke field Q " << endl;
    }
  else
    {
      cout << "Hecke field Q(alpha) with defining polynomial "<<f<<" of degree "<<deg;
      if (deg==2)
        {
          cout << ", discriminant "<<discriminant(f);
        }
      cout << endl;
      cout << "Basis for Hecke eigenvalues, in terms of powers of alpha:\n";
      for(int i=1; i<=dimS[j-1]; i++)
        {
          scalar n=Wdetnum[j-1], d=Wdetdenom[j-1];
          if (d>1 || n>1)
            {
              if (n>1)
                cout << "(" << d << "/" << n << ")";
              else
                cout << d;
              cout << "*";
            }
          cout << Winv_scaled[j-1].col(i)<<endl;
        }
    }
}

mat nfd::heckeop(Quadprime& P)
{
  return H1->calcop(HeckePOp(P, N), 0, 1, 0); // 1 cuspidal, 1 transpose, 0 display
}

mat nfd::heckeop_S(Quadprime& P, const subspace& S)
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
