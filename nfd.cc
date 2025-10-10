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
  Ndivs = alldivs(N);

  if(verbose && dH>1)
    cout<<"H has dimension "<<dimH<<", cuspidal dimension "<<cdimH<<", denominator "<<dH<<endl;
  Hscales.resize(dimH+1);
  Hscales[0]=1;
  for(int i=1; i<=dimH; i++) Hscales[i]=Hscales[i-1]*dH;
  if (verbose>2)
    cout << "Hscales = "<<Hscales<<endl;
}

// Compute T, either one T_P or a linear combination of T_P, and its
// char poly and the irreducible factors of multiplicity 1:
void nfd::find_T()
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
      T_mat = heckeop(P); // not cuspidal or dual
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
      T_mat = mat::scalar_matrix(dimH, cP);
      cout<<"Now enter the number of P: "; cin>>nP;
      for (int iP=0; iP<nP; iP++)
        {
          cout<<"Enter a prime P (label or generator): ";
          cin>>P;
          cout<<"Enter the coefficient of "<<opname(P,N)<<": ";
          cin>>cP;
          if(verbose)
            cout << "Computing "<<opname(P,N)<<" for P = " << ideal_label(P) << "..." << flush;
          mat TP = heckeop(P); // not cuspidal or dual
          if(verbose)
            cout<<"done."<<endl;
          T_mat += cP*TP;
        }
    }
  factor_T();
  return;
}

void nfd::factor_T()
{
  if (verbose)
    cout<<"Computing charpoly(T)..."<<flush;
  // Compute scaled char poly of T ( = char poly of T/dH, monic in ZZ[X])
  ZZX cpT = scaled_charpoly(mat_to_mat_ZZ(T_mat), to_ZZ(dH), hmod);
  if (verbose)
    {
      cout << "done (degree = "<<deg(cpT)<<").";
      if (verbose>1)
        cout<<" scaled char poly = "<<polynomial_string(cpT);
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
  factor(cont,NTL_factors,factors_with_multiplicities[0].a);
  ::sort(NTL_factors.begin(), NTL_factors.end(), fact_cmp);
  nfactors = NTL_factors.length();
  cout<<nfactors<<" irreducible factors of multiplicity 1 (may include non-cuspidal):"<<endl;
  for(int i=0; i<nfactors; i++)
    {
      ZZX fi = NTL_factors[i].a;
      cout<<(i+1)<<":\t"<<polynomial_string(fi)<<"\t(degree "<<deg(fi)<<")"<<endl;
      factors.push_back(fi);
    }
}

// compute T=T_P, trying all good P with N(P)<=maxnormP
int nfd::find_T_auto(INT maxnormP, Quadprime& P0, int verb)
{
  ZZX f;
  for ( auto& P : Quadprimes::list)
    {
      if (P.divides(N))
        continue; // to next prime
      if (P.norm() > maxnormP)
        {
          if (verb)
            cout << "No suitable splitting prime P found of norm up to "<<maxnormP<<endl;
          return 0; // give up
        }
      T_op = AutoHeckeOp(P,N);
      string T_name = T_op.name();
      if (verb)
        cout << "Trying P = " << ideal_label(P) << ", using " << T_name << "..." << flush;
      f = get_new_poly(N, T_op, H1->modulus);
      if (!IsSquareFree(f))
        {
          if (verb)
            cout << " NO: Hecke polynomial is not squarefree" << endl;
          continue; // to next prime
        }
      // We still need to check that f is coprime to the new polys
      // for P at lower levels (which have been computed and
      // cached so this is cheap)
      if (verb)
        {
          cout << " OK so far: new cuspidal Hecke polynomial for "<< T_name
               <<" is "<<polynomial_string(f)<<", which is squarefree." << endl;
          cout << " Now checking whether it is coprime to old polys..."<<endl;
        }
      int ok = 1;
      for( auto D : Ndivs)
        {
          if (D==N)
            continue;
          ZZX f_D = get_new_poly(D, AutoHeckeOp(P,D), H1->modulus); // from cache
          if (!AreCoprime(f, f_D))
            {
              if (verb)
                {
                  cout << " No: shares factors with poly at level "<<ideal_label(D) << endl;
                }
              ok = 0;
              break; // out of loop over divisors
            }
        }
      if (ok) // all proper divisors' new polys are coprime to f
        {
          if (verb)
            {
              cout << " OK: coprime to all polys at proper divisor levels" << endl;
            }
          T_mat = heckeop(P, 0, 1); // not cuspidal,  dual
          P0 = P;
          break; // out of loop over primes
        }
    }
  factors.clear();
  if (verb)
    {
      cout << " OK: new cuspidal Hecke polynomial for P="<<ideal_label(P0)
           <<" is "<<polynomial_string(f)<<", which is squarefree" << endl;
    }
  vec_ZZX NTL_factors= SFFactor(f);
  ::sort(NTL_factors.begin(), NTL_factors.end(), poly_cmp);
  nfactors = NTL_factors.length();
  if (verb)
    cout<<"Its irreducible factors of multiplicity 1 are:"<<endl;
  for(int i=0; i<nfactors; i++)
    {
      ZZX fi = NTL_factors[i];
      if (verb)
        cout<<(i+1)<<":\t"<<polynomial_string(fi)<<"\t(degree "<<deg(fi)<<")"<<endl;
      factors.push_back(fi);
    }
  //factor_T();
  return 1;
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
        cout << "Factor "<<j+1<<" is f = "<<polynomial_string(fj)<<" of degree "<<dj<<endl;

      // Compute f(T); since T is scaled by dH and f(X) is not, we
      // evaluate dH^d*f(X/dH) at T; that is, we scale the coefficient of
      // X^i by dH^(d-i):
      mat fT = evaluate(scale_poly_up(fj, to_ZZ(dH)), T_mat);
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
      if(dimSj!=dj)
        {
          cout<<"Problem: eigenspace has wrong dimension "<<dimSj<<endl;
          exit(1);
        }
      vector<scalar> Sscalesj(dimSj+1);
      Sscalesj[0]=1;
      for(int i=1; i<=dimSj; i++)
        Sscalesj[i]=Sscalesj[i-1]*dSj;
      if (verbose>2)
        cout << "Sscales = "<<Sscalesj<<endl;
      Sscales.push_back(Sscalesj);
      if (verbose)
        {
          cout<<"Finished constructing S of dimension "<<dimSj<<endl;
          cout<<"Computing A, the restriction of T to S..." <<flush;
        }

      mat Aj = transpose(restrict_mat(T_mat,Sj)); // matrix of T on chosen irreducible subspace of dual space

      if(verbose)
        cout<<"done."<<endl;
      A.push_back(Aj);

      // Check that (scaled) charpoly(A) = fT

      ZZX cpA = scaled_charpoly(mat_to_mat_ZZ(Aj), to_ZZ(dHSj), hmod);
      if (cpA!=fj)
        {
          cout<<"Error: f(X) =            "<<polynomial_string(fj)<<endl;
          cout<<"but scaled_charpoly(A) = "<<polynomial_string(cpA)<<endl;
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

      mat m1 = H1->FR.get_coord();
      //cout<<"H1->FR.get_coord() has size "<<m1.nrows()<<" x "<<m1.ncols()<<endl;
      mat m2 = basis(Sj);
      //cout<<"basis(Sj)          has size "<<m2.nrows()<<" x "<<m2.ncols()<<endl;
      mat projcoordj = transpose(m1*m2);
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
          cout<<"W^(-1)= (1/"<<Wdetnumj<<") * Winv"<<endl;
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

  // loop over unramified quadratic characters, finding the sign of each space for each
  int n2r = Quad::class_group_2_rank;
  if (n2r)
    {
      if (verbose)
        cout << "Unramified quadratic character values" << endl;
      vector<Qideal> t2ideals = make_nulist(N);
      for ( auto& A : t2ideals)
        {
          matop chi = CharOp(A, N);
          string chiname = chi.name();
          vector<scalar> epsvec = eps(chi);
          epsvecs.push_back(epsvec);
          epsnames.push_back(chiname);
          if (verbose)
            cout<<"chi(" << chiname <<") : " <<epsvec<<endl;
        }
    }
}

// ap_vec has length dim(S)
vector<vec> nfd::ap(Quadprime& P)
{
  return eig(AutoHeckeOp(P,N));
}

vector<vec> nfd::eig(const matop& T)
{
  vector<vec> ans;
  for (int j=0; j<nfactors; j++)
    {
      H1-> projcoord = projcoord[j];
      ans.push_back(H1->applyop(T, H1->freemods[pivots(S[j])[1] -1], 1)); // 1: proj to S
    }
  return ans;
}

vector<scalar> nfd::eps(const matop& T) // T should be a scalar
{
  vector<vec> epsvec = eig(T);
  auto e = [this, epsvec](int i) {return epsvec[i][1] / Sscales[i][1];};
  vector<scalar> ans(nfactors);
  int i=0;
  std::generate(ans.begin(), ans.end(), [&i, e]{return e(i++);});
  return ans;
}

void nfd::display_basis(int j) const // output basis info for subspace j (1<=j<=nfactors)
{
  cout << "Factor "<<j<<": ";
  int deg = dimS[j-1];
  ZZX f = factors[j-1];
  string fpol = polynomial_string(f);
  string tab = "          "; // enough for single digit j
  int n2r = epsvecs.size();
  vector<scalar> epsvec;
  for(int i=0; i<n2r; i++)
    epsvec.push_back(epsvecs[i][j-1]);
  auto showchar = [n2r, epsvec, tab] {
    if (n2r==1)
      {
        cout << tab << "Character: " <<  (epsvec[0]>0? "+1": "-1") << endl;
      }
    else
      if (n2r>1)
        {
          cout << tab << "Characters: " <<  epsvec << endl;
        }
  };
  if (deg==1)
    {
      cout << "Hecke field Q" << endl;
      if (n2r)
        showchar();
    }
  else
    {
      cout << "Hecke field Q(alpha) with defining polynomial "<<fpol<<" of degree "<<deg;
      if (deg==2)
        {
          cout << ", discriminant "<<discriminant(f);
        }
      cout << endl;
      showchar();
      cout << endl;
      cout << tab << "Basis for Hecke eigenvalues, in terms of powers of alpha:\n";
      for(int i=1; i<=dimS[j-1]; i++)
        {
          cout << tab;
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

mat nfd::heckeop(const matop& T, int cuspidal, int dual)
{
  return H1->calcop(T, cuspidal, dual, 0); // 0 display
}

mat nfd::heckeop(Quadprime& P, int cuspidal, int dual)
{
  return heckeop(AutoHeckeOp(P, N), cuspidal, dual);
}

int is_class_number_one(long d)
{
  return std::find(class_number_one_fields.begin(), class_number_one_fields.end(), d)
    != class_number_one_fields.end();
}
