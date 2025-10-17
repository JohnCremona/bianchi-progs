// FILE nfd.cc: implementation of class Newforms for newforms of any dimension
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

newform_comparison newform_cmp;

Newform::Newform(Newforms* x, const ZZX& f, int verbose)
  :nf(x), minpoly(f)
{
  d = deg(minpoly);
  if (verbose)
    cout << "Constructing Newform from factor f = "<<polynomial_string(minpoly)<<" of degree "<<d<<endl;

  // Compute f(T); since T is scaled by dH and f(X) is not, we
  // evaluate dH^d*f(X/dH) at T; that is, we scale the coefficient of
  // X^i by dH^(d-i):
  scalar dH = nf->H1->h1cdenom();
  mat fT = evaluate(scale_poly_up(minpoly, to_ZZ(dH)), nf->T_mat);
  if (verbose)
    cout << "Computed f(T), finding its kernel..."<<flush;
  S = kernel(fT);
  int dimS=dim(S);
  denom_rel=denom(S);
  denom_abs=dH*denom_rel;
  if(dimS!=d)
    {
      cout<<"Problem: eigenspace has wrong dimension "<<dimS<<", not "<<d<<endl;
      exit(1);
    }
  scales.resize(d+1);
  scales[0]=1;
  for(int i=1; i<=d; i++)
    scales[i]=scales[i-1]*denom_rel;
  if (verbose)
    {
      cout<<"Finished constructing subspace S of dimension "<<d
          <<", denom "<<denom_rel<<", absolute denom = "<<denom_abs<<endl;
      cout<<"Computing A, the restriction of T..." <<flush;
    }
  A = transpose(restrict_mat(nf->T_mat,S));
  if(verbose)
    cout<<"done."<<endl;

  // Check that (scaled) charpoly(A) = fT
  ZZX cpA = scaled_charpoly(mat_to_mat_ZZ(A), to_ZZ(denom_abs), nf->H1->hmod);
  if (cpA!=minpoly)
    {
      cout<<endl;
      cout<<"Error: f(X) =            "<<polynomial_string(minpoly)<<endl;
      cout<<"but scaled_charpoly(A) = "<<polynomial_string(cpA)<<endl;
      exit(1);
    }

  if(verbose)
    {
      cout<<"done."<<endl;
      cout<<"A (the matrix of T restricted to S) = ";
      output_flat_matrix(A);
      if(denom_abs>1)
        cout<<" / " << denom_abs;
      cout<<endl;
      cout<<"f(X) is the min poly of A"<<endl;
    }

  // compute projcoord, precomputed projections the basis of S - the
  // double transpose is only because eclib doesn't have column
  // content and divide functions
  projcoord = transpose(nf->H1->FR.get_coord() * basis(S));
  for (int i=1; i<=projcoord.nrows(); i++)
    {
      scalar ci = projcoord.row_content(i);
      contents.push_back(ci);
      projcoord.divrow(i,ci);
    }
  projcoord = transpose(projcoord);

  // Compute basis matrix, expressing the basis on which we will
  // express eigenvalues w.r.t. the power basis on the roots of f

  mat W(d,d);
  mat Winv(d,d);
  vec v(d);  v[1]=1; // so v=[1,0,...,0]
  W.setcol(1,v);
  for(int i=2; i<=d; i++)
    {
      v = A*v;
      W.setcol(i,v);
    }
  scalar nfbasis_num = inverse(W,Winv); // so W*Winv = nfbasis_num*identity
  if(verbose>1)
    {
      cout<<"W     = ";
      output_flat_matrix(W);
      cout<<endl;
      cout<<"Winv  = ";
      output_flat_matrix(Winv);
      cout<<endl;
      cout<<"W^(-1)= (1/"<<nfbasis_num<<") * Winv"<<endl;
    }

  nfbasis = Winv;
  for(int i=0; i<d; i++)
    nfbasis.multrow(i+1,nf->Hscales[i] * scales[i]);

  // scale *columns* by contents
  for(int k=1; k<=d; k++)
    {
      scalar c = contents[k-1];
      for(int i=1; i<=d; i++)
        nfbasis(i,k) *= c;
    }
  scalar nfbasis_den = nfbasis.content();
  nfbasis /= nfbasis_den;
  if(verbose>1)
    {
      cout << "nfbasis = ";
      output_flat_matrix(nfbasis);
      cout<<endl;
    }
  nfbasis_num *= denom_abs;
  scalar g = gcd(nfbasis_num, nfbasis_den);
  nfbasis_num /= g;
  nfbasis_den /= g;
  // yes this is supposed to be den/num, those names were badly chosen!
  nfbasis_factor = bigrational(bigint(nfbasis_den), bigint(nfbasis_num));

  if (verbose)
    {
      cout<<"Basis for Hecke eigenvalues, in terms of powers of alpha:"<<endl;
      if (nfbasis_factor != bigrational(1))
        cout << "(" << nfbasis_factor << ") * ";
      output_flat_matrix(transpose(nfbasis));
      cout<<endl;
      if (d==1 && nfbasis_factor != bigrational(1))
        cout << "(eigenvalues will be scaled by "<<nfbasis_factor<<" so this can be ignored)"<<endl;
    }

  // Compute character values
  epsvec.resize(nf->eps_ops.size());
  std::transform(nf->eps_ops.begin(), nf->eps_ops.end(), epsvec.begin(),
                 [this](const matop& op){return (eps(op)>0?+1:-1);});
  int n2r = Quad::class_group_2_rank;
  if (n2r==0)
    {
      genus_char_disc = INT(1);
      return;
    }
  if (verbose)
    cout<<"genus char values: "<<epsvec<<endl;

#if(0) // this needs more work to do anything sensible
  // Compute genus character discriminant. We cannot just multiply the
  // prime discriminants for which epsvec has entry -1 since the
  // 2-torsion in the class group has (possibly) a different basis.

  // First take product over those t2ideals whose eps is -1
  Qideal I(INT(1));
  auto epsveci = epsvec.begin();
  for (auto A : nf->t2ideals)
    {
      if (*epsveci++==-1)
        I *= A;
    }
  if (verbose)
    {
      cout<<"genus char ideal: "<<ideal_label(I)<<endl;
      cout<<"prime_disc_factors: "<<Quad::prime_disc_factors<<endl;
    }
  vector<int> genus_char = I.genus_character();
  if (verbose)
    cout<<"genus char: "<<genus_char<<endl;
  genus_char_disc = discchar(genus_char);
  if (genus_char_disc>0)
    genus_char_disc = Quad::disc/genus_char_disc;
  if (verbose)
    cout<<"genus char disc: "<<genus_char_disc<<endl;
#endif
}

int Newform::trivial_char() // 1 iff all  unramified quadratic character values (if any) are +1
{
  return std::all_of(epsvec.cbegin(), epsvec.cend(), [](int i) { return i == +1; });
}

Newforms::Newforms(homspace* h1, int maxnp, int maxc, int verb)
  : verbose(verb), H1(h1)
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

  // Make the unramified quadratic character operators
  t2ideals = make_nulist(N);
  eps_ops.resize(t2ideals.size());
  std::transform(t2ideals.begin(), t2ideals.end(), eps_ops.begin(),
                 [this](Qideal& A){return CharOp(A, N);});
  if (verbose&& t2ideals.size())
    {
      cout << "Unramified quadratic character operators: ";
      for (auto T: eps_ops) cout << T.name() << " ";
      cout << endl;
    }
  // Find the splitting operator
  find_T(maxnp, maxc);
  // Construct the newforms if that succeeded
  if (split_ok)
    for (auto f: factors)
      newforms.push_back(Newform(this, f, verbose));
  else
    // abort if not
    cout << "Unable to find a suitable splitting operator!" << endl;
  // Sort the newforms (by character, dimension, polynomial)
  std::sort(newforms.begin(), newforms.end(), newform_cmp);
}

// compute T=T_P, trying all good P with N(P)<=maxnormP
void Newforms::find_T(int maxnp, int maxc)
{
  split_ok = 0;
  gmatop T_op;
  ZZX f;
  if (maxc==0)
    {
      int np = 0;
      QuadprimeLooper Pi(N); // loop over primes not dividing N
      while(np<maxnp && Pi.ok())
        {
          Quadprime P = Pi;
          if (verbose>1)
            cout << "In loop over primes, P = " << P << " with norm " << P.norm() << endl;
          np++;
          T_op = gmatop(AutoHeckeOp(P,N));
          T_name = T_op.name();
          split_ok = test_splitting_operator(N, T_op, H1->modulus, verbose);
          if (split_ok)
            break; // out of loop over primes
        }
    } // end of simple case (using a single prime)
  else // use a linear combination
    {
      vector<matop> ops = eps_ops;
      vector<Quadprime> Plist = make_goodprimes1(N, maxnp, 1); // only_one_conj=1
      vector<matop> TPlist(maxnp);
      std::transform(Plist.begin(), Plist.end(), TPlist.begin(),
                     [this](Quadprime P){return AutoHeckeOp(P,N);});
      ops.insert(ops.end(), TPlist.begin(), TPlist.end());
      // NB should enhance all_linear_combinations so that the first
      // eps_ops.size() entries are only in {0,1} as only the parity
      // matters for involutions
      vector<vector<int>> lincombs = all_linear_combinations(ops.size(), maxc);
      if (verbose)
        {
          cout << "Trying linear combinations with coefficients up to "<<maxc
               <<" of " << ops.size() << " operators";
          if (verbose>1)
            {
              cout << " (";
              // for (int i=0; i<ops.size(); i++)
              //   cout << " " <<i<<" "<< ops[i].name();
              for (auto T: ops)
                cout << " " << T.name();
              cout << ")";
            }
          cout << endl;
        }
      for (auto lc: lincombs)
        {
          vector<scalar> ilc(lc.size());
          std::transform(lc.begin(), lc.end(), ilc.begin(), [](int c){return scalar(c);});
          T_op = gmatop(ops, ilc);
          T_name = T_op.name();
          if (verbose)
            cout << "Trying "<<lc<<": "<<T_name<<"..."<<flush;
          split_ok = test_splitting_operator(N, T_op, H1->modulus, verbose>1);
          if (split_ok)
            {
              if (verbose)
                cout<<"OK!"<<endl;
              break;
            }
          if (verbose)
            cout << " no good, continuing..." << endl;
        }
    } // end if not-simple case
  if (split_ok)
    {
      if (verbose)
        cout << " OK: using operator " << T_name << " to split off newforms" << endl;
      T_mat = heckeop(T_op, 0, 1); // not cuspidal,  dual
      f = get_new_poly(N, T_op, 1, H1->modulus); // cuspidal=1
    }
  else
    return;

  factors.clear();
  if (verbose)
    {
      cout << " New cuspidal Hecke polynomial for operator" << T_name
           <<" is "<<polynomial_string(f)<<endl;
    }
  vec_ZZX NTL_factors= SFFactor(f);
  ::sort(NTL_factors.begin(), NTL_factors.end(), poly_cmp);
  int nfactors = NTL_factors.length();
  if (verbose)
    cout<<"Irreducible factors:"<<endl;
  for(int i=0; i<nfactors; i++)
    {
      ZZX fi = NTL_factors[i];
      if (verbose)
        cout<<(i+1)<<":\t"<<polynomial_string(fi)<<"\t(degree "<<deg(fi)<<")"<<endl;
      factors.push_back(fi);
    }
  return;
}

vec Newform::eig(const matop& T) const
{
  nf->H1->projcoord = projcoord;
  vec ap = nf->H1->applyop(T, nf->H1->freemods[pivots(S)[1] -1], 1); // 1: proj to S
  // if the Hecke field is Q we apply the scale factor here (if any)
  // as we don't want to display a basis.
  if (d==1 && nfbasis_factor != bigrational(1))
    {
      bigint n = num(nfbasis_factor);
      bigint d = den(nfbasis_factor);
#if (SCALAR_OPTION==3)
      ap[1] = (ap[1]*n)/d;
#else
      ap[1] = (ap[1]*I2long(n))/I2long(d);
#endif
    }
  return ap;
}

vec Newform::ap(Quadprime& P) const
{
  return eig(AutoHeckeOp(P,nf->N));
}

scalar Newform::eps(const matop& T) const // T should be a scalar
{
  scalar ans = contents[0] * eig(T)[1];
  scalar den = denom_abs;
  return ans/den;

  // if (abs(ans)==den)
  //   return scalar(ans*den>0 ? +1 : -1);

  // cout<<"Problem in eps:  unscaled value "<<ans<<" but denom_rel = "<<denom_rel
  //     <<", denom_abs = "<<denom_abs<<endl;
  // cout<<"Recomputing another way..."<<endl;
  // mat E = nf->H1->calcop_restricted(T, S, 0, 0); // dual=0, display=0
  // cout<<"matrix of "<<T.name()<<" = \n"<<E<<endl;
  // scalar s = scalar(E(1,1)>0 ? +1 : -1);
  // cout<<"so eps = "<<s<<endl;
  // return s;
}

vector<vec> Newforms::eig(const matop& T) const
{
  vector<vec> ans(newforms.size());
  std::transform(newforms.begin(), newforms.end(), ans.begin(),
                 [T](const Newform f){return f.eig(T);});
  return ans;
}

vector<vec> Newforms::ap(Quadprime& P)
{
  return eig(AutoHeckeOp(P,N));
}

// output basis for the Hecke field and character of one newform
void Newform::display(int j) const
{
  int n2r = Quad::class_group_2_rank;
  cout << "Newform " << j << endl;
  if (n2r==1)
    cout << " - Genus character value: " <<  (epsvec[0]>0? "+1": "-1") << endl;
  if (n2r>1)
    cout << " - Genus character values: " <<  epsvec << endl;
  cout << " - Dimension: "<<d<<endl;
  cout << " - Hecke field: ";
  string fpol = polynomial_string(minpoly);
  if (d==1)
    {
      cout << "Q" << endl;
    }
  else
    {
      cout << "Q(alpha) with defining polynomial "<<fpol<<" of degree "<<d;
      if (d==2)
        cout << ", discriminant "<<discriminant(minpoly);
      cout << endl;
      cout << " - Basis for Hecke eigenvalues, with respect to alpha-power basis:\n   ";
      if (nfbasis_factor != bigrational(1))
        cout << "(" << nfbasis_factor << ") * ";
      output_flat_matrix(transpose(nfbasis));
      cout<<endl;
    }
}

// output basis for the Hecke field and character of all newforms
void Newforms::display_newforms(int triv_char_only) const
{
  int j=1;
  for ( auto F : newforms)
    {
      if ((!triv_char_only) || F.trivial_char())
        {
          F.display(j);
          cout<<endl;
        }
      j++;
    }
}

vector<int> Newforms::dimensions() const
{
  vector<int> dims(newforms.size());
  std::transform(newforms.begin(), newforms.end(), dims.begin(),
                 [](const Newform& F){return F.dimension();});
  return dims;
}

mat Newforms::heckeop(const gmatop& T, int cuspidal, int dual) const
{
  return H1->calcop(T, cuspidal, dual, 0); // 0 display
}

mat Newforms::heckeop(const matop& T, int cuspidal, int dual) const
{
  return H1->calcop(T, cuspidal, dual, 0); // 0 display
}

mat Newforms::heckeop(Quadprime& P, int cuspidal, int dual)
{
  return heckeop(AutoHeckeOp(P, N), cuspidal, dual);
}

// same as m.output(cout) except no newlines between rows
void output_flat_matrix(const mat& m)
{
  vector<scalar> entries = m.get_entries();
  auto mij = entries.begin();
  cout << "[";
      cout<<flush;
  long nr = m.nrows();
  while(nr--)
    {
      long nc = m.ncols();
      cout<<"[";
      cout<<flush;
      while(nc--)
        {
          cout<<(*mij++);
          cout<<flush;
          if(nc)
            cout<<",";
          cout<<flush;
        }
      cout<<"]";
      cout<<flush;
      if(nr)
        cout<<",";
      cout<<flush;
    }
  cout << "]";
}


#if(0)
// Compute T, either one T_P or a linear combination of T_P, and its
// char poly and the irreducible factors of multiplicity 1:
void Newforms::find_T_manual()
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

void Newforms::factor_T()
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
#endif
