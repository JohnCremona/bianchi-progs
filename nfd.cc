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

HeckeField::HeckeField(const ZZX& p)
{
  if (IsMonic(p) && IsIrreducible(p))
    {
      // Set A to be the companion matrix of p.  The function
      // CompanionMatrix(p) returns a mat_ZZ so we do this manually.
      d = deg(p);
      mat_m m(d,d);
      ZZ one(1);
      for(int i=1; i<d; i++)
        {
          m(i+1,i) = one;
          m(i,d) = -coeff(p, i-1);
        }
      m(d,d) = -coeff(p, d-1);
      // Finally call the other constructor
      *this = HeckeField(m);
    }
}

HeckeField::HeckeField() // defaults to Q
{
  d=1;
  ZZ one(1);
  denom=one;
  A.init(1,1); // zero matrix
  B.init(1,1); // zero matrix
  C.init(1,1); // zero matrix
  B(1,1) = Bdet = Bfactor = one;
  Binv = B;
}

HeckeField::HeckeField(const mat_m& m, const ZZ& den, int verb)
  : d(m.nrows()), denom(den), A(m)
{
  if (verb)
    {
      cout << "----------------------------"<<endl;
      cout << "In HeckeField constructor" << endl;
    }
  minpoly = scaled_charpoly(mat_to_mat_ZZ(A), denom);

  // Compute change of basis matrix B, with column j equal to
  // denom^(n-j)*A^(j-1)v for j from 1 to d
  vec_m v(d);
  v[1] = pow(denom,d-1); // so v=[1,0,...,0]*denom^(d-1)
  B.init(d,d);
  B.setcol(1,v);
  for(int i=2; i<=d; i++)
    {
      v = A*v / denom;
      B.setcol(i,v);
    }
  B /= B.content();
  Binv.init(d,d);
  Bdet = inverse(B,Binv); // so B*Binv = Bdet*identity
  Bfactor = denom*Bdet;
  // Now we should have Binv*A*B = Bfactor * companion matrix of minpoly
  C = Binv*A*B;
  assert (Bfactor*CompanionMatrix(minpoly) == mat_to_mat_ZZ(C));
  C /= (denom*Bdet); // now C == companion matrix
  if(verb)
    {
      cout<<"inverse basis = ";
      output_flat_matrix(B);
      cout<<endl;
      cout<<"scaled basis  = ";
      output_flat_matrix(Binv);
      cout<<endl;
      cout << "basis factor  = " << Bfactor <<endl;
      if (verb>1)
        {
          cout<<"companion matrix  = ";
          output_flat_matrix(C);
          cout<<endl;
        }
      cout << "Leaving HeckeField constructor" << endl;
      cout << "----------------------------"<<endl;
    }
}

// Linear combinarion of n>0 matrices, all dxd
mat_m lin_comb_mats(const vec_m& co, const vector<mat_m>& mats)
{
  int n = mats.size(), d = mats[0].nrows();
  mat_m a(d,d);
  for (int i=0; i<n; i++)
    {
      ZZ c = co[i+1];
      if (c!=0)
        a += c*mats[i];
    }
  return a;
}

// Linear combinarion of n>0 matrices, all dxd
mat_m lin_comb_mats(const vector<ZZ>& co, const vector<mat_m>& mats)
{
  int n = mats.size(), d = mats[0].nrows();
  mat_m a(d,d);
  for (int i=0; i<n; i++)
    {
      ZZ c = co[i];
      if (c!=0)
        a += c*mats[i];
    }
  return a;
}

void HeckeField::display_bases(ostream&s) const
{
  mat_m I(mat_m::identity_matrix(d));

  s << "Powers of A (i.e. powers of alpha in A-embedding):\n";
  vector<mat_m> Apowers(d);
  Apowers[0] = I;
  for (int i=1; i<d; i++)
    Apowers[i] = A*Apowers[i-1];
  for (auto Apow: Apowers)
    {
      Apow.output(s);
      s<<endl;
      s<<endl;
    }
  ZZ fac = pow(denom, d-1) * Bfactor;
  s << "Basis in A-embedding, scaled by "<< fac <<":\n";
  s << "(first columns should be standard basis vectors * "<< (fac/denom) <<")\n";
  for(int i=1; i<=d; i++)
    {
      vec_m coli = Binv.col(i);
      for(int j=1; j<=d; j++)
        coli[j] *= pow(denom, d-j);
      mat_m M = lin_comb_mats(coli, Apowers);
      M.output(s);
      s<<endl;
      s<<endl;
      assert (denom*M.col(1)==fac*I.col(i));
    }
  s << "Powers of C (i.e. powers of alpha in C-embedding):\n";
  vector<mat_m> Cpowers(d);
  Cpowers[0] = I;
  for (int i=1; i<d; i++)
    Cpowers[i] = C*Cpowers[i-1];
  for (auto Cpow: Cpowers)
    {
      Cpow.output(s);
      s<<endl;
      s<<endl;
    }
  s << "Basis in C-embedding, scaled by "<<Bfactor<<":\n";
  s << "(first columns should be columns of basis())\n";
  for(int i=1; i<=d; i++)
    {
      mat_m M = lin_comb_mats(Binv.col(i), Cpowers);
      M.output(s);
      s<<endl;
      s<<endl;
      assert (M.col(1) == Binv.col(i));
    }
}

void HeckeField::display(ostream&s) const
{
  string fpol = polynomial_string(minpoly);
  if (d==1)
    {
      s << "Q" << endl;
      return;
    }
  s << "Q(a) with defining polynomial "<< fpol <<" of degree "<<d;
  if (d==2)
    s << ", discriminant "<<discriminant(minpoly);
  s << endl;
  s << "   Raw basis with respect to alpha-power basis:\n   ";
  s<<"[";
  for(int i=1; i<=d; i++)
    {
      if(i>1) s << ", ";
      s << polynomial_string(Binv.col(i), "a");
    }
  s<<"]";
  // s<<"\t";
  // output_flat_matrix(transpose(Binv), s);
  if (Bfactor != scalar(1))
    s << " / " << Bfactor;
  s<<endl;
}

Newform::Newform(Newforms* x, const ZZX& f, int verbose)
  :nf(x)
{
  d = deg(f);
  string fstring = polynomial_string(f);
  if (verbose)
    cout << "Constructing Newform from factor f = "<< fstring <<" of degree "<<d<<endl;

  // Compute f(T); since T is scaled by dH and f(X) is not, we
  // evaluate dH^d*f(X/dH) at T; that is, we scale the coefficient of
  // X^i by dH^(d-i):
  scalar dH = nf->H1->h1cdenom();
  mat fT = to_mat(evaluate(scale_poly_up(f, to_ZZ(dH)), nf->T_mat));
  if (verbose)
    cout << "Computed f(T), finding its kernel..."<<endl;
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
  mat A = transpose(restrict_mat(to_mat(nf->T_mat),S));
  mat_m mA = to_mat_m(A);
  if(verbose)
    cout<<"done."<<endl;

  // Check that (scaled) charpoly(A) = fT
  ZZX cpA = scaled_charpoly(mat_to_mat_ZZ(A), to_ZZ(denom_abs), nf->H1->hmod);
  if (cpA != f)
    {
      cout<<endl;
      cout<<"Error: f(X) =            "<<fstring<<endl;
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

  // compute projcoord, precomputed projections the basis of S
  projcoord = nf->H1->FR.get_coord() * basis(S);

  // Compute Hecke field basis (expressing the basis on which we will
  // express eigenvalues w.r.t. the power basis on the roots of f)

  F = HeckeField(mA, to_ZZ(denom_abs), verbose>1);

  if (verbose)
    {
      cout <<"Hecke field data:" << endl;
      F.display();
      if (verbose>1)
        F.display_bases();
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

vec_m Newform::eig(const matop& T, basis_type bt) const
{
  nf->H1->projcoord = projcoord;
  vec_m ap = to_vec_m(nf->H1->applyop(T, nf->H1->freemods[pivots(S)[1] -1], 1)); // 1: proj to S
  // if the Hecke field is Q we apply the scale factor here (if any)
  // as we don't want to display a basis.
  if (d==1 && F.Bfactor != 1)
    {
      ap /= F.Bfactor;
    }
  if (bt==basis_type::powers)
    {
      // cout << "Changing basis for eig of " << T.name() << endl;
      // cout << "Before: " << ap << endl;
      ap = F.Binv * ap;
      // cout << "After:  " << ap << endl;
    }
  return ap;
}

vec_m Newform::ap(Quadprime& P, basis_type bt) const
{
  return eig(AutoHeckeOp(P,nf->N), bt);
}

ZZ Newform::eps(const matop& T) const // T should be a scalar
{
  return eig(T)[1]/denom_abs;
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

vector<vec_m> Newforms::eig(const matop& T, basis_type bt) const
{
  vector<vec_m> ans(newforms.size());
  std::transform(newforms.begin(), newforms.end(), ans.begin(),
                 [T, bt](const Newform f){return f.eig(T, bt);});
  return ans;
}

vector<vec_m> Newforms::ap(Quadprime& P, basis_type bt)
{
  return eig(AutoHeckeOp(P,N), bt);
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
  F.display();
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

mat_m Newforms::heckeop(const gmatop& T, int cuspidal, int dual) const
{
  return to_mat_m(H1->calcop(T, cuspidal, dual, 0)); // 0 display
}

mat_m Newforms::heckeop(const matop& T, int cuspidal, int dual) const
{
  return to_mat_m(H1->calcop(T, cuspidal, dual, 0)); // 0 display
}

mat_m Newforms::heckeop(Quadprime& P, int cuspidal, int dual)
{
  return heckeop(AutoHeckeOp(P, N), cuspidal, dual);
}

// same as m.output(cout) except no newlines between rows
template<class T>
void output_flat_matrix(const Zmat<T>& m, ostream&s)
{
  vector<T> entries = m.get_entries();
  auto mij = entries.begin();
  s << "[";
  long nr = m.nrows();
  while(nr--)
    {
      long nc = m.ncols();
      s<<"[";
      while(nc--)
        {
          s<<(*mij++);
          if(nc)
            s<<",";
        }
      s<<"]";
      if(nr)
        s<<",";
    }
  s << "]";
}

template void output_flat_matrix<int>(const Zmat<int>& m, ostream&s);
template void output_flat_matrix<long>(const Zmat<long>& m, ostream&s);
template void output_flat_matrix<bigint>(const Zmat<bigint>& m, ostream&s);


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
