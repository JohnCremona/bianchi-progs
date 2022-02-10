// FILE HECKE.CC: Implemention of Hecke operators for class homspace

#include <eclib/msubspace.h>
#include <eclib/xmod.h>
#include "homspace.h"

// Strategy for operators used to automatically cut out 1-dimensional
// rational eigenspaces, in characteristic 0:

// The first n2r (=Quad::class_group_2_rank) operators are T_{A,A}
// where A runs over n2r ideals coprime to N whose classes generate
// the 2-torsion in the class group. These are involutions, but we
// only consider the +1-eigenspace since we only want to find
// eigenspaces with trivial unramified character.  When the class
// number is odd, then n2r=0 so there are none of these.  The list of
// ideals A used here is homspace::nulist, filled on homspace
// construction.

// Next come nwq operators T_{A,A}*W_Q, where Q is the i'th prime
// dividing the level to power e and [Q^e] is square (so either [W] is
// square or e is even).  By the usual abuse of notation we write W_Q
// when we really mean W_{Q^e}, and A is coprime to N such that A^2Q^e
// is principal.  The list of primes Q is in homspace::badprimes, of
// length nwq, filled on homspace construction.  When the class number
// is odd, nwq is the number of prime factors of N, but in general it
// may be smaller, even 0.  For each of these we consider +1,-1 as
// eigenvalues.

// Finally come nap operators for good primes P, where the constructor
// sets nap=20 (by default) and fills the array goodprimes with the
// first nap primes not dividing N.  The operator for P is *either*
// T_{A,A}*T_P, when the class [P] is square and A^2*P is principal;
// *or* T_{A,A}*T_{P^2} when [P] is not square, and A*P is principal.
//
// In the first case the eigenvalues considered are integers a with
// |a|<=2*sqrt(N(P)).  In the second case we use the identity
//
// T_{P^2} = (T_P)^2 + N(P)T_{P,P}
//
// to deduce that when the central character is trivial, the
// eigenvalues satisfy
//
// a_{P^2} = (a_P)^2 + N(P)
//
// so the eigenvalues we consider for T_{A,A}T_{P^2} are
// {a^2+N(P) : 0<=a<=4N(P)}.

matop homspace::h1matop(int i) // return the list of matrices defining the i'th operator
{
  assert (i>=0);
  if (i<n2r) // then we yield T_{A,A} where A is the i'th generator of the class group mod squares
    return CharOp(nulist[i], N);
  i -= n2r;
  if (i<nwq) // then we yield T_{A,A}*W_QQ where QQ is the power of Q exactly dividing N and A^2*QQ is principal
    return AtkinLehnerOp(badprimes[i], N);
  // else we yield, for P the i'th good prime,
  // either T_{A,A}*T_P if [P] is square with A^2*P principal,
  // or     T_{A,A}*T_{P^2} if [P] is not square, where A*P is principal
  i -= nwq;
  Quadprime P = goodprimes[i];
  if (P.has_square_class())
    return HeckeOp(P, N);
  else
    return HeckeSqOp(P, N);
}

// the list of possible (integer) eigenvalues for the i'th operator:
vector<long> homspace::eigrange(long i)
{
  vector<long> ans;
  assert (i>=0);

  if (i<n2r)
    {
      ans = {1};
      return ans;
    }
  i -= n2r;

  if (verbose)
    cout << "eigrange for P = " << badprimes[i] << ":\t";

  if (i<nwq)
    {
      if (characteristic==2)
        ans = {1};
      else
        ans = {-1, 1};
      if (verbose)
	cout << ans << endl;
      return ans;
    }
  i -= nwq;

  if (characteristic>0)
    {
      ans = vector<long>(characteristic);
      std::iota(ans.begin(), ans.end(), 0);
      if (verbose)
        cout << ans << endl;
      return ans;
    }

  Quadprime P = goodprimes[i];
  long normp = I2long(P.norm());
  long aplim=2;
  while (aplim*aplim<=4*normp) aplim++;
  aplim--;
  if (P.has_square_class())
    {
      ans = vector<long>(2*aplim+1);
      std::iota(ans.begin(), ans.end(), -aplim);
    }
  else // want eigs of T_{P^2} such that T_P has integral eig
    {
      ans = vector<long>(aplim+1);
      std::iota(ans.begin(), ans.end(), 0);
      for (vector<long>::iterator ai = ans.begin(); ai!=ans.end(); ++ai)
        *ai = (*ai)*(*ai) + normp;
    }
  if (verbose)
    cout << ans << endl;
  return ans;
}

vec homspace::applyop(const matop& T, const RatQuad& alpha, int proj)
{ vec ans(rk);
  if (proj) ans.init(projcoord.ncols());
  for (vector<mat22>::const_iterator mi = T.mats.begin(); mi!=T.mats.end(); ++mi)
    ans = reduce_modp(ans + chain((*mi)(alpha), proj), hmod);
  return ans;
}

vec homspace::applyop(const matop& T, const modsym& m, int proj)
{ vec ans(rk);
  if (proj) ans.init(projcoord.ncols());
  for (vector<mat22>::const_iterator mi = T.mats.begin(); mi!=T.mats.end(); ++mi)
    {
      mat22 M = *mi;
      ans = reduce_modp(ans + chain(M(m.alpha()), M(m.beta()), proj), hmod);
    }
  return ans;
}

mat homspace::calcop(const matop& T, int dual, int display)
{
  mat m(rk,rk);
  for (long j=0; j<rk; j++) if (needed[j])
     { vec colj = applyop(T,freemods[j]);
       m.setcol(j+1,colj);
     }
  if(cuspidal) m = restrict_mat(smat(m),kern).as_mat();
  if(dual) m = transpose(m);
  if (display) cout << "Matrix of " << T.name() << " = " << m;
  if (display && (dimension>1)) cout << endl;
  return m;
}

mat homspace::calcop_cols(const matop& T, const vec& jlist)
{
  int i, d = dim(jlist);
  mat m(d,rk);
  for (i=1; i<=d; i++)
    {
      int j = jlist[i];
      vec colj = applyop(T,freemods[j-1]);
      m.setcol(i,colj);
     }
  return m;
}

smat homspace::s_calcop_cols(const matop& T, const vec& jlist)
{
  int i, d = dim(jlist);
  smat m(d,rk);
  for (i=1; i<=d; i++)
    {
      int j = jlist[i];
      svec colj = applyop(T,freemods[j-1]);
      m.setrow(i,colj);
     }
  return m;
}

smat homspace::s_calcop(const matop& T, int dual, int display)
{
  smat m(rk,rk);
  for (long j=0; j<rk; j++) if (needed[j])
     { svec colj = applyop(T,freemods[j]);
       m.setrow(j+1,colj);
     }
  if(cuspidal)
    {
      m = restrict_mat(transpose(m),kern);
      if(dual) m = transpose(m);
    }
  else if(!dual) {m=transpose(m);}
  if (display)
    {
      cout << "Matrix of " << T.name() << " = ";
      if (dimension>1) cout << "\n";
      cout<<m.as_mat();
    }
  return m;
}

mat homspace::calcop_restricted(const matop& T, const subspace& s, int dual, int display)
{
  long d=dim(s);
  mat m(d,rk);
  for (long j=0; j<d; j++)
     {
       long jj = pivots(s)[j+1]-1;
       vec colj = applyop(T,freemods[jj]);
       m.setrow(j+1,colj);
     }
  if(hmod)
    m = matmulmodp(m,basis(s),hmod);
  else
    m = m*basis(s);
  if(!dual) m=transpose(m); // dual is default for restricted ops
  if (display) cout << "Matrix of " << T.name() << " = " << m;
  if (display && (dimension>1)) cout << endl;
  return m;
}

smat homspace::s_calcop_restricted(const matop& T, const ssubspace& s, int dual, int display)
{
  long d=dim(s);
  smat m(d,rk);
  for (long j=1; j<=d; j++)
     {
       long jj = pivots(s)[j];
       svec colj = applyop(T,freemods[jj-1]);
       m.setrow(j,colj);
     }
  m = mult_mod_p(m,basis(s),MODULUS);
  m.reduce_mod_p();
  if(!dual) m=transpose(m); // dual is default for restricted ops
  if (display)
    {
      cout << "Matrix of " << T.name() << " = " << m.as_mat();
      if (dimension>1) cout << endl;
    }
  return m;
}

mat homspace::opmat(int i, int dual, int verb)
{
  matop op = h1matop(i);
  if(verb)
    cout<<"Computing " << op.name() <<"...";
  mat ans = calcop(op,dual,verb);
  if (verb)
    cout<<"done."<<endl;
  return ans;
}

vec homspace::opmat_col(int i, int j, int verb)
{
  matop op = h1matop(i);
  if(verb)
    cout<<"Computing " << op.name() <<"...";
  vec ans = calcop_col(op,j);
  if (verb)
    cout<<"done."<<endl;
  return ans;
}

mat homspace::opmat_cols(int i, const vec& jlist, int verb)
{
  matop op = h1matop(i);
  if(verb)
    cout<<"Computing " << op.name() <<"...";
  mat ans = calcop_cols(op,jlist);
  if (verb)
    cout<<"done."<<endl;
  return ans;
}

smat homspace::s_opmat_cols(int i, const vec& jlist, int verb)
{
  matop op = h1matop(i);
  if(verb)
    cout<<"Computing " << op.name() <<"...";
  smat ans = s_calcop_cols(op,jlist);
  if (verb)
    cout<<"done."<<endl;
  return ans;
}

mat homspace::opmat_restricted(int i, const subspace& s, int dual, int verb)
{
  matop op = h1matop(i);
  if(verb)
    cout<<"Computing " << op.name()
        <<" restricted to subspace of dimension "<<dim(s)<<" ..."<<flush;
  mat ans = calcop_restricted(op,s,dual,verb);
  if (verb)
    cout<<"done."<<endl;
  return ans;
}

smat homspace::s_opmat(int i, int dual, int verb)
{
  matop op = h1matop(i);
  if(verb)
    cout<<"Computing " << op.name() <<"..."<<flush;
  smat ans = s_calcop(op,dual,0);
  if(verb)
    cout<<"done."<<endl;
  return ans;
}

smat homspace::s_opmat_restricted(int i, const ssubspace& s, int dual, int verb)
{
  matop op = h1matop(i);
  if(verb)
    cout<<"Computing " << op.name() <<"..."
        <<" restricted to subspace of dimension "<<dim(s)<<" ..."<<flush;
  smat ans = s_calcop_restricted(op,s,dual,0); // Automatically chooses W or T
  if(verb)
    cout<<"done."<<endl;
  return ans;
}

vec homspace::maninvector(Quadprime& P, int proj)
{
  int t = P.is_principal();
  assert (t && "P should be principal in maninvector()");
  vector<Quad> resmodp=P.residues();
  Quad p = P.gen();
  vec ans = chain(Quad::zero,p, proj), part;             // =0, but sets the right length.
  vector<Quad>::const_iterator res=resmodp.begin()+1;
  while(res!=resmodp.end())
    ans = reduce_modp(ans + chain(*res++,p, proj), hmod);
  return ans;
}

vec homspace::manintwist(const Quad& lambda, const vector<Quad>& res, vector<int> chitable, int proj)
{
  vec ans = chain(Quad::zero,lambda, proj), part;          // =0, but sets the right length.
  vector<int>::const_iterator chi=chitable.begin();
  vector<Quad>::const_iterator r=res.begin();
  while(r!=res.end())
    ans = reduce_modp(ans + (*chi++)*chain(*r++,lambda, proj), hmod);
 return ans;
}
