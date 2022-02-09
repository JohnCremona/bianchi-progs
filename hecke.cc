// FILE HECKE.CC: Implemention of Hecke operators for class homspace

#include <eclib/msubspace.h>
#include <eclib/xmod.h>
#include "homspace.h"

// When the class number is odd, the i'th operator for 0 <= i < nap
// means either W_Q where Q is the i'th prime dividing the level, if i
// < nbp = #{bad primes}, or T_P where P is the (i-nbp)'th good prime.
// The list of these Q then P is in homspace::primelist.

// In general, (only relevant when the class number is even) we have 0
// <= i <= n2r+nap where n2r is the 2-rank of the class group; the
// first n2r operators are the involutions T_{A,A} for the n2r ideals
// A in homspace::nulist (coprime to N and generating the 2-torsion in
// the class group)

// the list of possible (integer) eigenvalues for the i'th operator:
vector<long> homspace::eigrange(long i)
{
  vector<long> ans;
  if((i<0)||(i>=nap+n2r))   // shouldn't happen
    {
      cerr << "Error in eigrange(): i="<<i<<" but should be between 0 and "<<nap+n2r<<endl;
      exit(1);
    }
  if (i<n2r)
    {
      if (characteristic==2)
        ans = {1};
      else
        ans = {-1, 1};
      return ans;
    }
  Quadprime P = primelist[i-n2r];
  long normp = I2long(P.norm());
  if (verbose)
    cout << "eigrange for P = " << P << ":\t";
  if(P.divides(N))
    {
      if (characteristic==2)
        ans = {1};
      else
        ans = {-1, 1};
      if (verbose)
	cout << ans << endl;
      return ans;
    }
  else
    {
      if (characteristic==0)
        {
          long aplim=2;
          while (aplim*aplim<=4*normp) aplim++;
          aplim--;
          if(verbose)
            cout << "|ap| up to "<<aplim<<":\t";
          ans = vector<long>(2*aplim+1);
          std::iota(ans.begin(), ans.end(), -aplim);
          if (verbose)
            cout << ans << endl;
          return ans;
        }
      else
        {
          ans = vector<long>(characteristic);
          std::iota(ans.begin(), ans.end(), 0);
          if (verbose)
            cout << ans << endl;
          return ans;
        }
    }
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
  if((i<0)||(i>=nap+n2r))   // shouldn't happen
    {
      cerr << "Error in opmat(): i="<<i<<" but should be between 0 and "<<nap+n2r<<endl;
      exit(1);
    }
  if (i<n2r)
    {
      return nu(i, dual, verb);
    }
  Quadprime P = primelist[i-n2r];
  if(verb)
    cout<<"Computing " << opname(P,N) <<"...";
  mat ans = heckeop(P,dual,verb); // Automatically chooses W or T
  if (verb)
    cout<<"done."<<endl;
  return ans;
}

vec homspace::opmat_col(int i, int j, int verb)
{
  if((i<0)||(i>=nap+n2r))   // shouldn't happen
    {
      cerr << "Error in opmat_col(): i="<<i<<" but should be between 0 and "<<nap+n2r<<endl;
      exit(1);
    }
  if (i<n2r)
    {
      return calcop_col(CharOp(nulist[i], N), j);
    }
  Quadprime P = primelist[i-n2r];
  if(verbose)
    cout<<"Computing " << opname(P,N) <<"...";
  vec ans = heckeop_col(P,j,verb); // Automatically chooses W or T
  if (verb)
    cout<<"done."<<endl;
  return ans;
}

mat homspace::opmat_cols(int i, const vec& jlist, int verb)
{
  if((i<0)||(i>=nap+n2r))   // shouldn't happen
    {
      cerr << "Error in opmat_cols(): i="<<i<<" but should be between 0 and "<<nap+n2r<<endl;
      exit(1);
    }
  if (i<n2r)
    {
      return calcop_cols(CharOp(nulist[i], N), jlist);
    }
  Quadprime P = primelist[i-n2r];
  if(verbose)
    cout<<"Computing " << opname(P,N) <<"...";
  mat ans = heckeop_cols(P,jlist,verb); // Automatically chooses W or T
  if (verb)
    cout<<"done."<<endl;
  return ans;
}

smat homspace::s_opmat_cols(int i, const vec& jlist, int verb)
{
  if((i<0)||(i>=nap+n2r))   // shouldn't happen
    {
      cerr << "Error in s_opmat_cols(): i="<<i<<" but should be between 0 and "<<nap+n2r<<endl;
      exit(1);
    }
  if (i<n2r)
    {
      return s_calcop_cols(CharOp(nulist[i], N), jlist);
    }
  Quadprime P = primelist[i-n2r];
  if(verb)
    cout<<"Computing " << opname(P,N) <<"...";
  smat ans = s_heckeop_cols(P,jlist,verb); // Automatically chooses W or T
  if (verb)
    cout<<"done."<<endl;
  return ans;
}

mat homspace::opmat_restricted(int i, const subspace& s, int dual, int verb)
{
  if((i<0)||(i>=nap+n2r))   // shouldn't happen
    {
      cerr << "Error in opmat_restricted(): i="<<i<<" but should be between 0 and "<<nap+n2r<<endl;
      exit(1);
    }
  if (i<n2r)
    {
      return calcop_restricted(CharOp(nulist[i], N), s, dual, verb);
    }
  Quadprime P = primelist[i-n2r];
  if(verb)
    cout<<"Computing " << opname(P,N)
        <<" restricted to subspace of dimension "<<dim(s)<<" ..."<<flush;
  return heckeop_restricted(P,s,dual,verb); // Automatically chooses W or T
}

smat homspace::s_opmat(int i, int dual, int v)
{
  if((i<0)||(i>=nap+n2r))   // shouldn't happen
    {
      cerr << "Error in s_opmat(): i="<<i<<" but should be between 0 and "<<nap+n2r<<endl;
      exit(1);
    }
  if (i<n2r)
    {
      return s_calcop(CharOp(nulist[i], N), dual, v);
    }
  Quadprime P = primelist[i-n2r];
  if(v)
    {
      cout<<"Computing " << opname(P,N) <<"...";
    }
  smat ans = s_heckeop(P,dual,0); // Automatically chooses W or T
  if(v)
    {
      cout<<"done."<<endl;
    }
  return ans;
}

smat homspace::s_opmat_restricted(int i, const ssubspace& s, int dual, int v)
{
  if((i<0)||(i>=nap+n2r))   // shouldn't happen
    {
      cerr << "Error in s_opmat_restricted(): i="<<i<<" but should be between 0 and "<<nap+n2r<<endl;
      exit(1);
    }
  if (i<n2r)
    {
      return s_calcop_restricted(CharOp(nulist[i], N), s, dual, v);
    }
  Quadprime P = primelist[i-n2r];
  if(v)
    {
      cout<<"Computing " << opname(P,N)
          <<" restricted to subspace of dimension "<<dim(s)<<" ..."<<flush;
    }
  smat ans = s_heckeop_restricted(P,s,dual,0); // Automatically chooses W or T
  if(v)
    {
      cout<<"done."<<endl;
    }
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
