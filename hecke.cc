// FILE HECKE.CC: Implemention of Hecke operators for class homspace

#include <eclib/msubspace.h>
#include <eclib/xmod.h>
#include "homspace.h"

vector<long> homspace::eigrange(long i)
{
  vector<long> ans;
  if((i<0)||(i>=nap)) return ans;  // shouldn't happen
  Quadprime P = primelist[i];
  long normp = I2long(P.norm());
  if (verbose)
    cout << "eigrange for P = " << P << ":\t";
  if(P.divides(N))
    {
      ans = {-1, 1};
      if (verbose)
	cout << ans << endl;
      return ans;
    }
  else
    {
      long aplim=2;
      while (aplim*aplim<=4*normp) aplim++;
      aplim--;
      if(verbose)
	cout << "|ap| up to "<<aplim<<":\t";
      long ap, l = 2*aplim+1;
      ans = vector<long>(l);
      ans[0]=0;
      for(ap=-aplim; ap<=aplim; ap++)
	ans[ap+aplim] = ap;
      if (verbose)
	cout << ans << endl;
      return ans;
    }
}

vec homspace::applyop(const matop& T, const RatQuad& alpha, int proj)
{ vec ans(rk);
  if (proj) ans.init(projcoord.ncols());
  for (vector<mat22>::const_iterator mi = T.mats.begin(); mi!=T.mats.end(); ++mi)
    ans = reduce_modp(ans + chain((*mi)(alpha), proj), hmod);
  return ans;
}

//#define DEBUG_APPLYOP
vec homspace::applyop(const matop& T, const modsym& m, int proj)
{ vec ans(rk);
  if (proj) ans.init(projcoord.ncols());
#ifdef DEBUG_APPLYOP
  cout<<"In applyop() with modular symbol "<<m<<" (proj = "<<proj<<")"<<endl;
#endif
  for (vector<mat22>::const_iterator mi = T.mats.begin(); mi!=T.mats.end(); ++mi)
    {
      mat22 M = *mi;
#ifdef DEBUG_APPLYOP
      cout<<"Applying matrix "<<M<<"..."<<flush;
#endif
      vec part = chain(M(m.alpha()), M(m.beta()), proj);
#ifdef DEBUG_APPLYOP
      cout<<" summand = "<<part<<"..."<<flush;
#endif
      ans = reduce_modp(ans + part, hmod);
#ifdef DEBUG_APPLYOP
      cout<<" partial sum = "<<ans<<"..."<<endl;
#endif
    }
#ifdef DEBUG_APPLYOP
  cout<<" final sum = "<<ans<<"..."<<endl;
#endif
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

vec homspace::calcop_col(const matop& T, int j)
{
  vec colj = applyop(T,freemods[j-1]);
  return colj;
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

mat homspace::heckeop(Quadprime& P, int dual, int display)
{
  return calcop(AtkinLehnerOrHeckeOp(P,N), dual, display);
}

vec homspace::heckeop_col(Quadprime& P, int j, int display)
{
  return calcop_col(AtkinLehnerOrHeckeOp(P,N), j);
}

mat homspace::heckeop_cols(Quadprime& P, const vec& jlist, int display)
{
  return calcop_cols(AtkinLehnerOrHeckeOp(P,N), jlist);
}

smat homspace::s_heckeop_cols(Quadprime& P, const vec& jlist, int display)
{
  return s_calcop_cols(AtkinLehnerOrHeckeOp(P,N), jlist);
}

smat homspace::s_heckeop(Quadprime& P, int dual, int display)
{
  return s_calcop(AtkinLehnerOrHeckeOp(P,N), dual, display);
}

mat homspace::heckeop_restricted(Quadprime& P, const subspace& s, int dual, int display)
{
  return calcop_restricted(AtkinLehnerOrHeckeOp(P,N), s, dual, display);
}

smat homspace::s_heckeop_restricted(Quadprime& P, const ssubspace& s, int dual, int display)
{
  return s_calcop_restricted(AtkinLehnerOrHeckeOp(P,N), s, dual, display);
}

mat homspace::wop(Quadprime& Q, int dual, int display)
{
  return calcop(AtkinLehnerOp(Q,N), dual,display);
}

mat homspace::fricke(int dual, int display)
{
  return calcop(FrickeOp(N), dual,display);
}

mat homspace::opmat(int i, int dual, int verb)
{
  if((i<0)||(i>=nap)) return mat(dimension);  // shouldn't happen
  Quadprime P = primelist[i];
  if(verbose)
    cout<<"Computing " << opname(P,N) <<"...";
  return heckeop(P,dual,verb); // Automatically chooses W or T
}

vec homspace::opmat_col(int i, int j, int verb)
{
  if((i<0)||(i>=nap)) return vec(dimension);  // shouldn't happen
  Quadprime P = primelist[i];
  if(verbose)
    cout<<"Computing " << opname(P,N) <<"...";
  return heckeop_col(P,j,verb); // Automatically chooses W or T
}

mat homspace::opmat_cols(int i, const vec& jlist, int verb)
{
  if((i<0)||(i>=nap)) return mat(dimension);  // shouldn't happen
  Quadprime P = primelist[i];
  if(verbose)
    cout<<"Computing " << opname(P,N) <<"...";
  return heckeop_cols(P,jlist,verb); // Automatically chooses W or T
}

smat homspace::s_opmat_cols(int i, const vec& jlist, int verb)
{
  if((i<0)||(i>=nap)) return smat(dimension);  // shouldn't happen
  Quadprime P = primelist[i];
  if(verbose)
    cout<<"Computing " << opname(P,N) <<"...";
  return s_heckeop_cols(P,jlist,verb); // Automatically chooses W or T
}

mat homspace::opmat_restricted(int i, const subspace& s, int dual, int verb)
{
  if((i<0)||(i>=nap)) return mat(dim(s));  // shouldn't happen
  Quadprime P = primelist[i];
  if(verbose)
    cout<<"Computing " << opname(P,N)
        <<" restricted to subspace of dimension "<<dim(s)<<" ..."<<flush;
  return heckeop_restricted(P,s,dual,verb); // Automatically chooses W or T
}

smat homspace::s_opmat(int i, int dual, int v)
{
  //  if(i==-1) return s_conj(dual,v);
  if((i<0)||(i>=nap))
    {
      return smat(dimension);  // shouldn't happen
    }
  Quadprime P = primelist[i];
  if(v)
    {
      cout<<"Computing " << opname(P,N) <<"...";
      smat ans = s_heckeop(P,dual,0); // Automatically chooses W or T
      cout<<"done."<<endl;
      return ans;
    }
  else return s_heckeop(P,dual,0); // Automatically chooses W or T
}

smat homspace::s_opmat_restricted(int i, const ssubspace& s, int dual, int v)
{
  if((i<0)||(i>=nap))
    {
      return smat(dim(s));  // shouldn't happen
    }
  Quadprime P = primelist[i];
  if(v)
    {
      cout<<"Computing " << opname(P,N)
          <<" restricted to subspace of dimension "<<dim(s)<<" ..."<<flush;
      smat ans = s_heckeop_restricted(P,s,dual,0); // Automatically chooses W or T
      cout<<"done."<<endl;
      return ans;
    }
  else return s_heckeop_restricted(P,s,dual,0); // Automatically chooses W or T
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
