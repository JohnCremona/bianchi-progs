// FILE HECKE.CC: Implemention of Hecke operators for class homspace

#include <eclib/msubspace.h>
#include <eclib/xmod.h>
#include "homspace.h"

const string W_opname("W");
const string T_opname("T");

string homspace::opname(const Quadprime& P)
{
  ostringstream ans;
  ans << (P.divides(N) ? W_opname : T_opname) << "(" << P << ")";
  return ans.str();
}

string homspace::opname(const Qideal& P)
{
  ostringstream ans;
  ans << (P.divides(N) ? W_opname : T_opname) << "(" << P << ")";
  return ans.str();
}

vector<long> homspace::eigrange(long i)
{
  vector<long> ans;
  if((i<0)||(i>=nap)) return ans;  // shouldn't happen
  Quadprime P = primelist[i];
  long normp = P.norm();
  if (verbose)
    cout << "eigrange for P = " << P << ":\t";
  if(P.divides(N))
    {
      vector<long> ans = {-1, 1};
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
      vector<long> ans(l);
      ans[0]=0;
      for(ap=-aplim; ap<=aplim; ap++)
	ans[ap+aplim] = ap;
      if (verbose)
	cout << ans << endl;
      return ans;
    }
}

vec homspace::applyop(const matop& mlist, const RatQuad& alpha, int proj)
{ vec ans(rk);
  if (proj) ans.init(projcoord.ncols());
  for (vector<mat22>::const_iterator mi = mlist.mats.begin(); mi!=mlist.mats.end(); mi++)
    {
      vec part = chain((*mi)(alpha), proj);
      if(hmod)
        ans.addmodp(part,hmod);
      else
        ans += part;
    }
  if(hmod) ans=reduce_modp(ans,hmod);
  return ans;
}

vec homspace::applyop(const matop& mlist, const modsym& m, int proj)
{ vec ans(rk);
  if (proj) ans.init(projcoord.ncols());
  for (vector<mat22>::const_iterator mi = mlist.mats.begin(); mi!=mlist.mats.end(); mi++)
    {
      mat22 M = *mi;
      vec part = chain(M(m.alpha()), M(m.beta()), proj);
      if(hmod)
        ans.addmodp(part,hmod);
      else
        ans += part;
    }
  if(hmod) ans=reduce_modp(ans,hmod);
  return ans;
}

mat homspace::calcop(const string opname, const matop& mlist, int dual, int display)
{
  mat m(rk,rk);
  for (long j=0; j<rk; j++) if (needed[j])
     { vec colj = applyop(mlist,freemods[j]);
       if(hmod) colj=reduce_modp(colj,hmod);
       m.setcol(j+1,colj);
     }
  if(cuspidal) m = restrict_mat(smat(m),kern).as_mat();
  if(dual) m = transpose(m);
  if (display) cout << "Matrix of " << opname << " = " << m;
  if (display && (dimension>1)) cout << endl;
  return m;
}

vec homspace::calcop_col(const matop& mlist, int j)
{
  vec colj = applyop(mlist,freemods[j-1]);
  if(hmod) colj=reduce_modp(colj,hmod);
  return colj;
}

mat homspace::calcop_cols(const matop& mlist, const vec& jlist)
{
  int i, j, d = dim(jlist);
  mat m(d,rk);
  for (i=1; i<=d; i++)
    {
      j = jlist[i];
      vec colj = applyop(mlist,freemods[j-1]);
      if(hmod) colj=reduce_modp(colj,hmod);
      m.setcol(i,colj);
     }
  return m;
}

smat homspace::s_calcop_cols(const matop& mlist, const vec& jlist)
{
  int i, j, d = dim(jlist);
  smat m(d,rk);
  for (i=1; i<=d; i++)
    {
      j = jlist[i];
      svec colj = applyop(mlist,freemods[j-1]);
      if(hmod) colj.reduce_mod_p(hmod);
      m.setrow(i,colj);
     }
  return m;
}

smat homspace::s_calcop(const string  opname, const matop& mlist, int dual, int display)
{
  smat m(rk,rk);
  for (long j=0; j<rk; j++) if (needed[j])
     { svec colj = applyop(mlist,freemods[j]);
       if(hmod) colj.reduce_mod_p(hmod);
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
      cout << "Matrix of " << opname << " = ";
      if (dimension>1) cout << "\n";
      cout<<m.as_mat();
    }
  return m;
}

mat homspace::calcop_restricted(const string opname, const matop& mlist, const subspace& s, int dual, int display)
{
  long d=dim(s);
  mat m(d,rk);
  for (long j=0; j<d; j++)
     {
       long jj = pivots(s)[j+1]-1;
       vec colj = applyop(mlist,freemods[jj]);
       if(hmod) colj=reduce_modp(colj,hmod);
       m.setrow(j+1,colj);
     }
  if(hmod)
    m = matmulmodp(m,basis(s),hmod);
  else
    m = m*basis(s);
  if(!dual) m=transpose(m); // dual is default for restricted ops
  if (display) cout << "Matrix of " << opname << " = " << m;
  if (display && (dimension>1)) cout << endl;
  return m;
}

smat homspace::s_calcop_restricted(const string opname, const matop& mlist, const ssubspace& s, int dual, int display)
{
  long d=dim(s);
  smat m(d,rk);
  for (long j=1; j<=d; j++)
     {
       long jj = pivots(s)[j];
       svec colj = applyop(mlist,freemods[jj-1]);
       if(hmod) colj.reduce_mod_p(hmod);
       m.setrow(j,colj);
     }
  m = mult_mod_p(m,basis(s),MODULUS);
  m.reduce_mod_p();
  if(!dual) m=transpose(m); // dual is default for restricted ops
  if (display)
    {
      cout << "Matrix of " << opname << " = " << m.as_mat();
      if (dimension>1) cout << endl;
    }
  return m;
}

mat homspace::heckeop(Quadprime& P, int dual, int display)
{
  return calcop(opname(P), matop(P,N), dual, display);
}

vec homspace::heckeop_col(Quadprime& P, int j, int display)
{
  return calcop_col(matop(P,N), j);
}

mat homspace::heckeop_cols(Quadprime& P, const vec& jlist, int display)
{
  return calcop_cols(matop(P,N), jlist);
}

smat homspace::s_heckeop_cols(Quadprime& P, const vec& jlist, int display)
{
  return s_calcop_cols(matop(P,N), jlist);
}

smat homspace::s_heckeop(Quadprime& P, int dual, int display)
{
  return s_calcop(opname(P), matop(P,N), dual, display);
}

mat homspace::heckeop_restricted(Quadprime& P, const subspace& s, int dual, int display)
{
  return calcop_restricted(opname(P), matop(P,N), s, dual, display);
}

smat homspace::s_heckeop_restricted(Quadprime& P, const ssubspace& s, int dual, int display)
{
  return s_calcop_restricted(opname(P), matop(P,N), s, dual, display);
}

mat homspace::wop(Quadprime& Q, int dual, int display)
{
  return calcop(opname(Q), matop(Q,N), dual,display);
}

mat homspace::fricke(int dual, int display)
{
  return calcop(opname(N), matop(N,N), dual,display);
}

mat homspace::opmat(int i, int dual, int verb)
{
  if((i<0)||(i>=nap)) return mat(dimension);  // shouldn't happen
  Quadprime P = primelist[i];
  if(verbose)
    cout<<"Computing " << (P.divides(N) ? W_opname : T_opname) <<"("<<P<<")...";
  return heckeop(P,dual,verb); // Automatically chooses W or T
}

vec homspace::opmat_col(int i, int j, int verb)
{
  if((i<0)||(i>=nap)) return vec(dimension);  // shouldn't happen
  Quadprime P = primelist[i];
  if(verbose)
    cout<<"Computing " << (P.divides(N) ? W_opname : T_opname) <<"("<<P<<")...";
  return heckeop_col(P,j,verb); // Automatically chooses W or T
}

mat homspace::opmat_cols(int i, const vec& jlist, int verb)
{
  if((i<0)||(i>=nap)) return mat(dimension);  // shouldn't happen
  Quadprime P = primelist[i];
  if(verbose)
    cout<<"Computing " << (P.divides(N) ? W_opname : T_opname) <<"("<<P<<")...";
  return heckeop_cols(P,jlist,verb); // Automatically chooses W or T
}

smat homspace::s_opmat_cols(int i, const vec& jlist, int verb)
{
  if((i<0)||(i>=nap)) return smat(dimension);  // shouldn't happen
  Quadprime P = primelist[i];
  if(verbose)
    cout<<"Computing " << (P.divides(N) ? W_opname : T_opname) <<"("<<P<<")...";
  return s_heckeop_cols(P,jlist,verb); // Automatically chooses W or T
}

mat homspace::opmat_restricted(int i, const subspace& s, int dual, int verb)
{
  if((i<0)||(i>=nap)) return mat(dim(s));  // shouldn't happen
  Quadprime P = primelist[i];
  if(verbose)
    cout<<"Computing " << (P.divides(N) ? W_opname : T_opname) <<"("<<P<<")"
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
      cout<<"Computing " << (P.divides(N) ? W_opname : T_opname) <<"("<<P<<")...";
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
      cout<<"Computing " << (P.divides(N) ? W_opname : T_opname) <<"("<<P
	  <<") restricted to subspace of dimension "<<dim(s)<<" ..."<<flush;
      smat ans = s_heckeop_restricted(P,s,dual,0); // Automatically chooses W or T
      cout<<"done."<<endl;
      return ans;
    }
  else return s_heckeop_restricted(P,s,dual,0); // Automatically chooses W or T
}

vec homspace::maninvector(Quadprime& P)
{
  assert (P.is_principal());
  vector<Quad> resmodp=P.residues();
  Quad p = P.gen();
  vec ans = chain(0,p), part;             // =0, but sets the right length.
  vector<Quad>::const_iterator res=resmodp.begin();
  while(res!=resmodp.end())
    {
      if (*res==0) res++;
      part = chain(*res++,p);
      if(hmod)
        ans.addmodp(part,hmod);
      else
        ans += part;
    }
  if(hmod) ans=reduce_modp(ans,hmod);
  return ans;
}

vec homspace::manintwist(const Quad& lambda, const vector<Quad>& res, vector<int> chitable)
{
  vec ans = chain(0,lambda), part;          // =0, but sets the right length.
  vector<int>::const_iterator chi=chitable.begin();
  vector<Quad>::const_iterator r=res.begin();
  while(r!=res.end())
   {
     part = (*chi++)*chain(*r++,lambda);
      if(hmod)
        ans.addmodp(part,hmod);
      else
        ans += part;
   }
  if(hmod) ans=reduce_modp(ans,hmod);
 return ans;
}

#if (0) // methods not used

vec homspace::projmaninvector(const Quadprime& P)    // Will only work after "proj"
{
  assert (P.is_principal());
  vector<Quad> resmodp=P.residues();
  Quad p = P.gen();
  vec ans = projchain(0,p), part;         // =0, but sets the right length.
  vector<Quad>::const_iterator res=resmodp.begin();
  while(res!=resmodp.end())
    {
      part = projchain(*res++,p);
      if(hmod)
        ans.addmodp(part,hmod);
      else
        ans += part;
    }
  if(hmod) ans=reduce_modp(ans,hmod);
  return ans;
}

vec homspace::newhecke(const Quad& p, const Quad& n, const Quad& d)
                                     // Will only work after "proj"
{
  vec ans = projchain(p*n,d), part;
  vector<Quad> resmodp=residues(p);  Quad dp = d*p;
  vector<Quad>::const_iterator res=resmodp.begin();
  while(res!=resmodp.end())
    {
      part = projchain(n+d*(*res++),dp);
      if(hmod)
        ans.addmodp(part,hmod);
      else
        ans += part;
    }
  if(hmod) ans=reduce_modp(ans,hmod);
  return ans;
}

#endif
