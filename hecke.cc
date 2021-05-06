// FILE HECKE.CC: Implemention of Hecke operators for class homspace

#include <eclib/msubspace.h>
#include <eclib/xmod.h>
#include "homspace.h"
#include <assert.h>

const string W_opname("W");
const string T_opname("T");

vector<long> homspace::eigrange(long i)  // implementing virtal function in matmaker
{
  vector<long> ans;
  if((i<0)||(i>=nap)) return ans;  // shouldn't happen
  Quad p = primelist[i];
  long normp = quadnorm(p);
  if (verbose) 
    cout << "eigrange for p = " << p << ":\t";
  if(div(p,modulus))
    {
      vector<long> ans(2);
      ans[0]=-1;
      ans[1]=1;
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

vec homspace::applyop(const matop& mlist, const RatQuad& q) const
{ vec ans(rk), part;
  long i=mlist.length();
  while (i--)
    {
      part = chain(mlist[i](q));
      if(hmod)
        ans.addmodp(part,hmod);
      else
        ans += part;
    }
  if(hmod) ans=reduce_modp(ans,hmod);
  return ans;
}

mat homspace::calcop(const string opname, const Quad& p, const matop& mlist, int dual, int display) const
{
  mat m(rk,rk);
  for (long j=0; j<rk; j++) if (needed[j])
     { vec colj = applyop(mlist,freemods[j]);
       if(hmod) colj=reduce_modp(colj,hmod);
       m.setcol(j+1,colj);
     }
  if(cuspidal) m = restrict_mat(smat(m),kern).as_mat();
  if(dual) m = transpose(m);
  if (display) cout << "Matrix of " << opname << "(" << p << ") = " << m;
  if (display && (dimension>1)) cout << endl;
  return m;
}

vec homspace::calcop_col(const string opname, const Quad& p, const matop& mlist, int j, int display) const
{
  vec colj = applyop(mlist,freemods[j-1]);
  if(hmod) colj=reduce_modp(colj,hmod);
  return colj;
}

mat homspace::calcop_cols(const string opname, const Quad& p, const matop& mlist, const vec& jlist, int display) const
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
 
smat homspace::s_calcop_cols(const string opname, const Quad& p, const matop& mlist, const vec& jlist, int display) const
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
 
smat homspace::s_calcop(const string  opname, const Quad& p, const matop& mlist, 
			int dual, int display) const
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
      cout << "Matrix of " << opname << "(" << p << ") = ";
      if (dimension>1) cout << "\n";
      cout<<m.as_mat();
    }
  return m;
}

mat homspace::calcop_restricted(const string opname, const Quad& p, const matop& mlist, const subspace& s, int dual, int display) const
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
  if (display) cout << "Matrix of " << opname << "(" << p << ") = " << m;
  if (display && (dimension>1)) cout << endl;
  return m;
}
 
smat homspace::s_calcop_restricted(const string opname, const Quad& p, const matop& mlist, const ssubspace& s, int dual, int display) const
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
      cout << "Matrix of " << opname << "(" << p << ") = " << m.as_mat();
      if (dimension>1) cout << endl;
    }
  return m;
}
 
mat homspace::heckeop(const Quad& p, int dual, int display) const
{
 matop matlist(p,modulus);
 string name = (div(p,modulus)) ? W_opname : T_opname;
 return calcop(name,p,matlist,dual,display);
}
 
vec homspace::heckeop_col(const Quad& p, int j, int display) const
{
 matop matlist(p,modulus);
 string name = (div(p,modulus)) ? W_opname : T_opname;
 return calcop_col(name,p,matlist,j,display);
}
 
mat homspace::heckeop_cols(const Quad& p, const vec& jlist, int display) const
{
 matop matlist(p,modulus);
 string name = (div(p,modulus)) ? W_opname : T_opname;
 return calcop_cols(name,p,matlist,jlist,display);
}
 
smat homspace::s_heckeop_cols(const Quad& p, const vec& jlist, int display) const
{
 matop matlist(p,modulus);
 string name = (div(p,modulus)) ? W_opname : T_opname;
 return s_calcop_cols(name,p,matlist,jlist,display);
}
 
smat homspace::s_heckeop(const Quad& p, int dual, int display) const
{
 matop matlist(p,modulus);
 string name = (div(p,modulus)) ? W_opname : T_opname;
 return s_calcop(name,p,matlist,dual,display);
}
 
mat homspace::heckeop_restricted(const Quad& p, const subspace& s, int dual, int display) const
{
 matop matlist(p,modulus);
 string name = (div(p,modulus)) ? W_opname : T_opname;
 return calcop_restricted(name,p,matlist,s,dual,display);
}
 
smat homspace::s_heckeop_restricted(const Quad& p, const ssubspace& s, int dual, int display) const
{
 matop matlist(p,modulus);
 string name = (div(p,modulus)) ? W_opname : T_opname;
 return s_calcop_restricted(name,p,matlist,s,dual,display);
}
 
mat homspace::wop(const Quad& q, int dual, int display) const
{
 matop matlist(q,modulus);
 return calcop(W_opname,q,matlist,dual,display);
}
 
mat homspace::fricke(int dual, int display) const
{
 matop frickelist(modulus,modulus);
 return calcop(W_opname,modulus,frickelist,dual,display);
}

mat homspace::opmat(int i, int dual, int verb)
{
  if((i<0)||(i>=nap)) return mat(dimension);  // shouldn't happen
  Quad p = primelist[i];
  if(verbose) 
      cout<<"Computing " << ((div(p,modulus)) ? W_opname : T_opname) <<"("<<p<<")...";
  return heckeop(p,dual,verb); // Automatically chooses W or T
}

vec homspace::opmat_col(int i, int j, int verb)
{
  if((i<0)||(i>=nap)) return vec(dimension);  // shouldn't happen
  Quad p = primelist[i];
  if(verbose) 
      cout<<"Computing " << ((div(p,modulus)) ? W_opname : T_opname) <<"("<<p<<")...";
  return heckeop_col(p,j,verb); // Automatically chooses W or T
}

mat homspace::opmat_cols(int i, const vec& jlist, int verb)
{
  if((i<0)||(i>=nap)) return mat(dimension);  // shouldn't happen
  Quad p = primelist[i];
  if(verbose) 
      cout<<"Computing " << ((div(p,modulus)) ? W_opname : T_opname) <<"("<<p<<")...";
  return heckeop_cols(p,jlist,verb); // Automatically chooses W or T
}

smat homspace::s_opmat_cols(int i, const vec& jlist, int verb)
{
  if((i<0)||(i>=nap)) return smat(dimension);  // shouldn't happen
  Quad p = primelist[i];
  if(verbose)
    cout<<"Computing " << ((div(p,modulus)) ? W_opname : T_opname) <<"("<<p<<")..."<<flush;
  return s_heckeop_cols(p,jlist,verb); // Automatically chooses W or T
}

mat homspace::opmat_restricted(int i, const subspace& s, int dual, int verb)
{
  if((i<0)||(i>=nap)) return mat(dim(s));  // shouldn't happen
  Quad p = primelist[i];
  if(verbose) 
      cout<<"Computing " << ((div(p,modulus)) ? W_opname : T_opname) <<"("<<p
	  <<") restricted to subspace of dimension "<<dim(s)<<" ..."<<flush;
  return heckeop_restricted(p,s,dual,verb); // Automatically chooses W or T
}

smat homspace::s_opmat(int i, int dual, int v)
{
  //  if(i==-1) return s_conj(dual,v);
  if((i<0)||(i>=nap)) 
    {
      return smat(dimension);  // shouldn't happen
    }
  Quad p = primelist[i];
  if(v) 
    {
      cout<<"Computing " << ((div(p,modulus)) ? W_opname : T_opname) <<"("<<p<<")...";
      smat ans = s_heckeop(p,dual,0); // Automatically chooses W or T
      cout<<"done."<<endl;
      return ans;
    }
  else return s_heckeop(p,dual,0); // Automatically chooses W or T
}

smat homspace::s_opmat_restricted(int i, const ssubspace& s, int dual, int v)
{
  if((i<0)||(i>=nap)) 
    {
      return smat(dim(s));  // shouldn't happen
    }
  Quad p = primelist[i];
  if(v) 
    {
      cout<<"Computing " << ((div(p,modulus)) ? W_opname : T_opname) <<"("<<p
	  <<") restricted to subspace of dimension "<<dim(s)<<" ..."<<flush;
      smat ans = s_heckeop_restricted(p,s,dual,0); // Automatically chooses W or T
      cout<<"done."<<endl;
      return ans;
    }
  else return s_heckeop_restricted(p,s,dual,0); // Automatically chooses W or T
}

vec homspace::maninvector(const Quad& p) const
{
  vector<Quad> resmodp=residues(p);
  vec ans = chain(0,p), part;             // =0, but sets the right length.
  vector<Quad>::const_iterator res=resmodp.begin();
  while(res!=resmodp.end())
    {
      part = chain(*res++,p);
      if(hmod)
        ans.addmodp(part,hmod);
      else
        ans += part;
    }
  if(hmod) ans=reduce_modp(ans,hmod);
  return ans;
}

vec homspace::manintwist(const Quad& lambda, const vector<Quad>& res, vector<int> chitable) const
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

vec homspace::projmaninvector(const Quad& p) const    // Will only work after "proj"
{
  vector<Quad> resmodp=residues(p);
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

vec homspace::newhecke(const Quad& p, const Quad& n, const Quad& d) const
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
