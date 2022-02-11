// FILE HECKE.CC: Implemention of Hecke operators for class homspace

#include <eclib/msubspace.h>
#include <eclib/xmod.h>
#include "homspace.h"

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
  if(display)
    cout<<"Computing " << T.name() <<"...";
  mat m(rk,rk);
  for (long j=0; j<rk; j++) if (needed[j])
     { vec colj = applyop(T,freemods[j]);
       m.setcol(j+1,colj);
     }
  if(cuspidal) m = restrict_mat(smat(m),kern).as_mat();
  if(dual) m = transpose(m);
  // if (display) cout << "Matrix of " << T.name() << " = " << m;
  // if (display && (dimension>1)) cout << endl;
  if (display)
    cout<<"done."<<endl;
  return m;
}

mat homspace::calcop_cols(const matop& T, const vec& jlist, int verb)
{
  if(verb)
    cout<<"Computing " << T.name() <<"...";
  int i, d = dim(jlist);
  mat m(d,rk);
  for (i=1; i<=d; i++)
    {
      int j = jlist[i];
      vec colj = applyop(T,freemods[j-1]);
      m.setcol(i,colj);
     }
  if (verb)
    cout<<"done."<<endl;
  return m;
}

smat homspace::s_calcop_cols(const matop& T, const vec& jlist, int verb)
{
  if(verb)
    cout<<"Computing " << T.name() <<"...";
  int i, d = dim(jlist);
  smat m(d,rk);
  for (i=1; i<=d; i++)
    {
      int j = jlist[i];
      svec colj = applyop(T,freemods[j-1]);
      m.setrow(i,colj);
     }
  if (verb)
    cout<<"done."<<endl;
  return m;
}

smat homspace::s_calcop(const matop& T, int dual, int display)
{
  if(display)
    cout<<"Computing " << T.name() <<"..."<<flush;
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
  // if (display)
  //   {
  //     cout << "Matrix of " << T.name() << " = ";
  //     if (dimension>1) cout << "\n";
  //     cout<<m.as_mat();
  //   }
  if(display)
    cout<<"done."<<endl;
  return m;
}

mat homspace::calcop_restricted(const matop& T, const subspace& s, int dual, int display)
{
  if(display)
    cout<<"Computing " << T.name()
        <<" restricted to subspace of dimension "<<dim(s)<<" ..."<<flush;
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
  // if (display) cout << "Matrix of " << T.name() << " = " << m;
  // if (display && (dimension>1)) cout << endl;
  if (display)
    cout<<"done."<<endl;
  return m;
}

smat homspace::s_calcop_restricted(const matop& T, const ssubspace& s, int dual, int display)
{
  if(display)
    cout<<"Computing " << T.name() <<"..."
        <<" restricted to subspace of dimension "<<dim(s)<<" ..."<<flush;
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
  // if (display)
  //   {
  //     cout << "Matrix of " << T.name() << " = " << m.as_mat();
  //     if (dimension>1) cout << endl;
  //   }
  if(display)
    cout<<"done."<<endl;
  return m;
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
