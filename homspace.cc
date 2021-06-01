// FILE HOMSPACE.CC: Implemention of class homspace

//#define USE_CRT // if using smats  mod MODULUS, try CRT-ing with another prime
                // NB this is experimental only

#include <eclib/msubspace.h>
#include <eclib/xmod.h>
#include "homspace.h"
#include "geometry.h"
#include "euclid.h"
#include <assert.h>

homspace::homspace(const Quad& n, int hp, int cuspid, int verb) :symbdata(n)
{
  verbose=verb;
  cuspidal=cuspid;
  hmod = 0;
  if (verbose) symbdata::display();
  plusflag=hp;                  // Sets static level::plusflag = hp
  ngens=0;
  coordindex.resize(nsymbx);
  gens.resize(nsymbx+1);
  //NB start of gens array is at 1 not 0

  edge_relations(); // sets coordindex

  if (verbose)
    {
      cout << "After 2-term relations, ngens = "<<ngens<<endl;
      cout << "gens = ";
      int i, j;
      for (i=1; i<=ngens; i++) cout << gens[i] << " ";
      cout << endl;
      cout << "coordindex = \n";
      for (j=0; j<n_alphas; j++)
        {
          int offset = j*nsymb;
          if(n_alphas>1)
            {
              mat22 M = M_alphas[j];
              RatQuad alpha(-M.entry(1,1), M.entry(1,0));
              cout << "Type " << j << ", alpha = "<< alpha<<", M_alpha = "<<M<<":\n";
            }
          for (i=0; i<nsymb; i++)
            cout << i<<":\t"<<symbol(i)<<"\t"<<coordindex[i+offset] << "\n";
        }
      cout << endl;
    }

  face_relations(); // fills relmat with the relations

  if(verbose)
    {
      cout << "Finished face relations: ";
      cout << "number of relations = " << numrel;
      cout << " (bound was "<<maxnumrel<<")"<<endl;
    }

  solve_relations();

  if (verbose)
    {
      cout << "dimension (relative to cusps) = " << rk << endl;
    }

  make_freemods();

  kernel_delta();

  if(verbose)
    {
      cout << "number of cusps = " << ncusps << endl;
      if (cuspidal)
        cout << "dimension = " << dimension << endl;
    }

  if (verbose) cout << "Finished constructing homspace.\n";
}

#ifdef USE_CRT
int liftmats_chinese(const smat& m1, scalar pr1, const smat& m2, scalar pr2, smat& m, scalar& dd);
#endif

void homspace::solve_relations()
{
   vec pivs, npivs;
   int i;
   if(verbose)
     {
       mat M = relmat.as_mat().slice(numrel,ngens);
       cout<<"relmat = "<<M<<endl;
       cout<<"rank(relmat) = "<<relmat.rank()<<", ngens = "<<ngens<<endl;
     }
#if(USE_SMATS)
   smat_elim sme(relmat);
   int d1;
   smat ker = sme.kernel(npivs,pivs), sp;
   int ok = liftmat(ker,MODULUS,sp,d1);
   if (!ok)
     {
       if(verbose)
         cout << "failed to lift modular kernel using modulus "
              << MODULUS << endl;
#ifdef USE_CRT
       int mod2 = 1073741783; // 2^30-41
       bigint mmod = to_ZZ(MODULUS)*to_ZZ(mod2);
       if(verbose)
         cout << "repeating kernel computation, modulo " << mod2 << endl;
       smat_elim sme2(relmat,mod2);
       vec pivs2, npivs2;
       smat ker2 = sme2.kernel(npivs2,pivs2), sp2;
       ok = (pivs==pivs2);
       if (!ok)
         {
           cout<<"pivs do not agree:\npivs  = "<<pivs<<"\npivs2 = "<<pivs2<<endl;
         }
       else
         {
           if(verbose) cout << " pivs agree" << endl;
           ok = liftmats_chinese(ker,MODULUS,ker2,mod2,sp,d1);
         }
       if (ok)
         {
           if(verbose)
             cout << "success using CRT, combined modulus = "<<mmod
                  <<", denominator= " << d1 << "\n";
         }
       else
         {
           if(verbose)
             cout << "CRT combination with combined modulus "<<mmod<<" failed\n" << endl;
         }
#endif
     }
   if (ok)
     {
       denom1 = d1;
     }
   else
     {
       hmod = MODULUS;
       denom1 = 1;
     }
   relmat=smat(0,0); // clear space
   if(verbose>1)
     {
       cout << "kernel of relmat = " << sp.as_mat() << endl;
       cout << "pivots = "<<pivs << endl;
       cout << "denom = "<<d1 << endl;
     }
   rk = sp.ncols();
   coord.init(ngens+1,rk); // 0'th is unused
   for(i=1; i<=ngens; i++)
     coord.setrow(i,sp.row(i).as_vec());
   // if hmod>0, coord is only defined modulo hmod
   sp=smat(0,0); // clear space
#else
  relmat = relmat.slice(numrel,ngens);
  if(verbose) 
    {
      cout << "relmat = "; relmat.output_pari(cout); cout << endl;
    }
  msubspace sp = kernel(mat_m(relmat),0);
  rk = dim(sp);
  if(verbose>2) cout<<"coord = "<<basis(sp)<<endl;
  coord = basis(sp).shorten((int)1);
  pivs = pivots(sp);
  denom1 = I2int(denom(sp));
  relmat.init(); sp.clear(); 
#endif
  if (verbose)
    {
      cout << "rk = " << rk << endl;
    }
  if (rk>0)
    {
      if (verbose)
        {
          cout << "coord:" << coord;
          if (hmod)
            cout << "failed to lift, coord is only defined modulo "<<hmod<<endl;
          else
            cout << "lifted ok, denominator = " << denom1 << endl;
          cout << "pivots = " << pivs <<endl;
        }
      freegens.resize(rk);
      for (i=0; i<rk; i++) freegens[i] = gens[pivs[i+1]];
      if (verbose)
        {
          cout << "freegens: ";
          for (i=0; i<rk; i++) cout << freegens[i] << " ";
          cout << endl;
        }
    }
}

void homspace::kernel_delta()
{
  if (verbose)
    cout<<"Computing boundary map"<<endl;
  cusplist cusps(2*rk,this);
  mat deltamat(2*rk,rk);
  int i;
  for (i=0; i<rk; i++)
    {
      modsym m = freemods[i];
      deltamat(cusps.index(m.beta())+1, i+1) += 1;  // N.B. offset of 1
      deltamat(cusps.index(m.alpha())+1, i+1) -= 1;
    }
  ncusps=cusps.count();

  kern = kernel(smat(deltamat));
  vec pivs, npivs;
  int d2;
  smat sk;
  int ok = liftmat(smat_elim(deltamat).kernel(npivs,pivs),MODULUS,sk,d2);
  if (!ok)
    cout << "**!!!** failed to lift modular kernel\n" << endl;

  tkernbas = transpose(kern.bas());         // dim(kern) x rank
  if(verbose>1)
    cout<<"tkernbas = "<<tkernbas.as_mat()<<endl;

  dimension = (cuspidal? dim(kern): rk);
  denom2 = d2;
  denom3 = denom1 * denom2;

  const smat& basiskern = basis(kern);
  if (verbose)
    {
      cout << "Basis of ker(delta):\n";
      cout << basiskern.as_mat();
      cout << "pivots: " << pivots(kern) << endl;
    }
  for (i=0; i<rk; i++)
    {
      int n = (cuspidal? ! trivial(basiskern.row(i+1).as_vec()) : 1);
      needed.push_back(n);
      if (verbose)
        {
          cout << "generator "<< i << ": " << freemods[i];
          if (!n) cout << " (not needed)";
          cout << endl;
        }
    }
}

void homspace::make_freemods()
{
  if (verbose)  cout << "Freemods:\n";
  int i,j,s,t=0;
  modsym m;
  for (i=0; i<rk; i++)
    {
      s = j = freegens[i];
      if (n_alphas>1)
        {
          //      cout<<"j = "<<j<<": ";
          std::div_t st = div(j, nsymb);
          s = st.rem;  // remainder gives (c:d) symbol number
          t = st.quot; // quotient gives symbol type
          //      cout<<"(s,t) = ("<<s<<","<<t<<")\n";
        }
      m = modsym(symbol(s), t);
      freemods.push_back(m);
      if (verbose) cout<<m<<endl;
    }
  if (!verbose) return;
  vec ei(rk);
  for (i=0; i<rk; i++)
    {
      m = freemods[i];
      vec v = chain(m.beta()) - chain(m.alpha());
      cout<< m << " --> " << v;
      ei[i+1] = denom1;
      if (v!=ei)
        cout<<" *** WRONG, should be "<<ei;
      cout<<endl;
      ei[i+1] = 0;
    }
}

vec homspace::chaincd(const Quad& c, const Quad& d, int type, int proj) const
{
 long i= coordindex[index2(c,d) + nsymb*type];
 const mat& co = (proj? projcoord: coord);
 vec ans(co.ncols());
 if (i) ans = sign(i)*(co.row(abs(i)));
 return ans;
}

//#define DEBUG_NON_EUCLID

vec homspace::chain(const RatQuad& alpha, const RatQuad& beta, int proj) const
// Instead of just
//  {return chain(beta, proj) - chain(alpha, proj);}
// we apply "Karim's trick"
{
  Quad a(num(alpha)), b(den(alpha)), x, y;
  quadbezout(a,b, x,y); // discard its value which is 1
  mat22 M(b,-a, x,y);    // det(M)=1 and M(alpha) = 0
  assert (M.det()==Quad::one);
#ifdef DEBUG_NON_EUCLID
  cout<<"Computing alpha->beta chain {"<<alpha<<","<<beta<<"}\n";
#endif
  return chain(M(beta), proj, x,-b);
}

vec homspace::chain(const Quad& aa, const Quad& bb, int proj, const Quad& cc, const Quad& dd) const
{
  Quad e, a(aa), b(bb), c(cc), d(dd), q, f;
  vec ans = chaincd(c,d,0,proj); // this is the path {0,oo} when (c:d)=(0:1) (the default)
  int t=0, u;
#ifdef DEBUG_NON_EUCLID
  //   if (!Quad::is_Euclidean)
  cout<<" INIT (c:d)_0=("<<c<<":"<<d<<")_0 = "<< modsym(symb(c,d,this),0)<<") AT "<< RatQuad(a,b,1) << endl;
#endif
   while (quadnorm(b))
     {
       pseudo_euclidean_step(a,b, t, c,d);
       assert (t!=-1);
       u = alpha_inv[t];
#ifdef DEBUG_NON_EUCLID
       //       if (!Quad::is_Euclidean)
         cout<<" STEP (t="<<t<<", t'="<<u<<", (c:d)_u=("<<c<<":"<<d<<")_"<<u<<" = "<< modsym(symb(c,d,this),u)<<") TO "<<RatQuad(a,b,1) << endl;
#endif
       // Look up this symbol, convert to a vector w.r.t. homology basis
       vec part = chaincd(c, d, u, proj);
       if(hmod)
         ans.addmodp(part,hmod);
       else
         ans += part;
     }
   if(hmod) ans=reduce_modp(ans,hmod);
#ifdef DEBUG_NON_EUCLID
       //       if (!Quad::is_Euclidean)
   cout<<" coordinate vector = "<<ans<<endl;
#endif
   return ans;
}



// The remaining methods for class homspace are implelemted in
// edge_relations.cc, face_relations.cc and hecke.cc.

// Below are implementations of linear algebra utility functions
// reduce_modp (for vec and mat) and liftmats_chinese.


#ifdef USE_CRT

//#define DEBUG_CHINESE

int liftmats_chinese(const smat& m1, scalar pr1, const smat& m2, scalar pr2,
                     smat& m, scalar& dd)
{
  long modulus=(long)pr1*(long)pr2,n,d,mij;
  long nr,nc,u,v;
  float lim=floor(sqrt(modulus/2.0));

  dd = bezout(pr1,pr2,u,v); //==1
  if (dd!=1) return 0;

  // First time through: compute CRTs, common denominator and success flag
  m = m1; // NB We assume that m1 and m2 have nonzero entries in the same places
  for(nr=0; nr<m1.nro; nr++)
    for(nc=0; nc<m1.col[nr][0]; nc++)
      {
        mij = mod(v*m1.val[nr][nc],pr1)*pr2 + mod(u*m2.val[nr][nc],pr2)*pr1;
        mij = mod(mij,modulus);
#ifdef DEBUG_CHINESE
        if (((mij-m1.val[nr][nc])%pr1)||((mij-m2.val[nr][nc])%pr2))
          {
            cout<< "bad CRT(["<<m1.val[nr][nc]<<","<<m2.val[nr][nc]<<"],["<<pr1<<","<<pr2<<"]) = "<<mij<<endl;
          }
#endif
        m.val[nr][nc] = mij;
	if (modrat(mij,modulus,lim,n,d))
          dd=lcm(d,dd);
        else
          {
#ifdef DEBUG_CHINESE
            cout<<"CRT("<<m1.val[nr][nc]<<","<<m2.val[nr][nc]<<")="<<mij<<" (mod "<<modulus<<") fails to lift (lim="<<lim<<")\n";
            cout << "Problems encountered in chinese lifting of smat modulo "<<pr1<<" and "<<pr2<< endl;
#endif
            return 0;
          }
      }
  dd=abs(dd);
#ifdef DEBUG_CHINESE
  cout << "Common denominator = " << dd << "\n";
#endif
  // Second time through: rescale
  for(nr=0; nr<m.nro; nr++)
    for(nc=0; nc<m.col[nr][0]; nc++)
      {
        m.val[nr][nc] = mod(xmodmul((dd/d),(long)m.val[nr][nc],modulus),modulus);
      }
  return 1;
}

#endif

vec reduce_modp(const vec& v, const scalar& p)
{
  int i, d=dim(v);
  scalar ai, p2 = p>>1;
  vec ans(d);
  for(i=1; i<=dim(v); i++)
    {
      ai = v[i]%p;
      while( ai>p2) ai-=p;
      while(-ai>p2) ai+=p;
      ans[i] = ai;
    }
  return ans;
}

mat reduce_modp(const mat& m, const scalar& p)
{
  int i, j, nr=m.nrows(), nc=m.ncols();
  scalar aij, p2 = p>>1;
  mat ans(nr,nc);
  for(i=1; i<=nr; i++)
    for(j=1; j<=nc; j++)
      {
        aij = m(i,j)%p;
        while( aij>p2) aij-=p;
        while(-aij>p2) aij+=p;
        ans(i,j) = aij;
      }
  return ans;
}
