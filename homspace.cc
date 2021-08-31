// FILE HOMSPACE.CC: Implemention of class homspace

//#define USE_CRT // if using smats  mod MODULUS, try CRT-ing with another prime
                // NB this is experimental only

#include <eclib/method.h>
#include <eclib/matrix.h>
#include "euclid.h"
#include "cusp.h"
#include "homspace.h"
#include <assert.h>

homspace::homspace(const Qideal& I, int hp, int cuspid, int verb)
  :N(I)
{
  P1 = P1N(N);
  nsymb = P1.size();
  verbose=verb;
  cuspidal=cuspid;
  hmod = 0;
  plusflag=hp;
  nap = 20;
  primelist = ::primelist(N, nap);
  for (vector<Quadprime>::const_iterator Pi = Quadprimes::list.begin();
       Pi != Quadprimes::list.end() && primelist.size()<(unsigned)nap;
       ++Pi)
    {
      Quadprime P = *Pi;
      if (!P.divides(N))
        primelist.push_back(P);
    }

  ER = edge_relations(&P1, hp, verb);
  ngens = ER.get_ngens();

  FR = face_relations(&ER, hp, verb); // fills relmat with the relations and solves
  denom1 = FR.get_denom();
  rk = FR.get_rank();
  hmod = FR.get_hmod();

  make_freemods();

  kernel_delta();

  if(verbose)
    {
      cout << "number of cusps = " << ncusps << endl;
      if (cuspidal)
        cout << "dimension = " << dimension << endl;
      cout<<"denom1 = "<<denom1<<endl;
      cout<<"denom2 = "<<denom2<<endl;
      cout<<"denom3 = "<<denom3<<endl;
      cout << "Finished constructing homspace.\n";
    }
}

void homspace::make_freemods()
{
  if (rk==0) return;

  long i;
  modsym m;

  freegens.resize(rk);
  for (i=0; i<rk; i++)
    freegens[i] = ER.gen(FR.gen(i+1));
  if (verbose)
    {
      cout << "freegens: ";
      for (i=0; i<rk; i++) cout << freegens[i] << " ";
      cout << endl;
      cout << "Freemods:\n";
    }

  for (i=0; i<rk; i++)
    {
      long j = freegens[i];
      long s = j;
      int t = 0;
      if (n_alphas>1)
        {
          pair<long, int> st = ER.symbol_number_and_type(j);
          s = st.first;  // (c:d) symbol number
          t = st.second; // symbol type (negative for singular edges)
        }
      mat22 U = P1.lift_to_SL2(s);
      m = modsym(U, t);
      freemods.push_back(m);
      if (verbose)
        cout<<"--lifting symbol #"<<s<<" to SL2: "<<U
            <<", type "<<t<<" --> "<<m<<"\n"
            <<i<<": "<<m<<endl;
    }

  if (verbose)
    {
      cout<<"Checking that freemods convert back to unit vectors:"<<endl;
    }
  vec ei(rk); // initially 0
  for (i=0; i<rk; i++)
    {
      m = freemods[i];
      if (verbose)
        cout<< m << " --> " << flush;
      //      vec v = chain(m);
      vec v = chain(m.beta()) - chain(m.alpha());
      ei[i+1] = denom1;
      if (verbose)
        cout << v << flush;
      if (v!=ei)
        {
          cerr<<endl;
          if (!verbose) cerr<< m << " --> " << v<<endl;
          cerr<<" *** WRONG, should be "<<ei<<endl;
          exit(1);
        }
      else
        {
          if (verbose) cout << " OK"<<endl;
        }
      ei[i+1] = 0;
    }
}

void homspace::kernel_delta()
{
  if (verbose)
    cout<<"Computing boundary map"<<endl;
  cusplist cusps(N, plusflag);
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

//#define DEBUG_CHAIN

vec homspace::chaincd(const Quad& c, const Quad& d, int type, int proj)
{
  long ind = P1.index(c,d);
  long i= ER.coords(ind, type);
#ifdef DEBUG_CHAIN
  cout<<"Symbol ("<<c<<":"<<d<<") has index "<<ind<<" plus offset "<< ER.offset(type) <<" = "<<ind+ER.offset(type)
       <<", giving coordindex "<<i;
#endif
  long n = (proj? projcoord.ncols(): rk);
  vec ans(n); // initialises to 0
  if (i)
    {
      ans = sign(i) * (proj? projcoord.row(abs(i)) : coords(abs(i)));
#ifdef DEBUG_CHAIN
      cout << ": coordinate vector "<<ans<<endl;
#endif
    }
  return ans;
}

vec homspace::chain(const RatQuad& alpha, const RatQuad& beta, int proj)
// Instead of just  {return chain(beta, proj) - chain(alpha, proj);}
// we apply a version of "Karim's trick"
{
  Quad a(alpha.num()), b(alpha.den()), x, y;
  Quad g = quadbezout(a,b, x,y);
  //  cout<<"gcd("<<a<<","<<b<<") = " << g <<endl;
  assert (g==1);
  mat22 M(b,-a, x,y);    // det(M)=1 and M(alpha) = 0
  assert (M.det()==Quad::one);
  Quad c = N.reduce(x), d = N.reduce(-b);
#ifdef DEBUG_CHAIN
  cout<<"Computing alpha->beta chain {"<<alpha<<","<<beta<<"}\n";
  cout<<"   translated to {0, "<<M(beta)<<"} with c="<<c<<", d="<<d<<"\n";
#endif
  return chain(M(beta), proj, c, d);
}

vec homspace::chain(const Quad& aa, const Quad& bb, int proj, const Quad& cc, const Quad& dd)
{
  Quad e, a(aa), b(bb), c(cc), d(dd), q, f;
  vec ans = chaincd(c,d,0,proj); // this is the path {0,oo} when (c:d)=(0:1) (the default)
  int t=0, u;
#ifdef DEBUG_CHAIN
  //   if (!Quad::is_Euclidean)
  cout<<" INIT (c:d)_0=("<<c<<":"<<d<<")_0 = "<< modsym(lift_to_SL2(N,c,d),0)<<") AT "<< RatQuad(a,b,1) << endl;
#endif
   while (quadnorm(b))
     {
       pseudo_euclidean_step(a,b, t, c,d);
       //c = N.reduce(c); d = N.reduce(d); // reduce modulo the level

       // either t>=0 and we have a standard edge:
       if (t>=0)
         {
           u = alpha_inv[t];
#ifdef DEBUG_CHAIN
           cout<<" STEP (t="<<t<<", t'="<<u<<", (c:d)_t'=("<<c<<":"<<d<<")_"<<u<<" = "<< modsym(lift_to_SL2(N,c,d),u)<<") TO "<<RatQuad(a,b,1) << endl;
#endif
           // Look up this symbol, convert to a vector w.r.t. homology basis
           vec part = chaincd(c, d, u, proj);
           if(hmod)
             {
               ans.addmodp(part,hmod);
               ans = reduce_modp(ans,hmod);
             }
           else
             ans += part;
         }
       else // t<0 means that the last step took us to sigma[|t|] via
            // a translation only.  We will not reach b=0; instead we
            // finish off by subtracting M{sigma[|t|],oo} where M has
            // second row (c,d). [See Lingham's thesis, p.77]
         {
           vec part = - chaincd(c, d, t, proj);
           if(hmod)
             {
               ans.addmodp(part,hmod);
               ans = reduce_modp(ans,hmod);
             }
           else
             ans += part;
#ifdef DEBUG_CHAIN
           cout<<" coordinate vector = "<<ans<<endl;
#endif
           return ans;
         }
     }

   // We get here when b=0, so no singular edge was used
#ifdef DEBUG_CHAIN
   cout<<" coordinate vector = "<<ans<<endl;
#endif
   return ans;
}

vec reduce_modp(const vec& v, const scalar& p)
{
  long i, d=dim(v);
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
  long i, j, nr=m.nrows(), nc=m.ncols();
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

// List of bad primes (dividing N) followed by good primes to length np:
vector<Quadprime> primelist(Qideal& N, int np)
{
  vector<Quadprime> ans = N.factorization().primes();
  vector<Quadprime>::const_iterator Pi = Quadprimes::list.begin();
  while (ans.size()<(unsigned)np)
    {
      Quadprime P = *Pi++;
      if (!P.divides(N))
        ans.push_back(P);
    }
  return ans;
}

// The remaining methods for class homspace are implemented in hecke.cc.
