//  newforms.cc

#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <functional>   // std::multiplies
#include <numeric>   // std::multiplies
#include "looper.h"
#include "newforms.h"
#include "eclib/curvesort.h" // for letter codes

// Notes on scaling:
//
// For a newform F (however normalised), the periods of F are the
// integrals of (the differential associated to) F along paths between
// Gamma_0(N)-equivalent cusps.  These are all integral multiples of a
// minimal positive real period P.  The map from integral homology to
// the associated period multiple takes gamma In H_1(-,Z) to
// n_F(gamma) = (integral of F over gamma)/P, and the image of this
// map is (by definiton) Z.

// Similarly the "cuspidal periods" of F are its integrals along all
// paths gamma between cusps, whether or not they are
// Gamma_0(N)-equivalent, and these are the integral multiples of a
// least positive real cuspidal period P'.  Now P*Z is a subgroup of
// P'*Z, and its index c is called the "cuspidal factor", also equal
// to c=P/P'.  Let n'_F(gamma) = (integral of F over gamma)/P' for gamma
// in H_1(-,Z; cusps).

// The map from gamma to n'_F(gamma) is modular and an eigenfunction
// for Hecke, so is a primitive basis vector for the associated dual
// eigenspace.  Since we do linear algebra on the dual (by transposing
// Hecke matrices), we can compute this map, up to scaling, by simply
// taking the dot product of our dual basis vector v with the coords
// of gamma with respect to the homology basis.  Regardless of whether
// our "freegens" generator the whole integral homology (w.r.t. cusps)
// or only a sublattice of finite index, i.e. whether denom1=1 or >1,
// this does not nmatter since the map n'_F is primitive, so to get
// the correct values (up to sign) we just have to divide out these
// dot products by their gcd.

// The map n'_F(gamma) is almost encoded in the matrix projcoord
// (created by calling make_projcoord()): its rows are indexed by
// edges modulo edge relations, so to get the value for the j'th
// newform on an edge (c:d)_t we find the edge number i =
// offset(t)+index(c,d)), use that to look up
// k=ER.coords(i)=ER.coords(index(c,d), t), and the value is
// sign(k)*projcoord[abs(k),j].

// For example n'_F({0,oo}) is obtained this way with (c,d,t)=(0,1,0).
// Now L(F,1) is the integral of F over {0,oo} (?times a normalising
// factor?), hence we obtain L(F,1)/P' = n'_F({0,oo}) (an integer).
// For L(F,1)/P we divide by the cuspidal factor for F, obtaining a
// rational.  [P=c*P' so L/P = L/(c*P') = (L/P')/c.]

// A second way to compute L(F,1)/P' which only involves integral
// periods is to use a "Manin vector" mvp for some good prime p, which
// is the sum over x mod p of {0,x/p}, since (1+N(p)-a_p)*n'_F({0,oo})
// = n'_F(mvp), hence L/P' = n'_F(mvp)/(1+N(p)-a_p).

// NB in the newform constructor we cannot use projcoord since that is
// only computed after all the newforms are found

newform::newform(newforms* nfs, const vec& v, const vector<long>& eigs)
  : eigs(eigs)
{
  //cout<<"Constructing newform with eigs "<<eigs<<" and basis "<<v<<endl;
  nf=nfs;

  // convert basis vector from coords w.r.t. face-gens to coords
  // w.r.t.edge-gens:
  if(nf->hmod)
    { // we don't have a mod p mat*vec
      mat vcol(dim(v),1);
      vcol.setcol(1,v);
      basis = reduce_modp(matmulmodp(nf->h1->FR.coord, vcol, nf->hmod).col(1));
    }
  else
    basis = (nf->h1->FR.coord)*v;
  if (nf->characteristic==0)
    {
      makeprimitive(basis); // this is now independent of h1's denom1
    }
  //  else we have just reduced it mod p
  if (nf->verbose >1)
  {
    cout << "short newform basis = "<<v<<endl;
    cout << "long  newform basis = "<<basis<<endl;
  }
  // Find the ratio of the least period w.r.t. integral homology
  // divided by the least period w.r.t. homology relative to cusps.
  // this uses the vector v so must be done now
  find_cuspidal_factor(v);
}

// fill in data for the j'th newform (j based at 1)

void newform::fill_in_data(int j)
{
  // extract Atkin-Lehner eigs from first entries of eigs list:
  copy(eigs.begin(), eigs.begin()+(nf->npdivs), back_inserter(aqlist));

  // Sign of functional equation = minus product of all A-L eigenvalues
  sfe = std::accumulate(aqlist.begin(),aqlist.end(),-1,std::multiplies<long>());

  // compute L/P as n_F({0,oo})
  int pdot0 = abs(nf->zero_infinity[j]);
  loverp =  rational(pdot0, (Quad::nunits) * cuspidalfactor);

  // compute L/P again using Manin vector
  dp0  =  1 + (nf->nP0) - nf->aP0[j-1];  // aP0 is based at 0
  pdot = abs(nf->mvp[j]);
  rational loverp_mvp(pdot, dp0 * (Quad::nunits) * cuspidalfactor);

  // Check they agree:
  if (nf->characteristic>0)
    return;

  if (pdot != dp0*pdot0)
    {
      cout << "Inconsistent values for L/P computed two ways!"<<endl;
      cout << "from {0,oo} directly: " << loverp <<endl;
      cout << "pdot0 = "<<pdot0<<endl;
      cout << "cuspidalfactor = "<<cuspidalfactor<<endl;
      cout << "from Manin vector:    " << loverp_mvp <<endl;
      cout << "pdot = "<<pdot<<endl;
      cout << "nP0 = "<<nf->nP0<<endl;
      cout << "iP0 = "<<nf->iP0<<endl;
      cout << "eigs (size "<<eigs.size()<<") = "<<eigs<<endl;
      cout << "ap0 = "<<nf->aP0[j]<<endl;
      cout << "dp0 = "<<dp0<<endl;
      cout << "cuspidalfactor = "<<cuspidalfactor<<endl;
    }

  // find (a,b,c,d) such that cusp b/d is equivalent to 0 and the
  // integral over {0,M(0)} = {0,b/d} with M = [a,b;N*c,d] is a
  // nonzero multiple "matdot" of the period P.

  // NB We do not currently use this, it should be further scaled
  find_matrix(j);
}

newform::newform(newforms* nfs,
                 const vector<int>& intdata, const vector<Quad>& Quaddata,
                 const vector<long>& aq, const vector<long>& ap)
{
  nf=nfs;
  sfe = intdata[0];
  //  cout<<"sfe from file = "<<sfe<<endl;
  pdot = intdata[1];
  dp0 = intdata[2];
  loverp = rational(abs(pdot),dp0*(Quad::nunits));
  cuspidalfactor = intdata[3];
  lambda = Quaddata[0];
  lambdadot = intdata[4];
  a = Quaddata[1];
  b = Quaddata[2];
  c = Quaddata[3];
  d = Quaddata[4];
  matdot = intdata[5];
  aqlist = aq;
  // Recompute sign of functional equation = minus product of all A-L eigenvalues
  int newsfe = std::accumulate(aqlist.begin(),aqlist.end(),-1,std::multiplies<long>());
  if (newsfe!=sfe)
    cout<<"Problem in data on file for level "<<nfs->N<<": sfe = "<<sfe<<" and aqlist = "<<aqlist<<", but minus product of latter is "<<newsfe<<endl;
  aplist = ap;
  // recreate eigs list (in case we need to recover basis vector):
  // start with Atkin-Lehner eigs aq, then Tp eigs ap for good p
  eigs = aq;
  vector<Quadprime>::const_iterator pr=Quadprimes::list.begin();
  vector<long>::const_iterator api=aplist.begin();
  Qideal N(nf->N);
  while (((int)eigs.size() < nf->nap) && (api!=aplist.end()))
    {
      while (pr->divides(N))  { ++pr; ++api; }
      eigs.push_back(*api);
      ++pr; ++api;
    }
}

void newform::find_cuspidal_factor(const vec& v)
{
  if(nf->h1->cuspidal || nf->characteristic)
    {
      cuspidalfactor=1;
    }
  else
    {
      cuspidalfactor = vecgcd((nf->h1->tkernbas)*v);
      if(nf->verbose)
        cout<<"cuspidalfactor = "<<cuspidalfactor<<endl;
    }
}

// find (a,b,c,d) such that cusp b/d is equivalent to 0 and the
// integral over {0,M(0)} = {0,b/d} with M = [a,b;N*c,d] is a
// nonzero multiple "matdot" of the period P.

//  for the j'th newform (j based at 1)

void newform::find_matrix(int j)
{
  if(nf->verbose>1)
    cout<<"computing integration matrix..."<<flush;
  matdot=0;
  Qideal N(nf->N);
  for (Quadlooper dl(2, 1000, 1); dl.ok()&&!matdot; ++dl)
    { d=(Quad)dl;
      Qideal D(d);
      if (N.is_coprime_to(D))
        {
          vector<Quad> reslist = residues(d);
          vector<Quad>::const_iterator res;
          for(res=reslist.begin(); res!=reslist.end() && !matdot; ++res)
            {
              b=*res;
              Qideal bN = b*N;
              if (D.is_coprime_to(bN, a, c))
                // found a candidate q=b/d: a+c=1 with d|a and b|c and c/b in N
                {
                  c /= -b;
                  a /= d; // now a*d-b*c=1 with c in N
                  assert (a*d-b*c==Quad::one);
                  matdot = abs((nf->h1->chain(b,d, 1))[j]);
                } // b coprime to d test
            } // loop over b
        } // d coprime to N test
    } // loop over d
  if(nf->verbose>1)
    cout<<"M = ["<<a<<","<<b<<";"<<c<<","<<d<<"] with factor "<<matdot<<endl;
}

int newform::is_base_change(void) const
{
  if(!(nf->N.is_Galois_stable()))
    return 0;
  vector<long>::const_iterator ap = aplist.begin();
  vector<Quadprime>::const_iterator pr=Quadprimes::list.begin();
  while(ap!=aplist.end())
    {
      long api = *ap++;
      Quadprime p0 = *pr++;
      //cout<<"p="<<p0<<" has ap="<<api<<endl;
      if(!p0.is_Galois_stable()) // this prime not inert or ramified
        {
          if (ap==aplist.end()) // the conjugate ap is not known
            return 1;
          long apj = *ap++;
          //cout<<"Next prime has ap="<<apj<<endl;
          if(api!=apj) // ap mismatch
            {
              //cout<<"Mismatch -- not base-change"<<endl;
              return 0;
            }
          ++pr;
        }
    }
  //cout<<"All OK -- base-change"<<endl;
  return 1;
}

int newform::is_base_change_twist(void) const
{
  vector<long>::const_iterator ap = aplist.begin();
  vector<Quadprime>::const_iterator pr=Quadprimes::list.begin();
  Qideal N(nf->N);
  while(ap!=aplist.end())
    {
      long api = *ap++;
      Quadprime p0 = *pr++;
      //cout<<"p="<<p0<<" has ap="<<api<<endl;
      if(!p0.is_Galois_stable()) // this prime not inert or ramified
        {
          if (ap==aplist.end()) // the conjugate ap is not known
            {
              //cout<<"All OK -- base-change up to twist"<<endl;
              return 1;
            }
            // read next (conjugate) prime and eigenvalue:
          long apj = *ap++;
          Quadprime P1 =  *pr++;
          // skip if either divides level:
          if((nf->P0).divides(N))
            continue;
          if(P1.divides(N))
            continue;
          // Check the ap agree up to sign:
          if(abs(api)!=abs(apj)) // ap mismatch
            {
              //cout<<"Mismatch -- not base-change-twist"<<endl;
              return 0;
            }
        }
    }
  //cout<<"All OK -- base-change up to twist"<<endl;
  return 1;
}

// if form is base-change, find the d s.t. the bc has eigenvalues in Q(sqrt(d))
int newform::base_change_discriminant(void) const
{
  if (is_base_change()==0) return 0;
  int bcd = 1;
  Qideal N(nf->N);
  vector<long>::const_iterator api = aplist.begin();
  vector<Quadprime>::const_iterator pr=Quadprimes::list.begin();
  while(api!=aplist.end())
    {
      long ap = *api++;
      Quadprime P = *pr++;
      if(!P.is_inert())
        continue;
      if(P.divides(N)) // this prime is bad
        continue;
      long dp = ap+2*P.prime();
      //cout<<"p="<<p<<" has ap="<<ap<<", disc = "<<dp;
      dp = squarefree_part(dp);
      //cout<<" with squarefree part "<<dp<<endl;
      if (dp==0) continue;
      if (bcd==1) // first one
        {
          bcd = dp;
        }
      else
        {
          if (dp!=bcd) // mismatch: not possible?
            {
              //cout<<"mismatch: bcd=0"<<endl;
              return 1;
            }
        }
    }
  return bcd;
}

// if form is twist of base-change, find the d s.t. the bc has eigenvalues in Q(sqrt(d))
int newform::base_change_twist_discriminant(void) const
{
  if (is_base_change_twist()==0) return 0;
  int bcd1, bcd2, cmd = is_CM();
  if (cmd!=0)
    {
      bcd1 = -cmd/(Quad::d);
      //cout << "base change twist and CM("<<cmd<<") so bcd = "<<bcd1<<endl;
      return bcd1;
    }
  bcd1 = bcd2 = 1;
  Qideal N(nf->N);
  vector<long>::const_iterator api = aplist.begin();
  vector<Quadprime>::const_iterator pr=Quadprimes::list.begin();
  while(api!=aplist.end())
    {
      long ap = *api++;
      Quadprime P = *pr++;
      if(!P.is_inert())
        continue;
      if(P.divides(N)) // this prime is bad
        continue;
      long dp1 = ap  + 2*P.prime();
      long dp2 = dp1 - 2*ap;
      //cout<<"p="<<p<<" has ap="<<ap<<", discs "<<dp1<<", "<<dp2;
      dp1 = squarefree_part(dp1);
      dp2 = squarefree_part(dp2);
      //cout<<" with squarefree parts "<<dp1<<", "<<dp2<<endl;
      if (dp1*dp2==0) continue;
      if (dp1==dp2) return dp1; // no ambiguity!
      if ((bcd1==1)&&(bcd2==1))  // first pair, store
        {
          bcd1 = dp1;
          bcd2 = dp2;
        }
      else // see if only one is a repeatl if so it's the value we want
        {
          if ((dp1!=bcd1)&&(dp1!=bcd2))
            {
              return dp2;
            }
          if ((dp2!=bcd1)&&(dp2!=bcd2))
            {
              return dp1;
            }
          // otherwise we have the same two values as before so we go on
        }
    }
  //cout<<"Warning from base_change_twist_discriminant(): unable to eliminate either of "<<bcd1<<" or "<<bcd2<<" so returning 1";
  return 1;
}

// Test if form is CM, return 0 or the CM disc
int newform::is_CM(void) const
{
  int cmd = 0;
  vector<long>::const_iterator api = aplist.begin();
  vector<Quadprime>::const_iterator pr=Quadprimes::list.begin();
  while(api!=aplist.end())
    {
      long ap = *api++;
      Quadprime P = *pr++;
      if (ap==0) continue;
      long dp = ap*ap-4*I2long(P.norm());
      //cout<<"p="<<p<<" has ap="<<ap<<", disc = "<<dp;
      dp = squarefree_part(dp);
      //cout<<" with squarefree part "<<dp<<endl;
      if (dp==0) continue;
      if (cmd==0) // first one
        {
          cmd = dp;
        }
      else
        {
          if (dp!=cmd) // mismatch: not CM
            {
              //cout<<"mismatch: CM=0"<<endl;
              return 0;
            }
        }
    }
  return cmd;
}

void newforms::makeh1plus(void)
{
  if(!h1)
    {
      h1 = new homspace(N,1,0,0, characteristic);
      nfhmod=hmod = h1->h1hmod();
    }
}

long newforms::matdim(void) {return h1->dimension;}
long newforms::matden(void) {return h1->denom3;}

newforms::newforms(const Qideal& iN, int disp, long ch)
  : N(iN), verbose(disp), characteristic(ch)
{
  init();
}

void newforms::init()
{
  is_square = N.is_square();
  Factorization F = N.factorization();
  // List of bad primes (dividing N) followed by good primes to length np:
  nap = 20;
  plist = make_primelist(N, nap, iP0); // shadows h1->primelist in case we do not construct h1
  P0 = plist[iP0];
  nP0 = I2long(P0.norm());
//
// P0 is the smallest good principal prime: and iP0 its index (in
// plist, which starts with the bad primes and then the good
// primes in order).  P0 must be principal since we have only
// implemented maninvector() for principal primes.
//
  if (verbose>1)
    cout << "Ordered list of primes (bad primes first): "<<plist<<endl;
  nwq = npdivs = F.size();
  ntp=nap-nwq;
  h1=0;
  of=0;
  nfhmod=0;
}

//#define DEBUG_LAMBDA
void newforms::find_lambdas()
{
  vector<int> gotlambda(n1ds);
  int i, nfound=0;
#ifdef DEBUG_LAMBDA
  if(verbose) cout<<"Looking for twisting primes.\n";
#endif
  for (i=0; i<n1ds; i++)
    if(nflist[i].pdot!=0)
      {
#ifdef DEBUG_LAMBDA
        if(verbose) cout<<"Newform "<<i<<": lambda=1 will do.\n";
#endif
        nflist[i].lambda=Quad::one;
        nflist[i].lambdadot=nflist[i].pdot;
        gotlambda[i]=1;
        nfound++;
      }
    else
      {
        nflist[i].lambda=Quad::zero; // indicates no lambda exists (yet)
        nflist[i].lambdadot=0;
	gotlambda[i]=0;
      }
#ifdef DEBUG_LAMBDA
  if(verbose)cout<<nfound<<" easy cases out of "<<n1ds<<endl;
#endif
  if (is_square || !N.is_principal())
    {
      return;
    }

  vector<Quadprime>::const_iterator Li;
  for(Li=Quadprimes::list.begin(); Li!=Quadprimes::list.end() && (nfound<n1ds); ++Li)
    {
      Quadprime L = *Li;
      if (L.divides(2)) continue;
      if (L.divides(N)) continue;
      if (!L.is_principal()) continue;
      Quad lam = L.gen();
#ifdef DEBUG_LAMBDA
      if(verbose)cout << "Testing lambda = " << L << " = ("<<lam<<"): principal, odd, good"<<endl;
#endif
      vector<Quad> lamres = L.residues();
      if(squaremod(fundunit,lam,lamres)==1)
        {
#ifdef DEBUG_LAMBDA
          if(verbose)cout<<"passed second test: fundamental unit is a square"<<endl;
#endif
          // TODO:  work out what to do here if N is not principal!
          int chimod  = squaremod(N.gen(),lam,lamres);
          vector<int> chitab = makechitable(lam,lamres);
          vec mvtw = h1->manintwist(lam,lamres,chitab, 1);
          for(int j=0; (j<n1ds)&&(nfound<n1ds); j++)
            {
              if(gotlambda[j]==0)
                {
#ifdef DEBUG_LAMBDA
                  if(verbose)cout<<"Newform # "<<j<<": ";
#endif
                  newform& nfj = nflist[j];
                  int dot = abs(mvtw[j+1]);  // j based at 0 but vec mvtw based at 1
                  if(dot&&((chimod*nfj.sfe)==+1))
                    {
#ifdef DEBUG_LAMBDA
                      if(verbose)cout<<"Success! ";
#endif
                      nfj.loverp = rational(dot, Quad::nunits * nfj.cuspidalfactor);
                      nfj.lambda = lam;
                      nfj.lambdadot = dot;
                      gotlambda[j] = 1;
                      nfound++;
                    }
                }
#ifdef DEBUG_LAMBDA
              if(verbose)cout<<endl;
#endif
            }
        }
    }
}

void newforms::createfromscratch()
{
  if(verbose)
    cout<<"Constructing homspace at level "<<ideal_label(N)<<" ...\n";
  makeh1plus();
  nfhmod=hmod = h1->h1hmod();
  int dimcusp = h1->h1cuspdim();
  int dimall = h1->h1dim();

  if(verbose)
    {
      cout<<"Dimension = "<<dimall<<" (cuspidal dimension = "<<dimcusp<<")\n";
      cout<<"Retrieving oldform data for level "<<N<<" (primelist="<<plist<<")...\n";
    }

  of = new oldforms(N, plist, verbose, characteristic);
  if(verbose)
    {
      of->display();
      cout<<"Finding rational newforms...\n";
    }

  maxdepth = nap;
  long mindepth = npdivs;
  dimsplit = n1ds = 0;
  long olddimall = (of->olddimall);
  nnflist = upperbound = (h1->h1cuspdim()) - olddimall;
  if(verbose)
    {
      cout<<"cuspidal dimension = "<<h1->h1cuspdim()<<", olddimall = "<<olddimall<<", so upper bound = "<<upperbound<<endl;
    }
  if(verbose>1)
    {
      cout<<"upperbound = "<<upperbound<<endl;
      cout<<"maxdepth = "<<maxdepth<<endl;
      cout<<"mindepth = "<<mindepth<<endl;
    }
  if(upperbound<0) // check for error condition
    {
      cout<<"Error:  total old dimension = "<<olddimall<<" as computed is greater than total cuspidal dimension "<<(h1->h1cuspdim())<<" -- aborting"<<endl;
      exit(1);
    }
  if(upperbound>0)  // Else no newforms certainly so do no work!
    {
      use_nf_number=-1; // flags to use() that the nfs found are new
      form_finder ff(this,1,maxdepth,mindepth,1,0,verbose);
      ff.find();
     }
  if(verbose>1) cout<<"n1ds = "<<n1ds<<endl;
  n2ds=upperbound-n1ds; // dimension of new, non-rational forms
  if(verbose>1) cout<<"n2ds = "<<n2ds<<endl;
  if(verbose)
    {cout << "Total dimension " << dimall << " made up as follows:\n";
     cout << "dim(newforms) = " << n1ds+n2ds << " of which " << n1ds << " is rational; \n";
     cout << "dim(oldforms) = " << olddimall << " of which " << (of->olddim1) << " is rational; \n";
     cout<<endl;
   }
  delete of;

  fill_in_newform_data();
  nap=0;
}

// fill in extra data in each newforms:
void newforms::fill_in_newform_data(int everything)
{
  if(n1ds==0) return; // no work to do

  make_projcoord();    // Compute homspace::projcoord before filling in newform data
  find_jlist();
  zero_infinity = h1->chaincd(Quad::zero, Quad::one, 0, 1); // last 1 means use projcoord
  mvp=h1->maninvector(P0, 1);              // last 1 means use projcoord
  aP0 = apvec(P0);                         // vector of ap for first good principal prime
  if (verbose>1) cout << "found eigenvalues for P0="<<P0<<": "<<aP0<<endl;
  if (everything)
    for (int j=0; j<n1ds; j++)
      nflist[j].fill_in_data(j+1);

// Find the twisting primes for each newform (more efficient to do
// this here instead of within the newform constructors, as one lambda
// might work for more than one newform). NB If the level is square
// SFE=-1 then no such lambda will exist.
  if (characteristic==0)
    find_lambdas();
}

void newforms::use(const vec& b1, const vec& b2, const vector<long> eigs)
{
  if (use_nf_number==-1)
    {
      if (n1ds<upperbound)
        {
          //cout<<"Constructing newform with eigs "<<eigs<<endl;
          nflist.push_back(newform(this,b1,eigs));
          n1ds++;
        }
      else
        {
          cout << "Error in splitting eigenspaces (level "<<ideal_label(N)<<"): apparently found more ";
          cout << "1D newforms ("<< n1ds+1 <<") than the total new-dimension ("
               <<upperbound<<").\n";
          //cout<<"Extra newform has eigs "<<eigs<<endl;
          exit(1);
        }
    }
  else // store eigs and basis
    {
      nflist[use_nf_number].eigs = eigs;
      // convert basis vector from coords w.r.t. face-gens to coords w.r.t.edge-gens:
      if(hmod)
        { // we don't have a mod p mat*vec
          mat vcol(dim(b1),1);
          vcol.setcol(1,b1);
          nflist[use_nf_number].basis = reduce_modp(matmulmodp(h1->FR.coord, vcol, hmod).col(1));
        }
      else
        nflist[use_nf_number].basis = (h1->FR.coord)*b1;
      makeprimitive(nflist[use_nf_number].basis); // this is now independent of h1's denom1
    }
}

void newforms::display(int detail)
{
 if (n1ds==0) {cout<<"No newforms."<<endl; return;}
 cout << "\n"<<n1ds<<" newform(s) at level " << ideal_label(N) << " = " << gens_string(N) << ":" << endl;
 if (detail)
   for(int i=0; i<n1ds; i++)
     {
       cout<<i+1<<":\t";
       nflist[i].display();
     }
}

void newform::display(void) const
{
  cout << "basis = " << basis
      << ";\taqlist = " << aqlist
      << ";\taplist = " << aplist << endl;
 cout << "Sign of F.E. = " << sfe << endl;
 cout << "Twisting prime lambda = " << lambda << ", factor = " << lambdadot << endl;
 cout << "L/P ratio    = " << loverp << ", cuspidal factor = " << cuspidalfactor << endl;
 cout << "Integration matrix = [" << a << "," << b << ";" << c << "," << d << "], factor   = " << matdot << endl;
}

void newforms::list(long nap)
{
  string idlabel = ideal_label(N), idgens = gens_string(N), flabel = field_label();
  for(int i=0; i<n1ds; i++)
    {
      cout << flabel << " " << idlabel << " " << codeletter(i) << " " << idgens << " 2 ";  // last is weight
      int bc = nflist[i].is_base_change();
      if (bc)
        {
          int bcd = nflist[i].base_change_discriminant();
          cout << bcd;
        }
      else
        {
          int bct = nflist[i].is_base_change_twist();
          if (bct)
            {
              // NB if we have not enough inert a_P we might not be
              // able to determine the discriminant; the following
              // will return 1 in this case.
              int bcd = nflist[i].base_change_twist_discriminant();
              cout << (-bcd);
            }
          else
            cout << "0";
        }
      cout << " ";
      cout << nflist[i].is_CM() << " ";
      nflist[i].list(nap);
      cout << endl;
    }
}

void newform::list(long nap) const
{
  if(nap==-1) nap=aplist.size();
  cout << sfe << " " << loverp << " ";
  cout << "[";
  vector<long>::const_iterator ai;
  for(ai=aqlist.begin(); ai!=aqlist.end(); ++ai)
    {
      if(ai!=aqlist.begin()) cout<<",";
      cout<<(*ai);
    }
  cout << "] ";
  // The x here is essentially a place-holder representing the Hecke
  // field defining polynomial so that the output here will be
  // consistent with newforms whose Hecke field has degree >1
  cout << "x ";
  cout << "[";
  for(ai=aplist.begin(); ai!=aplist.begin()+nap && ai!=aplist.end(); ++ai)
    {
      if(ai!=aplist.begin()) cout<<",";
      cout<<(*ai);
    }
  cout <<"]";
}

// Sorting functions

// Compare two integers, using the natural order ...-2,-1.0,1,2,...
int less_ap(long a, long b)
{
  return sign(a-b);
}

// Compare two integer vectors lexicographically, using less_ap():
int less_apvec(const vector<long>& v, const vector<long>& w)
{
  vector<long>::const_iterator vi=v.begin(), wi=w.begin();
  while(vi!=v.end())
    {
      int s = less_ap(*vi++,*wi++);
      if(s) return s;
    }
  return 0;
}

struct newform_eigs_comparer {
  bool operator()(const newform& f, const newform& g)
  {
    return less_apvec(f.eigs,g.eigs)==-1;
  }
}
  less_newform_eigs;

struct newform_aplist_comparer {
  bool operator()(const newform& f, const newform& g)
  {
    return less_apvec(f.aplist,g.aplist)==-1;
  }
}
  less_newform_lmfdb;

void newforms::sort_eigs(void)
{
  ::sort(nflist.begin(),nflist.end(),less_newform_eigs);
}

void newforms::sort_lmfdb(void)
{
  ::sort(nflist.begin(),nflist.end(),less_newform_lmfdb);
}

// Replaces matrix "coord" in homspace with "projcoord"
//
// Each has ngens rows (ngens = number of edges modulo
// edge-relations).  coord has rk columns, and the i'th row of coord
// gives the coordinates of the i'th generating edge with respect to
// the basis modulo face-realations, these begin implicitly scaled by
// denom1.  projcoord has n1ds columns, and the i'th row gives the
// coords with respect to the (partial) basis of eigenforms; its
// columns are primitive and uniquely determined (up to sign),
// independent of the homology basis, and in particular of denom1.

void newforms::make_projcoord()
{
  h1->projcoord.init(h1->ngens,n1ds);
  for (int j=1; j<=n1ds; j++)
    h1->projcoord.setcol(j, nflist[j-1].basis);
}

// Second constructor, in which data is read from file in directory
// newforms. If not, recreates it from eigs

void newforms::createfromdata()
{
  if(verbose)
    cout << "Retrieving newform data for N = " << ideal_label(N) << endl;

// Read newform data from file into eigdata structure.

  eigdata filedata(N, N, -1, verbose>1, characteristic);  // neigs=-1 means get ALL from file

// Extract number of newforms and their eigenvalues from this.

  nnflist=n1ds=filedata.nforms;
  n2ds=filedata.nforms2;
  if(verbose>1) cout << " found "<<n1ds << " newforms for N = " << ideal_label(N) << endl;

 // construct the newforms from this data
  for(int i=0; i<n1ds; i++)
    {
      if(verbose>1)
        {
          cout << " constructing newform # " << i << endl;
          cout << " intdata  = "<< filedata.intdata[i] <<endl;
          cout << " Quaddata = "<< filedata.Quaddata[i] <<endl;
          cout << " aqs = "<< filedata.aqs[i] <<endl;
          cout << " aps = "<< filedata.aps[i] <<endl;
        }
      nflist.push_back(newform(this,filedata.intdata[i],filedata.Quaddata[i],
                               filedata.aqs[i],filedata.aps[i]));
    }
  nap=filedata.nap;
}

void newforms::makebases()
{
  if(!h1) makeh1plus();  // create the homology space
  sort_eigs();   // sort the newforms by their eigs list for efficient basis recovery
  form_finder splitspace(this, 1, nap, 0, 1, 0, verbose);
  if(verbose) cout<<"About to recover "<<n1ds<<" newform bases (nap="<<nap<<")"<<endl;
  for (use_nf_number=0; use_nf_number<n1ds; use_nf_number++)
    {
      if (verbose) cout<<"Recovering newform #"<<(use_nf_number+1)
                       <<", eigs "<<nflist[use_nf_number].eigs<< "...\n";
      splitspace.splitoff(nflist[use_nf_number].eigs);
    }
  if(verbose) cout<<"Finished recovering newform bases, resorting back into lmfdb order..."<<endl;
  sort_lmfdb();
  if(verbose>1) cout<<"Filling in newform data..."<<endl;
  fill_in_newform_data(0);
  if(verbose) cout<<"Finished makebases()"<<endl;
}

void newforms::getoneap(Quadprime& P, int verbose, int store)
{
  vector<long> apv=apvec(P);
  int vp = val(P, N);

  if(verbose)
    {
      string PQ = (vp>0? "Q": "P");
      cout<<PQ<<" = "<<P<<" = "<<gens_string(P)<<"\tN("<<PQ<<") = "<<P.norm()<<"\t";
    }
  for (int i=0; i<n1ds; i++)
    {
      int ap = apv[i];
      int cp = (vp==0?ap:(vp==1?-ap:0));
      if(verbose)
        cout<<setw(5)<<ap<<" ";
      if(store)
        nflist[i].aplist.push_back(cp);
    }
  if(verbose)
    cout << endl;
}

void newforms::getap(int first, int last, int verbose)
{
  if (n1ds==0) return;
  int nQP = Quadprimes::list.size();
  if(last>nQP)
    {
      last=nQP;
      cout<<"Cannot compute more than "<<nQP
          <<" ap since we only have that many primes precomputed"<<endl;
    }
  if(last<=nap)
    {
      cout<<"Already have "<<nap <<" ap " << "at level "<<N<<" so no need to compute more"<<endl;
    }
  // now nap < last <= nQP
  vector<Quadprime>::iterator pr = Quadprimes::list.begin()+first-1;
  while((pr!=Quadprimes::list.end()) && (nap<last))
    {
      getoneap(*pr++,verbose);
      nap++;
    }
}

void newforms::output_to_file(string eigfile) const
{
  int echo=0;
  ofstream out;
  out.open(eigfile.c_str());
  // Line 1
  out<<n1ds<<" "<<n2ds<<endl;
  if(echo) cout<<n1ds<<" "<<n2ds<<endl;
  if(!n1ds)
    {
      out.close();
      return;
    }
  // Line 2
  out<<nap<<endl;  if(echo) cout<<nap<<endl;
  vector<newform>::const_iterator f;
  vector<long>::const_iterator ap;
  // Line 3: SFEs
  for(f=nflist.begin(); f!=nflist.end(); ++f)
    {
      out<<setw(5)<<(f->sfe)<<" ";
      if(echo) cout<<setw(5)<<(f->sfe)<<" ";
    }
  out<<endl;  if(echo) cout<<endl;
  // Line 4: pdots
  for(f=nflist.begin(); f!=nflist.end(); ++f)
    {
      out<<setw(5)<<(f->pdot)<<" ";
      if(echo) cout<<setw(5)<<(f->pdot)<<" ";
    }
  out<<endl;  if(echo) cout<<endl;
  // Line 5: dp0s
  for(f=nflist.begin(); f!=nflist.end(); ++f)
    {
      out<<setw(5)<<(f->dp0)<<" ";
      if(echo) cout<<setw(5)<<(f->dp0)<<" ";
    }
  out<<endl;  if(echo) cout<<endl;
  // Line 6: cuspidal factors
  for(f=nflist.begin(); f!=nflist.end(); ++f)
    {
      out<<setw(5)<<(f->cuspidalfactor)<<" ";
      if(echo) cout<<setw(5)<<(f->cuspidalfactor)<<" ";
    }
  out<<endl;  if(echo) cout<<endl;
  // Line 7: lambdas
  for(f=nflist.begin(); f!=nflist.end(); ++f)
    {
      Quad lambda = f->lambda;
      out<<setw(5)<< lambda.re()<<" "<< lambda.im()<<" ";
      if(echo) cout<<setw(5)<< lambda.re()<<" "<< lambda.im()<<" ";
    }
  out<<endl;  if(echo) cout<<endl;
  // Line 8: lambdadots
  for(f=nflist.begin(); f!=nflist.end(); ++f)
    {
      out<<setw(5)<<(f->lambdadot)<<" ";
      if(echo) cout<<setw(5)<<(f->lambdadot)<<" ";
    }
  out<<endl;  if(echo) cout<<endl;
  // Lines 9,10,11,12: a,b,c,d:
  Quad a;
  for(f=nflist.begin(); f!=nflist.end(); ++f)
    {
      a = f->a;
      out<<setw(5)<< a.re()<<" "<< a.im()<<" ";
      if(echo) cout<<setw(5)<< a.re()<<" "<< a.im()<<" ";
    }
  out<<endl;  if(echo) cout<<endl;
  for(f=nflist.begin(); f!=nflist.end(); ++f)
    {
      a = f->b;
      out<<setw(5)<< a.re()<<" "<< a.im()<<" ";
      if(echo) cout<<setw(5)<< a.re()<<" "<< a.im()<<" ";
    }
  out<<endl;  if(echo) cout<<endl;
  for(f=nflist.begin(); f!=nflist.end(); ++f)
    {
      a = f->c;
      out<<setw(5)<< a.re()<<" "<< a.im()<<" ";
      if(echo) cout<<setw(5)<< a.re()<<" "<< a.im()<<" ";
    }
  out<<endl;  if(echo) cout<<endl;
  for(f=nflist.begin(); f!=nflist.end(); ++f)
    {
      a = f->d;
      out<<setw(5)<< a.re()<<" "<< a.im()<<" ";
      if(echo) cout<<setw(5)<< a.re()<<" "<< a.im()<<" ";
    }
  out<<endl;  if(echo) cout<<endl;
  // Line 13: matdots
  for(f=nflist.begin(); f!=nflist.end(); ++f)
    {
      out<<setw(5)<<(f->matdot)<<" ";
      if(echo) cout<<setw(5)<<(f->matdot)<<" ";
    }
  out<<endl;  if(echo) cout<<endl;
  out<<endl;  if(echo) cout<<endl;
  for(int i=0; i<nwq; i++)
    {
      for(f=nflist.begin(); f!=nflist.end(); ++f)
	{
	  out<<setw(5)<<(f->aqlist)[i]<<" ";
	  if(echo) cout<<setw(5)<<(f->aqlist)[i]<<" ";
	}
      out<<endl;
      if(echo) cout<<endl;
    }
  out<<endl;  if(echo) cout<<endl;
  for(int i=0; i<nap; i++)
    {
      for(f=nflist.begin(); f!=nflist.end(); ++f)
	{
	  out<<setw(5)<<(f->aplist)[i]<<" ";
	  if(echo) cout<<setw(5)<<(f->aplist)[i]<<" ";
	}
      out<<endl;      if(echo) cout<<endl;
    }
  out.close();
}

// For each newform we want a pivotal index j in [1,ngens] such that
// the j'th coordinate is nonzero, so that we can compute the Hecke
// eigenvalue a_p by only computing the image of the j'th
// edge-generator.

void newforms::find_jlist()
{
  if(verbose>1)
    cout<<"Finding pivotal indices..."<<flush;

  int i, j, ok=0; j0=0;

  // First we see whether a single j works for all the newforms:

  for(j=1; (!ok)&&(j<=h1->ngens); j++)
    {
      ok=1;
      for (i=0; (i<n1ds)&&ok; i++)
        ok=(nflist[i].basis[j]!=0);
      if(ok) j0=j;
    }
  if(ok)
    {
      jlist.insert(j0);
      for (i=0; i<n1ds; i++)
	{
          newform& nfi = nflist[i];
	  nfi.j0 = j0;
          nfi.fac = nfi.basis[j0];
          if(nfhmod) nfi.facinv=invmod(nfi.fac,nfhmod);
	}
      if(verbose>1)
        cout<<"index j0="<<j0<<" works as a pivot for all newforms"<<endl;
      return;
    }

  if(verbose>1)
    cout<<"...failed to find a pivotal index which works for all newforms..." <<flush;

  // Instead, find a set of pivots to use:
  for (i=0; i<n1ds; i++)
    {
      newform& nfi = nflist[i];
      vec& bas = nfi.basis;
      j=1; while(bas[j]==0) j++;
      jlist.insert(j); // jlist is a set so this will do nothing if it is already there
      nfi.j0 = j;
      nfi.fac = bas[j];
      if(nfhmod) nfi.facinv=invmod(nfi.fac,nfhmod);
    }
  if(verbose>1)
    cout<<"set of pivotal indices = "<<jlist<<endl;
}

//#define DEBUG_APVEC

// The following utility does the following.  Given an integer ind:
// - if ind>0 it adds the ind'th row of pcd to imagej;
// - if ind<0 it subtracts the |ind|'th row of pcd from imagej;
// - if ind=0 it leaves imagej unchaged.
// if hmod is nonzero the vector addition is done modulo hmod.

void update(const mat& pcd, vec& imagej, long ind, long hmod)
{
  if (ind==0) return;
  vec part = (ind>0? pcd.row(ind): -pcd.row(-ind));
#ifdef DEBUG_APVEC
  cout<<"--adding "<<part<<" (ind="<<ind<<") to imagej, ";
#endif
  imagej = reduce_modp(imagej + part, hmod);
#ifdef DEBUG_APVEC
  cout<<"updated imagej is "<<imagej<<endl;
#endif
}

//#define DEBUG_APVEC
vector<long> newforms::apvec(Quadprime& P)  // computes a[P] for each newform, for principal P
{
  Quad p = P.gen();
#ifdef DEBUG_APVEC
  cout<<"In apvec with P = "<<P<<endl;
#endif
  vector<long> apv(n1ds);
  long ap, normp=I2long(P.norm());

  int vp = val(P,N);

  if (vp>0) // bad prime, we already know the eigenvalues
    {
      for (int i=0; i<n1ds; i++)
        {
          int ip = find(plist.begin(),plist.end(),P) - plist.begin();
          long aq = nflist[i].aqlist[ip];
          if(!((aq==1)||(aq==-1)))
            {
              cout<<"Error: Atkin-Lehner eigenvalue "<<aq<<" for Q="<<P
                  <<" for form # "<<(i+1)<<" is neither +1 nor -1"<<endl;
              cout<<"------------------------"<<endl;
              nflist[i].display();
              cout<<"------------------------"<<endl;
              exit(1);
            }
          apv[i] = aq;
        }
#ifdef DEBUG_APVEC
      cout<<"Bad prime: aq list = "<<apv<<endl;
#endif
      return apv;
    }

  // now P is a good prime

  long maxap=(long)(2*sqrt((double)normp)); // for validity check

  matop Tp;
  if (!Quad::is_Euclidean)  // not needed in Euclidean case where we
    Tp = HeckeOp(P, N);       // use Manin-Heilbronn matrices instead

  map<int,vec> images; // [j,v] stores image of j'th M-symbol in v
                       // (so we don't compute any more than once)
  Quad a,b,c,q,u1,u2,u3;

  // Compute the image of the necessary M-symbols (hopefully only one)
#ifdef DEBUG_APVEC
  cout<<"Computing images of M-symbols"<<endl<<flush;
  cout<<"jlist = "<<jlist<<endl;
  cout<<"j's for each newform: ";
  for (i=0; i<n1ds; i++) cout<<nflist[i].j0<<" ";
  cout<<endl;
  cout<<"factors for each newform: ";
  for (i=0; i<n1ds; i++) cout<<nflist[i].fac<<" ";
  cout<<endl;
#endif

  for(std::set<long>::const_iterator jj=jlist.begin(); jj!=jlist.end(); ++jj)
    {
      vec imagej=vec(n1ds); // initialised to 0
      long j=*jj; // from 1
      pair<long, int> st = h1->ER.symbol_number_and_type(h1->ER.gen(j));
      long s_number = st.first;  // (c:d) symbol number
      int s_type    = st.second; // symbol type (negative for singular edges)

      Quad u, v;
      h1->P1.make_symb(s_number, u, v);

#ifdef DEBUG_APVEC
      cout<<"Computing image under T("<<p<<") of "<<j<<"'th M-symbol"
          <<" = ("<<u<<":"<<v<<")_"<<s_type<<" ..."<<flush;
#endif

      // Now we compute the projected image of symbol s=(u:v)_t under
      // T_P.

      if (Quad::is_Euclidean)
        {

          // This code is for Euclidean fields only, using Manin-Heilbronn
          // matrices: Loop over residues res mod P and for each res
          // compute several M-symbol image parts (u1:v1).  Accumulate the
          // associated vectors in vec imagej using the utility
          // update(projcoord, imagej, ind, nfhmod), where (u1:v1) is the
          // ind'th symbol.

          // Since this code is only used in the Euclidean case,
          // s_type=0 and can be ignored.

          mat& pcd = h1->projcoord;

          // Matrix [1,0;0,p]
          long ind = h1->ER.coords(h1->index(u,p*v));
#ifdef DEBUG_APVEC
          cout<<"u1="<<u<<", u2="<<p*v<<", ind="<<ind<<endl;
#endif
          update(pcd,imagej,ind,nfhmod);

          // Matrix [p,0;0,1]
          ind = h1->ER.coords(h1->index(p*u,v));
#ifdef DEBUG_APVEC
          cout<<"u1="<<p*u<<", u2="<<v<<", ind="<<ind<<endl;
#endif
          update(pcd,imagej,ind,nfhmod);

          // Other matrices, several for each nonzero residue b mod p
          vector<Quad> resmodp = P.residues();
          vector<Quad>::const_iterator res=resmodp.begin();
          while(res!=resmodp.end())
            {
              b = *res++;
              if(b.is_zero()) continue; // handled above as special case
              a = -p;
              u1=u*p; u2=v-u*b;
              ind = h1->ER.coords(h1->index(u1,u2));
#ifdef DEBUG_APVEC
              cout<<"Residue class "<<b<<": ";
              cout<<"a="<<a<<", b="<<b<<", u1="<<u1<<", u2="<<u2<<", ind="<<ind<<endl;
#endif
              update(pcd,imagej,ind,nfhmod);
              while(!b.is_zero())
                {
                  q=a/b; c=a-b*q; u3=q*u2-u1;
                  a=-b; b=c; u1=u2; u2=u3;
                  ind = h1->ER.coords(h1->index(u1,u2));
#ifdef DEBUG_APVEC
                  cout<<"a="<<a<<", b="<<b<<", u1="<<u1<<", u2="<<u2<<", ind="<<ind<<endl;
#endif
                  update(pcd,imagej,ind,nfhmod);
                }
#ifdef DEBUG_APVEC
              cout<<" partial image after term is "<<imagej<<endl;
#endif
            }
          images[j]=imagej;

        } // end of Euclidean case

      else // non-Euclidean case, apply T_p directly to the modular symbol represented by j

        {
          mat22 M = h1->P1.lift_to_SL2(s_number);
          modsym m(M, s_type);
#ifdef DEBUG_APVEC
          cout<<"\n coset rep = " << M << "\n modular symbol m = "<<m<<endl;
          cout<<"chain(m, 0) = "<<h1->chain(m,0)<<endl;
          cout<<"chain(m, 1) = "<<h1->chain(m,1)<<endl;
#endif
          images[j] = h1->applyop(Tp, m, 1);
        }

#ifdef DEBUG_APVEC
      cout<<" image " << j << " is "<<images[j]<<endl;
#endif

    }

// recover eigenvalues:

  for (int i=0; i<n1ds; i++)
    {
      int j0 = nflist[i].j0;
      int fac = nflist[i].fac;
      int top = images[j0][i+1];
      // The eigenvalue is now top/fac (which should divide exactly)
      if(nfhmod)
        ap=mod(xmodmul(top,nflist[i].facinv,nfhmod), nfhmod);
      else
        {
#ifdef DEBUG_APVEC
          cout << "ap   = " << top << "/" << fac << " = " <<top/fac<<endl;
#endif
          if (top%fac !=0)
            {
              cout<<"Problem in apvec: for newform #"<<(i+1)<<", with pivotal index "<<j0<<" and pivot "<<fac<<endl;
              cout<<"\timage list = "<<images[j0]<< " has "<<(i+1)<<" entry "<<top<<" which is not divisible by pivot "<<fac<<endl;
              cout<<flush;
            }
           ap = top/fac;
        }
      apv[i]=ap;
// #ifdef DEBUG_APVEC
//       cout << "ap = " << ap << endl;
// #endif
// check it is in range (in characteristic 0 only):
      if (characteristic==0)
        if((ap>maxap)||(-ap>maxap))
          {
            cout<<"Error:  eigenvalue "<<ap<<" for P="<<P
                <<" for form # "<<(i+1)<<" is outside valid range "
                <<-maxap<<"..."<<maxap<<endl;
            exit(1);
          }
    }
#ifdef DEBUG_APVEC
      cout << "Good prime: ap list = " << apv << endl;
#endif
  return apv;
}

//end of newforms.cc
