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

#define MAXDEPTH 20 // maximum depth for splitting off eigenspaces

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
    {
      basis = (nf->h1->FR.coord)*v;
      makeprimitive(basis); // this is now independent of h1's denom1
    }

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

// In detail: this assumes we have the list of eigs; it extracts
// aqlist from this (but this may not be complete for even class
// number) and the extra data (sfe, integration data etc) but *not*
// aplist:

//  - aq and sfe from eigs
//  - L/P
//  - manin vector data
//  - integration matrix  // using find_matrix()

// The reason for the parameter j is to access data computed and
// stored in the newforms class

void newform::data_from_eigs(int j)
{
  sfe = 0;

  // Atkin-Lehner eigenvalues

  // TODO (for even h)
  // Does not yet account for the fact that in even class number we
  // may have nwq<npdivs and need to work harder! Also in this case
  // the first n2r eigs will be +1 and can be ignored.
  if (nf->characteristic==0)
    {
      // extract Atkin-Lehner eigs from first entries of eigs list:
      copy(eigs.begin(), eigs.begin()+(nf->nwq), back_inserter(aqlist));

      // Sign of functional equation = minus product of all A-L eigenvalues
      sfe = std::accumulate(aqlist.begin(),aqlist.end(),-1,std::multiplies<long>());
    }

  // compute L/P as n_F({0,oo})
  int pdot0 = abs(nf->zero_infinity[j]);
  loverp =  rational(pdot0, (Quad::nunits) * cuspidalfactor);

  // compute L/P again using Manin vector
  dp0  =  1 + (nf->nP0) - nf->aP0[j-1];  // aP0 is based at 0
  pdot = abs(nf->mvp[j]);
  rational loverp_mvp(pdot, dp0 * (Quad::nunits) * cuspidalfactor);

  if (nf->characteristic>0)
    return;

  // Check they agree:

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

// When a newform has been read from file, we have the aqlist and
// aplist but not the sequence of eigs in order.  This is needed
// both for recovering the basis vector from the h1 (in case we want
// to compute more ap), and for computing oldform multiplcities.

void newform::eigs_from_data()
  // recreate eigs list (in case we need to recover basis vector):
  // start with unramified char eigs (all +1), then Tp eigs ap for
  // good p
{
  // cout<<"In eigs_from_data, aplist = "<<aplist<<endl;
  int ch(nf->characteristic);
  if (ch == 0)      // the first n2r eigs are all +1
    eigs.resize(nf->n2r, +1);
  else
    eigs.resize(0, +1);

  // Get a_P or a_{P^2} from the a_P in aplist, for good P
  vector<Quadprime>::const_iterator pr=Quadprimes::list.begin();
  vector<long>::const_iterator api=aplist.begin();
  while (((int)eigs.size() < nf->nap+nf->n2r) && (api!=aplist.end()))
    {
      Quadprime P = *pr;
      QUINT normP = P.norm();
      while ((P.divides(nf->N)) || (ch>0 && (normP%ch==0)))
        {
          // cout<<" - P = "<<P<<": bad prime, skipping"<<endl;
          ++pr;
          ++api;
          P = *pr;
          normP = P.norm();
        }
      long ap = *api;
      if (!P.has_square_class()) // eigenvalues of T_{P^2}, not T_P
        {
          ap = ap*ap + normP;
        }
      eigs.push_back(ap);
      // cout<<" - P = "<<P<<": eig = "<<ap<<endl;
      ++pr; ++api;
    }
  // cout<<" eigs_from_data produced eigs = "<<eigs<<endl;
}

// For M a *multiple* of this level N, make the list of eigs
// appropriate for the higher level, taking into account the primes
// P (if any) dividing M but not N. For such P we delete the a_P
// from the sublist of T_P eigenvalues and insert a 0 into the W_Q
// eigenvalues. The oldform constructor will deal with this.
vector<long> newform::oldform_eigs(Qideal& M)
{
  assert (nf->N.divides(M));

  eigs_from_data();
  vector<long> M_eigs;

  if (nf->verbose)
    {
      cout<<"Making oldform eigs at level "<<ideal_label(M)<<" from eigs at level "<<ideal_label(nf->N)<<endl;
      cout<<" - input eigs: "<<eigs<<endl;
    }
  // insert eigs for central characters:
  if (nf->characteristic == 0)  // the first n2r eigs are all +1
    {
      M_eigs.resize(nf->n2r, +1);
    }

  vector<long>::const_iterator ei = eigs.begin() + (nf->n2r);
  if (nf->nwq>0)
    {
      vector<Quadprime> allMprimes = M.factorization().sorted_primes();
      vector<Quadprime> Mprimes = make_badprimes(M, allMprimes);
      for (vector<Quadprime>::iterator Qi = Mprimes.begin(); Qi!=Mprimes.end(); ++Qi)
        {
          Quadprime Q = *Qi;
          if (Q.divides(nf->N))
            {
              if (nf->verbose>1)
                cout << " keeping W eigenvalue "<<(*ei)<< " at "<<Q<<endl;
              M_eigs.push_back(*ei++);
            }
          else
            {
              if (nf->verbose>1)
                cout << " dummy W eigenvalue 0 at "<<Q<<endl;
              M_eigs.push_back(0); // dummy value, oldforms class will handle this
            }
        }
    }

  for (vector<Quadprime>::const_iterator Pi = Quadprimes::list.begin();
       Pi != Quadprimes::list.end() && ei!=eigs.end(); ++Pi)
    {
      Quadprime P = *Pi;
      if (!P.divides(nf->N))
        {
          if (!P.divides(M)) // else this T_P eigenvalue is ignored
            {
              if (nf->verbose>1)
                cout << " keeping eigenvalue "<<(*ei)<< " at "<<P<<endl;
              M_eigs.push_back(*ei);
            }
          ++ei;
        }
    }
  if (nf->verbose)
    {
      cout<<" - output eigs: "<<M_eigs<<endl;
    }
  return M_eigs;
}


newform::newform(newforms* nfs,
                 const vector<int>& intdata, const vector<Quad>& Quaddata,
                 const vector<long>& aq, const vector<long>& ap)
{
  nf=nfs;
  Qideal N(nf->N);
  int ch(nf->characteristic);

  if (ch == 0)
    {
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
        cout<<"Problem in data on file for level "<<N<<": sfe = "<<sfe<<" and aqlist = "<<aqlist<<", but minus product of latter is "<<newsfe<<endl;
    }

  aplist = ap;

  // recreate eigs list (in case we need to recover basis vector):

  eigs_from_data();
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

newforms::newforms(const Qideal& iN, int disp, long ch)
  : maxdepth(MAXDEPTH), N(iN), verbose(disp), n2r(Quad::class_group_2_rank), characteristic(ch)
{
  is_square = N.is_square();

  // nulist is a list of n2r ideals coprime to N whose classes generate the 2-torsion
  if ((characteristic==0) && (Quad::class_group_2_rank > 0))
    nulist = make_nulist(N);

  // badprimes is a list of all primes Q|N
  allbadprimes = N.factorization().sorted_primes();

  // badprimes is a list of primes Q|N such that [Q^e] is square
  if (characteristic==0)
    badprimes = make_badprimes(N, allbadprimes);
  // nwq = badprimes.size();
  nwq = 0; // prevents any W_Q being used for splitting

  // goodprimes is a list of at least nap good primes (excluding those
  // dividing characteristic if >0), includinge at least one principal
  // one which has index iP0;

  nap = 20;
  goodprimes = make_goodprimes(N, nap, iP0, characteristic);
  nap = goodprimes.size(); // it may be > original nap
  if (nap!=20)
    cout<<" nap changed to "<<nap<<" since goodprimes = "<<goodprimes<<endl;
  P0 = goodprimes[iP0];
  nP0 = I2long(P0.norm());

// P0 is the smallest good principal prime: and iP0 its index (in
// plist, which starts with the bad primes and then the good
// primes in order).  P0 must be principal since we have only
// implemented maninvector() for principal primes.

  if (verbose>1)
    {
      // if (characteristic==0)
      //   cout << "bad primes used: "<< badprimes<<endl;
      cout << "good primes used: "<<goodprimes<<endl;
    }

  h1=0;
  of=0;
  nfhmod=0;
}

// instantiations of virtual functions required by the splitter_base class:
mat newforms::opmat(int i, int dual, int verb)
{
  return h1->calcop(h1matops[i],dual,verb);
}

vec newforms::opmat_col(int i, int j, int verb)
{
  return h1->calcop_col(h1matops[i],j, verb);
}

mat newforms::opmat_cols(int i, const vec& jlist, int verb)
{
  return h1->calcop_cols(h1matops[i],jlist, verb);
}

mat newforms::opmat_restricted(int i, const subspace& s, int dual, int verb)
{
  return h1->calcop_restricted(h1matops[i],s,dual,verb);
}

smat newforms::s_opmat(int i, int dual, int verb)
{
  return h1->s_calcop(h1matops[i],dual, verbose);
}

smat newforms::s_opmat_cols(int i, const vec& jlist, int verb)
{
  return h1->s_calcop_cols(h1matops[i],jlist, verbose);
}

smat newforms::s_opmat_restricted(int i, const ssubspace& s, int dual, int verb)
{
  return h1->s_calcop_restricted(h1matops[i],s,dual,0);
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

void newforms::find()
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
      if (characteristic==0)
        cout<<"Retrieving oldform data for level "<<N<<"...\n";
    }
  of = new oldforms(this);
  long olddimall = (of->olddimall);
  if(verbose)
    {
      if (characteristic==0)
        of->display();
      cout<<"Finding rational newforms...\n";
    }

  long mindepth = (characteristic==0? n2r+nwq+iP0: nap);
  n1ds = 0;
  upperbound = (characteristic==0? (h1->h1cuspdim()) - olddimall: h1->h1dim());
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
      for (int i=0; i<maxdepth; i++)
        {
          h1matops.push_back(h1matop(i));
          eigranges.push_back(eigrange(i));
        }
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
     if (characteristic==0)
       cout << "dim(oldforms) = " << olddimall << " of which " << (of->olddim1) << " is rational; \n";
     cout<<endl;
   }
  delete of;

  fill_in_newform_data();
  nap=0;
}

// fill in extra data in each newform:
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
    {
      // compute A-L eigenvalues
      for (vector<Quadprime>::iterator Qi = allbadprimes.begin(); Qi!=allbadprimes.end(); ++Qi)
        {
          vector<long> apv = apvec(AtkinLehnerOp(*Qi,N), 1);
          for (int j=0; j<n1ds; j++)
            nflist[j].aqlist.push_back(apv[j]);
        }
      for (int j=0; j<n1ds; j++)
        nflist[j].data_from_eigs(j+1);
    }

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
      if (characteristic==0)
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
  cout << "basis = " << basis;
  if (nf->characteristic==0)
    cout << ";\taqlist = " << aqlist;
  cout << ";\taplist = " << aplist << endl;
  if (nf->characteristic==0)
    {
      cout << "Sign of F.E. = " << sfe << endl;
      cout << "Twisting prime lambda = " << lambda << ", factor = " << lambdadot << endl;
      cout << "L/P ratio    = " << loverp << ", cuspidal factor = " << cuspidalfactor << endl;
      cout << "Integration matrix = [" << a << "," << b << ";" << c << "," << d << "], factor   = " << matdot << endl;
    }
}

void newforms::list(long nap)
{
  string idlabel = ideal_label(N), idgens = gens_string(N), flabel = field_label();
  for(int i=0; i<n1ds; i++)
    {
      cout << flabel << " " << idlabel << " " << codeletter(i) << " " << idgens << " 2 ";  // last is weight
      if (characteristic>0)
        cout << characteristic << " ";
      else
        {
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
        }
      nflist[i].list(nap);
      cout << endl;
    }
}

void newform::list(long nap) const
{
  if(nap==-1) nap=aplist.size();
  vector<long>::const_iterator ai;

  if (nf->characteristic==0)
    {
      cout << sfe << " " << loverp << " ";
      cout << "[";
      for(ai=aqlist.begin(); ai!=aqlist.end(); ++ai)
        {
          if(ai!=aqlist.begin()) cout<<",";
          cout<<(*ai);
        }
      cout << "] ";
    }
  // The x here is essentially a place-holder representing the Hecke
  // field defining polynomial so that the output here will be
  // consistent with newforms whose Hecke field has degree >1
  cout << "x ";
  cout << "[";
  for(ai=aplist.begin(); ai!=aplist.begin()+nap && ai!=aplist.end(); ++ai)
    {
      if(ai!=aplist.begin()) cout<<",";
      long ap = *ai;
      if ((nf->characteristic>0) && (ap<0)) // we use -999 for omitted eigs
        cout << "*";
      else
        cout << ap;
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

// try to read from file, and if no data file exists, finds from scratch and stores
void newforms::read_from_file_or_find()
{
  if (verbose>1)
    cout << " - reading newform data for level "<<ideal_label(N)<<endl;
  int ok = read_from_file();
  if (ok)
    {
      if (verbose>1)
        cout << " - successfully read newform data for level "<<ideal_label(N)<<endl;
      return;
    }
  if (verbose)
    cout << " - no newform data for level "<<ideal_label(N)<<" exists, finding newforms..."<<endl;
  find();
  if (verbose)
    cout << " - found "<<n1ds<<" newforms for level "<<ideal_label(N)<<endl;
  if (n1ds>0)
    {
      if (verbose)
        cout << " - computing eigenvalues numbers 1 to "<<nap<<"... "<<endl;
      getap(1, max(nap,25), 0);
    }
  string eigfilename = (Quad::class_number==1? eigfile(N.gen(), characteristic): eigfile(N, characteristic));
  output_to_file(eigfilename);
  if(verbose)
    {
      cout << "  finished creating and storing newforms at level " << N << endl;
      if (verbose>1)
        display();
    }
}

int newforms::read_from_file()
{
  if(verbose)
    cout << "Retrieving newform data for N = " << ideal_label(N) << endl;

// Read newform data from file

  if(verbose>1) cout << "Getting newform data for " << N << endl;
  string eigfilename = (Quad::class_number==1? eigfile(N.gen(), characteristic): eigfile(N, characteristic));
  ifstream data(eigfilename.c_str());
  if (!data)
    {
      if(verbose)
        {
          cout << "No data file for level " << ideal_label(N);
          if (characteristic) cout << " mod " << characteristic;
          cout << endl;
        }
      return 0;
    }
  data >> n1ds >> n2ds >> nap;
  if(verbose>1)
    {
      cout<<" read data for "<<n1ds<<" newforms at level "<<N
          <<", total new dimension = "<<(n1ds+n2ds)<<", nap = "<<nap<<endl;
    }
  if (n1ds==0)
    return 1;

  vector<vector<long> > aqs(n1ds), aps(n1ds), eigs(n1ds);
  vector<vector<int> > intdata(n1ds);   // sfe, pdot, dp0, cuspidalfactor, lambdadot, matdot
  vector<vector<Quad> > Quaddata(n1ds); // lambda, a, b, c, d

  int i;
  for(i=0; i<n1ds; i++)
    {
      eigs[i].resize(nap);
      aps[i].resize(nap);
      if (characteristic==0)
        {
          aqs[i].resize(allbadprimes.size());
          intdata[i].resize(6);
          Quaddata[i].resize(5);
        }
    }

  vector<vector<long> >::iterator f; long eig;

  if (characteristic==0)
    {
      // Read the auxiliary data (unless in positive characteristic):
      for (i=0; i<n1ds; i++) data>>intdata[i][0];  // sfe
      for (i=0; i<n1ds; i++) data>>intdata[i][1];  // pdot
      for (i=0; i<n1ds; i++) data>>intdata[i][2];  // dp0
      for (i=0; i<n1ds; i++) data>>intdata[i][3];  // cuspidalfactor
      for (i=0; i<n1ds; i++) data>>Quaddata[i][0]; // lambda
      for (i=0; i<n1ds; i++) data>>intdata[i][4];  // lambdadot
      for (i=0; i<n1ds; i++) data>>Quaddata[i][1]; // a
      for (i=0; i<n1ds; i++) data>>Quaddata[i][2]; // b
      for (i=0; i<n1ds; i++) data>>Quaddata[i][3]; // c
      for (i=0; i<n1ds; i++) data>>Quaddata[i][4]; // d
      for (i=0; i<n1ds; i++) data>>intdata[i][5];  // matdot

      //  Read the W-eigenvalues at level M into aqs:
      for(i=0; i<(int)allbadprimes.size(); i++)
        for(f=aqs.begin(); f!=aqs.end(); ++f)
          {
            data>>eig;
            (*f)[i]=eig;
          }
    }

  // Next read the coefficients at level M into aps:
  for(i=0; i<nap; i++)
    for(f=aps.begin(); f!=aps.end(); ++f)
      {
        data>>eig;
        (*f)[i]=eig;
      }

  data.close();

  if(verbose>1)
    {
      cout << "Finished reading newform data for level " << N << endl;
      if (characteristic==0)
        {
          cout << "aqs = " << endl;
          for(i=0; i<n1ds; i++) cout<<i<<": "<<aqs[i]<<endl;
        }
      cout << "aps = " << endl;
      for(i=0; i<n1ds; i++) cout<<i<<": "<<aps[i]<<endl;
      if (characteristic==0)
        {
          cout << "intdata = " << endl;
          for(i=0; i<n1ds; i++) cout<<i<<": "<<intdata[i]<<endl;
          cout << "Quaddata = " << endl;
          for(i=0; i<n1ds; i++) cout<<i<<": "<<Quaddata[i]<<endl;
        }
    }

// Extract number of newforms and their eigenvalues from this.

  if(verbose>1) cout << " read "<<n1ds << " newforms for N = " << ideal_label(N) << endl;

 // construct the newforms from this data
  for(int i=0; i<n1ds; i++)
    {
      if(verbose>1)
        {
          cout << " constructing newform # " << i << endl;
          cout << " intdata  = "<< intdata[i] <<endl;
          cout << " Quaddata = "<< Quaddata[i] <<endl;
          cout << " aqs = "<< aqs[i] <<endl;
          cout << " aps = "<< aps[i] <<endl;
        }
      // the constructor here calls eigs_from_data() so these newforms have their eigs lists
      nflist.push_back(newform(this,intdata[i],Quaddata[i], aqs[i],aps[i]));
    }
  return 1;
}

void newforms::makebases()
{
  if(!h1) makeh1plus();  // create the homology space
  sort_eigs();   // sort the newforms by their eigs list for efficient basis recovery
  for (int i=0; i<maxdepth; i++)
    {
      h1matops.push_back(h1matop(i));
      eigranges.push_back(eigrange(i));
    }
  form_finder splitspace(this, 1, maxdepth, 0, 1, 0, verbose);
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
      int cp = ((vp==0)||(characteristic>0)? ap : (vp==1? -ap : 0));
      if(verbose)
        {
          if ((characteristic>0) && (ap==-999))
            cout<<setw(5)<<"?"<<" ";
          else
            cout<<setw(5)<<ap<<" ";
        }
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
      getoneap(*pr++, verbose);
      nap++;
    }
}

void newforms::output_to_file(string eigfile) const
{
  int echo=0; // Set to 1 to echo what is written to the file for debugging
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

  if (characteristic==0)
    {

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
  for(int i=0; i<(int)allbadprimes.size(); i++)
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
    } // end of if (characteristic==0) block
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

// compute eigenvalues given the image images[j] for each j in jlist
vector<long> newforms::apvec_from_images(map<int,vec> images, long maxap, const string& name)
{
  vector<long> apv(n1ds);

  for (int i=0; i<n1ds; i++)
    {
      int j0 = nflist[i].j0;
      int fac = nflist[i].fac;
      int top = images[j0][i+1];
      long ap;
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
      if (characteristic>0)
        apv[i] = posmod(ap, characteristic);
      else
        {
          apv[i]=ap;
          // check it is in range (in characteristic 0 only):
          if (abs(ap)>maxap)
            {
              cout<<"Error:  eigenvalue "<<ap<<" for operator "<<name
                  <<" for form # "<<(i+1)<<" is outside valid range "
                  << -maxap<<"..."<<maxap<<endl;
              exit(1);
            }
        }
    }
  return apv;
}


// compute eigenvalue of op for each newform and check that it is <= maxap
vector<long> newforms::apvec(const matop& op, long maxap)
{
#ifdef DEBUG_APVEC
  cout<<"In apvec with operator "<<op.name()<<endl;
#endif
  // Compute the image images[j] of the j'th symbol under op, for all necessary j.
  map<int,vec> images;
  for(std::set<long>::const_iterator jj=jlist.begin(); jj!=jlist.end(); ++jj)
    {
      long j=*jj; // from 1
      pair<long, int> st = h1->ER.symbol_number_and_type(h1->ER.gen(j));
      long s_number = st.first;  // (c:d) symbol number
      int s_type    = st.second; // symbol type (negative for singular edges)

      // Now we compute the projected image of symbol under op

      modsym m(h1->P1.lift_to_SL2(s_number), s_type);
      images[j] = h1->applyop(op, m, 1);
    }

  vector<long> apv = apvec_from_images(images, maxap, op.name());
#ifdef DEBUG_APVEC
  cout << "eigenvalue list = " << apv << endl;
#endif
  return apv;
}

// Special code for T_P for Euclidean fields, for good P only

// The following utility does the following.  Given an integer ind:
// - if ind>0 it adds the ind'th row of pcd to imagej;
// - if ind<0 it subtracts the |ind|'th row of pcd from imagej;
// - if ind=0 it leaves imagej unchaged.
// if hmod is nonzero the vector addition is done modulo hmod.

void update(const mat& pcd, vec& imagej, long ind, long hmod)
{
  if (ind==0) return;
  vec part = (ind>0? pcd.row(ind): -pcd.row(-ind));
  imagej = reduce_modp(imagej + part, hmod);
}

// compute eigenvalue at P for each newform (good P, Euclidean) and check that it is <= maxap
vector<long> newforms::apvec_euclidean(Quadprime& P, long maxap)
{
  assert (Quad::is_Euclidean && "field must be Euclidean in apvec_euclidean()");
  assert (val(P,N)==0 && "P must be good in apvec_euclidean()");
  Quad p = P.gen();

  //images[j] is the image of the j'th M-symbol
  map<int,vec> images;
  Quad a,b,c,q,u1,u2,u3;

  // Compute the image of the necessary M-symbols (hopefully only one)

  for(std::set<long>::const_iterator jj=jlist.begin(); jj!=jlist.end(); ++jj)
    {
      vec imagej=vec(n1ds); // initialised to 0
      long j=*jj; // from 1
      // Since this code is only used in the Euclidean case,
      // all symbols have type 0
      long s_number = h1->ER.gen(j);  // (c:d) symbol number
      Quad u, v;
      h1->P1.make_symb(s_number, u, v);

      // Now we compute the projected image of symbol s=(u:v)_t under
      // T_P.

      // This code is for Euclidean fields only, using Manin-Heilbronn
      // matrices: Loop over residues res mod P and for each res
      // compute several M-symbol image parts (u1:v1).  Accumulate the
      // associated vectors in vec imagej using the utility
      // update(projcoord, imagej, ind, nfhmod), where (u1:v1) is the
      // ind'th symbol.

      mat& pcd = h1->projcoord;

      // Matrix [1,0;0,p]
      long ind = h1->ER.coords(h1->index(u,p*v));
      update(pcd,imagej,ind,nfhmod);

      // Matrix [p,0;0,1]
      ind = h1->ER.coords(h1->index(p*u,v));
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
          update(pcd,imagej,ind,nfhmod);
          while(!b.is_zero())
            {
              q=a/b; c=a-b*q; u3=q*u2-u1;
              a=-b; b=c; u1=u2; u2=u3;
              ind = h1->ER.coords(h1->index(u1,u2));
              update(pcd,imagej,ind,nfhmod);
            }
        }
      images[j]=imagej;
    }

  vector<long> apv = apvec_from_images(images, maxap, opname(P,N));
#ifdef DEBUG_APVEC
  cout << "eigenvalue list = " << apv << endl;
#endif
  return apv;
}

// compute a[P] for each newform
vector<long> newforms::apvec(Quadprime& P)
{
#ifdef DEBUG_APVEC
  cout<<"In apvec with P = "<<P<<endl;
#endif
  vector<long> apv(n1ds);
  long normp=I2long(P.norm());
  int i, vp = val(P,N);

  if ((characteristic>0) && ((vp>0) || (normp%characteristic ==0)))
    {
      for (i=0; i<n1ds; i++) apv[i] = -999;
#ifdef DEBUG_APVEC
      cout << "ignored prime: ap list = " << apv << endl;
#endif
      return apv;
    }

  if (vp>0) // bad prime
    {
      apv = apvec(AtkinLehnerOp(P,N), 1);

      // int ip = std::find(badprimes.begin(),badprimes.end(),P) - badprimes.begin();
      // long aq;
      // for (i=0; i<n1ds; i++)
      //   {
      //     apv[i] = aq = nflist[i].aqlist[ip];
      //     if(!((aq==1)||(aq==-1)))
      //       {
      //         cout<<"Error: Atkin-Lehner eigenvalue "<<aq<<" for Q="<<P
      //             <<" for form # "<<(i+1)<<" is neither +1 nor -1"<<endl;
      //         cout<<"------------------------"<<endl;
      //         nflist[i].display();
      //         cout<<"------------------------"<<endl;
      //         exit(1);
      //       }
      //   }
#ifdef DEBUG_APVEC
      cout<<"Bad prime: aq list = "<<apv<<endl;
#endif
      return apv;
    }

  // now P is a good prime

  long maxap = max_T_P_eigenvalue(P);

  if (Quad::is_Euclidean)
    apv = apvec_euclidean(P, maxap);
  else
    apv = apvec(HeckeOp(P,N), maxap);
#ifdef DEBUG_APVEC
  cout<<"Good prime: ap list = "<<apv<<endl;
#endif
  return apv;
}

// Strategy for operators used to automatically cut out 1-dimensional
// rational eigenspaces, in characteristic 0:

// The first n2r (=Quad::class_group_2_rank) operators are T_{A,A}
// where A runs over n2r ideals coprime to N whose classes generate
// the 2-torsion in the class group. These are involutions, but we
// only consider the +1-eigenspace since we only want to find
// eigenspaces with trivial unramified character.  When the class
// number is odd, then n2r=0 so there are none of these.  The list of
// ideals A used here is newforms::nulist.

// Next come nwq operators T_{A,A}*W_Q, where Q is the i'th prime
// dividing the level to power e and [Q^e] is square (so either [W] is
// square or e is even).  By the usual abuse of notation we write W_Q
// when we really mean W_{Q^e}, and A is coprime to N such that A^2Q^e
// is principal.  The list of primes Q is in badprimes, of length nwq.
// When the class number is odd, nwq is the number of prime factors of
// N, but in general it may be smaller, even 0.  For each of these we
// consider +1,-1 as eigenvalues.

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

matop newforms::h1matop(int i) // return the list of matrices defining the i'th operator
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

long max_T_P_eigenvalue(Quadprime& P)
{
  long normp = I2long(P.norm());
  long aplim=2;
  while (aplim*aplim<=4*normp) aplim++;
  aplim--;
  return aplim;
}

// Return list of integers between -2*sqrt(N(P)) and +2*sqrt(N(P))
vector<long> good_eigrange(Quadprime& P)
{
  long normp = I2long(P.norm());
  long aplim=max_T_P_eigenvalue(P);
  if (P.has_square_class())
    {
      return range(-aplim, aplim);
    }
  else // want eigs of T_{P^2} such that T_P has integral eig
    {
      vector<long> ans = range(0,aplim);
      for (vector<long>::iterator ai = ans.begin(); ai!=ans.end(); ++ai)
        *ai = (*ai)*(*ai) + normp;
      return ans;
    }
}

// the list of possible (integer) eigenvalues for the i'th operator:
vector<long> newforms::eigrange(int i)
{
  vector<long> ans;
  assert (i>=0);

  if (i<n2r)
    {
      ans = {1};
      return ans;
    }
  i -= n2r;

  if (i<nwq)
    {
      if (characteristic==2)
        ans = {1};
      else
        ans = {-1, 1};
      if (verbose>1)
        cout << "eigrange for Q = " << badprimes[i] << " (norm "<<badprimes[i].norm()<<"):\t" << ans << endl;
      return ans;
    }
  i -= nwq;

  if (characteristic>0)
    {
      ans = range(0,characteristic-1);
      if (verbose>1)
        cout << ans << endl;
      return ans;
    }

  ans = good_eigrange(goodprimes[i]);
  if (verbose>1)
    cout << "eigrange for P = " << goodprimes[i] << " (norm "<<goodprimes[i].norm()<<"):\t" << ans << endl;
  return ans;
}

// List of bad primes (dividing N) followed by good primes to length
// at least np, making sure that the list includes at least one good
// principal prime.  iP0 is set to the index in the list of the first
// good principal prime.
vector<Quadprime> make_primelist(Qideal& N, int np, int& iP0, int p)
{
  vector<Quadprime> ans;
  if (p==0)
    ans = N.factorization().sorted_primes();
  vector<Quadprime>::const_iterator Pi = Quadprimes::list.begin();
  iP0 = -1;
  while ((ans.size()<(unsigned)np) || (iP0<0))
    {
      Quadprime P = *Pi++;
      if (P.divides(N))
        continue;
      if ((p>0) && (P.norm()%p==0))
        continue;
      ans.push_back(P);
      if (P.is_principal() && iP0==-1)
        {
          iP0 = ans.size()-1; // index of current last in list, indexed from 0
        }
    }
  // cout<<"make_primelist(N="<<N<<", np="<<np<<") returns "<<ans<<" (size "<<ans.size()<<") with P0="<<ans[iP0]<<", iP0="<<iP0<<endl;
  return ans;
}


// compute a list of ideals coprime to N whose classes generate the 2-torsion
vector<Qideal> make_nulist(Qideal& N)
{
  vector<Qideal> nulist;
  for (vector<Qideal>::iterator Ai = Quad::class_group_2_torsion_gens.begin();
       Ai!=Quad::class_group_2_torsion_gens.end(); ++Ai)
    nulist.push_back(Ai->equivalent_coprime_to(N));
  return nulist;
}

// compute a list of primes Q dividing N with Q^e||N such that [Q^e] is square
vector<Quadprime> make_badprimes(Qideal& N, const vector<Quadprime>& allbadprimes)
{
  vector<Quadprime> badprimes;
  for (vector<Quadprime>::const_iterator Qi = allbadprimes.begin(); Qi!=allbadprimes.end(); ++Qi)
    {
      Quadprime Q = *Qi;
      long e = val(Q,N);
      if (Quad::class_group_2_rank==0 || e%2==0 || Q.has_square_class()) // then [Q^e] is a square
        badprimes.push_back(Q);
    }
  return badprimes;
}

// compute a list of at least nap good primes (excluding those
// dividing characteristic if >0), to include at least on principal
// one which has index iP0;
vector<Quadprime> make_goodprimes(Qideal& N,  int np, int& iP0, int p)
{
  vector<Quadprime> goodprimes;
  QuadprimeLooper L(p==0? N : p*N);
  iP0=-1;
  for (int i=0; (i<np) || (iP0<0); i++, ++L)
    {
      Quadprime P = L;
      goodprimes.push_back(P);
      if (P.is_principal() && iP0==-1)
        iP0 = goodprimes.size()-1;
    }
  return goodprimes;
}

//end of newforms.cc
