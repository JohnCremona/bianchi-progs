//  newforms.cc

#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <functional>   // std::multiplies
#include <numeric>   // std::multiplies
#include "looper.h"
#include "oldforms.h"
#include "newforms.h"
#include "eclib/curvesort.h" // for letter codes

scalar dotmodp(const vec& v1, const vec& v2, scalar pr)
{
  scalar ans=0;
  for(long i=1; i<=dim(v1); i++) ans=xmod(ans+xmodmul(v1[i],v2[i],pr),pr);
  return mod(ans,pr);
}

newform::newform(newforms* nfs, const vec& v, const vector<long>& eigs)
  :basis(v), eigs(eigs)
{
  //cout<<"Constructing newform with eigs "<<eigs<<endl;
  nf=nfs;
  copy(eigs.begin(), eigs.begin()+(nf->npdivs), back_inserter(aqlist));
  dp0    =  1 + quadnorm(nf->p0) - eigs[nf->npdivs];
  if(nf->hmod)
    pdot = dotmodp((nf->mvp),basis,nf->hmod);
  else
    pdot = (nf->mvp)*basis;
  //cout<<"dp0 = "<<dp0<<", pdot="<<pdot<<endl;
  //No division by h1denom() here:
  loverp = rational(abs(pdot),dp0*(Quad::nunits));
  sfe = std::accumulate(aqlist.begin(),aqlist.end(),-1,std::multiplies<long>());
  // Find the ratio of the least period w.r.t. integral homology
  // divided by the least period w.r.t. homology relative to cusps.
  find_cuspidal_factor();
  find_matrix();
}

newform::newform(newforms* nfs,
                 const vector<int>& intdata, const vector<Quad>& Quaddata,
                 const vector<long>& aq, const vector<long>& ap)
{
  nf=nfs;
  sfe = intdata[0];
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
  aplist = ap;
  // recreate eigs list (in case we need to recover basis vector):
  // start with Atkin-Lehner eigs aq, then Tp eigs ap for good p
  eigs = aq;
  vector<Quad>::const_iterator pr=quadprimes.begin();
  vector<long>::const_iterator api=aplist.begin();
  while ((eigs.size()<20) && (api!=aplist.end()))
    {
      while (div(*pr,nf->modulus))  { pr++; api++; }
      eigs.push_back(*api++);
      pr++;
    }
}

void newform::find_cuspidal_factor(void)
{
  vec bc;
  int verbose = nf->verbose;
  cuspidalfactor=1;

  if(!(nf->h1->cuspidal))
    {
      bc=(nf->h1->tkernbas)*basis;
      cuspidalfactor = vecgcd(bc);
      bc /= cuspidalfactor;
    }
  if(verbose&&(cuspidalfactor>1))
    {
      cout<<"cuspidalfactor = "<<cuspidalfactor<<endl;
      if(verbose>2) cout<<"bc = "<<bc<<endl;
    }
}

void newform::find_matrix()
{
  int verbose=(nf->verbose);
  if(verbose) cout<<"computing a,b,c,d..."<<flush;
  Quad N = nf->modulus;
  matdot=0;
   //  Look for a QuadRational q=b/d for which {0,q}={0,g(0)} is
   //  nontrivial, and record the matrix g and matdot = the multiple
   //  of the fundamental period which the integral over {0,g(0)} is.
  for (Quadlooper dl(Quad::d, 2, 1000, 1); dl.ok()&&!matdot; ++dl)
    { d=(Quad)dl;
      if (coprime(d,N))
        {
          //cout<<"d="<<d<<endl;
          vector<Quad> reslist = residues(d);
          vector<Quad>::const_iterator res;
          for(res=reslist.begin(); res!=reslist.end() && !matdot; res++)
            {
              b=*res;
              Quad g = quadbezout(N*b,d,c,a);  // b*N*c+a*d=1
              if (quadnorm(g)==1)
                {   // found a candidate q=b/d
                  //cout<<"b="<<b<<endl;
                  vec v=nf->h1->chain(b,d);  //starts at 1
                  //cout<<"v="<<v<<endl;
                  if(nf->hmod)
                    matdot = dotmodp(v,basis,nf->hmod);
                  else
                    matdot = v*basis;
                  //cout<<"matdot="<<matdot<<endl;
                  if (matdot)
                    {
                      c = -N*c;
                      if (!(a*d-b*c==1))
                        cout<<"Error computing matrix"<<endl;
                    }
                } // b coprime to d
            } // loop over b
        } // d coprime to modulus
    } // loop over d
}

int newform::is_base_change(void) const
{
  if(!(nf->is_Galois_stable))
    return 0;
  vector<long>::const_iterator ap = aplist.begin();
  vector<Quad>::const_iterator pr=quadprimes.begin();
  while(ap!=aplist.end())
    {
      long api = *ap++;
      Quad p0 = *pr++;
      //cout<<"p="<<p0<<" has ap="<<api<<endl;
      if(!is_ideal_Galois_stable(p0)) // this prime not inert or ramified
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
          pr++;
        }
    }
  //cout<<"All OK -- base-change"<<endl;
  return 1;
}

int newform::is_base_change_twist(void) const
{
  vector<long>::const_iterator ap = aplist.begin();
  vector<Quad>::const_iterator pr=quadprimes.begin();
  while(ap!=aplist.end())
    {
      long api = *ap++;
      Quad p0 = *pr++;
      //cout<<"p="<<p0<<" has ap="<<api<<endl;
      if(!is_ideal_Galois_stable(p0)) // this prime not inert or ramified
        {
          if (ap==aplist.end()) // the conjugate ap is not known
            {
              //cout<<"All OK -- base-change up to twist"<<endl;
              return 1;
            }
            // read next (conjugate) prime and eigenvalue:
          long apj = *ap++;
          Quad p1 =  *pr++;
          //cout<<"p'="<<p1<<" has ap="<<apj<<endl;
          // skip if either divides level:
          if(div(p0,nf->modulus) || div(p1,nf->modulus))
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

long squarefree_part(long d)
{
  if (d==0) return d;
  vector<long> sd = sqdivs(d);
  long maxd = sd[sd.size()-1];
  long ans = d/(maxd*maxd);
  //cout << "d has max square divisor "<<maxd<<"^2"<<" and squarefree part "<<ans<<endl;
  return ans;
}

// if form is base-change, find the d s.t. the bc has eigenvalues in Q(sqrt(d))
int newform::base_change_discriminant(void) const
{
  if (is_base_change()==0) return 0;
  int bcd = 1;
  long ap, dp;
  Quad p;
  vector<long>::const_iterator api = aplist.begin();
  vector<Quad>::const_iterator pr=quadprimes.begin();
  while(api!=aplist.end())
    {
      ap = *api++;
      p = *pr++;
      if(imag(p)!=0) // this prime is not inert
        continue;
      if(div(p,nf->modulus)) // this prime is bad
        continue;
      dp = ap+2*real(p);
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
  long ap, dp1, dp2;
  Quad p;
  vector<long>::const_iterator api = aplist.begin();
  vector<Quad>::const_iterator pr=quadprimes.begin();
  while(api!=aplist.end())
    {
      ap = *api++;
      p = *pr++;
      if(imag(p)!=0) // this prime is not inert
        continue;
      if(div(p,nf->modulus)) // this prime is bad
        continue;
      dp1 =  ap+2*real(p);
      dp2 = -ap+2*real(p);
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
  long ap, dp;
  Quad p;
  vector<long>::const_iterator api = aplist.begin();
  vector<Quad>::const_iterator pr=quadprimes.begin();
  while(api!=aplist.end())
    {
      ap = *api++;
      p = *pr++;
      if (ap==0) continue;
      dp = ap*ap-4*quadnorm(p);
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
      h1 = new homspace(modulus,1,0);
      nfhmod=hmod = h1->h1hmod();
    }
}

long newforms::matdim(void) {return h1->dimension;}
long newforms::matden(void) {return h1->denom3;}

newforms::newforms(const Quad& n, int disp)
 :level(n), verbose(disp)
{
  nap=level::nap;
  nwq=level::npdivs;
  ntp=nap-nwq;
  h1=0;
  of=0;
  nfhmod=0;
//
// get smallest "good" prime:
//
  vector<Quad>::const_iterator pr=quadprimes.begin();
  p0 = *pr;
  while (div(p0,modulus)) p0=*pr++;     // First "good" prime
}

void newforms::get_lambda()
{
//#define DEBUG_LAMBDA
  int* gotlambda = new int[n1ds];
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
        nflist[i].lambda=Quad(1);
        nflist[i].lambdadot=nflist[i].pdot;
        gotlambda[i]=1;
        nfound++;
      }
    else
      {
        nflist[i].lambda=Quad(0); // indicates no lambda exists
        nflist[i].lambdadot=0;
	gotlambda[i]=0;
      }
#ifdef DEBUG_LAMBDA
  if(verbose)cout<<nfound<<" easy cases out of "<<n1ds<<endl;
#endif
  if (is_square)
    {
      delete[] gotlambda;
      return;
    }

  vector<Quad>::const_iterator pr;
  for(pr=quadprimes.begin(); pr!=quadprimes.end() && (nfound<n1ds); pr++)
    { Quad lam = *pr;
#ifdef DEBUG_LAMBDA
      if(verbose)cout << "Testing lambda = " << lam << endl;
#endif
      if(!div(lam,2*modulus))
        {
#ifdef DEBUG_LAMBDA
          if(verbose)cout<<"passed first general test"<<endl;
#endif
          vector<Quad> lamres = residues(lam);
          if(squaremod(fundunit,lam,lamres)==1)
            {
#ifdef DEBUG_LAMBDA
              if(verbose)cout<<"passed second general test"<<endl;
#endif
              int chimod  = squaremod(modulus,lam,lamres);
              int* chitab = makechitable(lam,lamres);
              vec mvtw = h1->manintwist(lam,lamres,chitab);
              delete[] chitab;
              for(int j=0; (j<n1ds)&&(nfound<n1ds); j++)
                {
                  if(gotlambda[j]==0)
                    {
#ifdef DEBUG_LAMBDA
                      if(verbose)cout<<"Newform # "<<j<<": ";
                      if(verbose)cout<<"trying: ";
#endif
                      newform& nfj = nflist[j];
                      int dot;
                      if(hmod)
                        dot = abs(dotmodp(mvtw,nfj.basis,hmod));
                      else
                        dot = abs(mvtw*nfj.basis);
                      if(dot&&((chimod*nfj.sfe)==+1))
                        {
#ifdef DEBUG_LAMBDA
                          if(verbose)cout<<"Success! ";
#endif
                          nfj.loverp = rational(dot,(Quad::nunits));
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
  delete[] gotlambda;
}

void newforms::createfromscratch()
{
  if(verbose)
    cout<<"Constructing homspace...\n";
  makeh1plus();
  nfhmod=hmod = h1->h1hmod();
  mvp=h1->maninvector(p0);
  int dimcusp = h1->h1cuspdim();
  int dimall = h1->h1dim();
  if(verbose)
      cout<<"Dimension = "<<dimall<<" (cuspidal dimension = "<<dimcusp<<")\n";

  if(verbose)
    cout<<"Retrieving oldform data...\n";
  of = new oldforms(this,verbose>1);
  if(verbose)
    of->display();

  if(verbose)
    cout<<"Finding rational newforms...\n";
  maxdepth = nap;
  long mindepth = npdivs;
  dimsplit = n1ds = 0;
  long olddimall = (of->olddimall);
  if(verbose>1) cout<<"olddimall = "<<olddimall<<endl;
  nnflist = upperbound = (h1->h1cuspdim()) - olddimall;
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

  if(n1ds==0) return; // else no work to do

// Find the twisting primes for each newform (more efficient to do
// this here instead of within the newform constructors, as one lambda
// might work for more than one newform). NB If the level is square
// SFE=-1 then no such lambda will exist.

  get_lambda();
  allproj();    // Compute homspace::projcoord, so projcycle can be used
  find_jlist();
  nap=0;
}

void newforms::use(const vec& b1, const vec& b2, const vector<long> eigs)
{
  if (use_nf_number==-1)
    {
      if (n1ds<upperbound)
        {
          nflist.push_back(newform(this,b1,eigs));
          n1ds++;
        }
      else
        {
          cout << "Error in splitting eigenspaces: apparently found more ";
          cout << "1D newforms ("<< n1ds+1 <<") than the total new-dimension ("
               <<upperbound<<").\n";
          exit(1);
        }
    }
  else // store eigs and basis
    {
      nflist[use_nf_number].eigs = eigs;
      nflist[use_nf_number].basis = b1;
    }
}

void newforms::display(void) const
{
 if (n1ds==0) {cout<<"No newforms."<<endl; return;}
 cout << "\n"<<n1ds<<" newform(s) at level " << ideal_label(modulus) << " = (" << modulus << "):" << endl;
 for(int i=0; i<n1ds; i++)
   {cout<<i+1<<":\t";
    nflist[i].display();
  }
}

void newform::display(void) const
{
 cout << "basis = " << basis 
   //      << ";\teigs = " << eigs 
      << ";\taqlist = " << aqlist 
      << ";\taplist = " << aplist << endl;
 cout << "Sign of F.E. = " << sfe << endl;
 cout << "Twisting prime lambda = " << lambda << ", factor = " << lambdadot << endl;
 cout << "L/P ratio    = " << loverp << ", cuspidal factor = " << cuspidalfactor << endl;
 cout << "Integration matrix = [" << a << "," << b << ";" << c << "," << d << "], factor   = " << matdot << endl;
}

void newforms::list(long nap) const
{
  string id = ideal_label(modulus);
  string flabel = field_label();
  int bc, bcd, bct;
  for(int i=0; i<n1ds; i++)
    {
      cout << flabel << " " << id << " " << codeletter(i) << " (" << modulus <<") ";
      // weight
      cout << "2 ";
      bc = nflist[i].is_base_change();
      if (bc)
        {
          bcd = nflist[i].base_change_discriminant();
          bct = 0;
          cout << bcd;
        }
      else
        {
          bcd = 0;
          bct = nflist[i].is_base_change_twist();
          if (bct)
            {
              // NB if we have not enough inert a_P we might not be
              // able to determine the discriminant; the following
              // will return 1 in this case.
              bcd = nflist[i].base_change_twist_discriminant();
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
  for(ai=aqlist.begin(); ai!=aqlist.end(); ai++)
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
  for(ai=aplist.begin(); ai!=aplist.begin()+nap && ai!=aplist.end(); ai++)
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

struct less_newform_eigs : public binary_function<newform, newform, bool> {
  bool operator()(const newform& f, const newform& g)
  {
    return less_apvec(f.eigs,g.eigs)==-1;
  }
};

struct less_newform_lmfdb : public binary_function<newform, newform, bool> {
  bool operator()(const newform& f, const newform& g)
  {
    return less_apvec(f.aplist,g.aplist)==-1;
  }
};

// struct less_apvec_function : public binary_function<const vector<long>&, const vector<long>&, bool> {
//   bool operator()(const vector<long>& f, const vector<long>& g)
//   {
//     return 1==less_apvec(f,g);
//   }
// };

void newforms::sort_eigs(void)
{
  ::sort(nflist.begin(),nflist.end(),less_newform_eigs());
}

void newforms::sort_lmfdb(void)
{
  ::sort(nflist.begin(),nflist.end(),less_newform_lmfdb());
}


// vec newforms::proj(const vec& v)  //returns vec of components in each eig-space
// {
//   vec ans(n1ds);
//   for (int i=0; i<n1ds; i++) ans[i+1]=(nflist[i].basis)*v;
//   return ans;
// }

void newforms::allproj() //Replaces "coord" member of homspace with projections
                      //onto eigenspaces, to save time
{
  int ncoord = (h1->coord).nrows(); long pcij;
  h1->projcoord.init(ncoord,n1ds);
  for (int i=1; i<=ncoord; i++)
    {
      vec coordi = h1->coord.row(i);
      for (int j=1; j<=n1ds; j++)
        {
	  if (hmod) pcij = dotmodp(coordi,nflist[j-1].basis, hmod);
          else      pcij = coordi * (nflist[j-1].basis);
          h1->projcoord.set(i,j, pcij);
        }
    }
  if(!hmod) return;

  // Check lifts of projcoord columns
  vec c1, c2(ncoord);
  int test_prim = 0;
  for (int j=1; j<=n1ds; j++)
    {
      c1 = h1->projcoord.col(j);
      if (lift(c1,hmod,c2))
        {
          if((c1==c2)||(c1==-c2))
            {
              if(test_prim) 
		cout<<"projcoord column "<<j<<" is already Z-primitive"<<endl;
            }
          else
            {
              if(test_prim) 
		cout<<"projcoord column "<<j<<"="<<c1<<" lifts ok to Z-primitive "<<c2<<endl;
              h1->projcoord.setcol(j,c2);
            }
          nfhmod=0;
        }
      else
        {
	  if(test_prim) 
	    cout<<"projcoord column "<<j<<" cannot be lifted to Z"<<endl;
        }
    }
}

// Second constructor, in which data is read from file in directory
// newforms. If not, recreates it from eigs

void newforms::createfromdata()
{
  if(verbose) cout << "Retrieving newform data for N = " << modulus << endl;

// Read newform data from file into eigdata structure.

  eigdata filedata(this,modulus,-1,verbose>1);  // neigs=-1 means get ALL from file

// Extract number of newforms and their eigenvalues from this.

  nnflist=n1ds=filedata.nforms;
  n2ds=filedata.nforms2;

 // construct the newforms from this data
  for(int i=0; i<n1ds; i++)
    nflist.push_back(newform(this,filedata.intdata[i],filedata.Quaddata[i],
                                  filedata.aqs[i],filedata.aps[i]));
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
  if(verbose>1) cout<<"Doing allproj() and find_jlist()..."<<endl;
  allproj();
  find_jlist();
  if(verbose) cout<<"Finished makebases()"<<endl;
}

void newforms::getoneap(const Quad& p, int verbose, int store)
{
  vector<long> apv=apvec(p);
  int vp = val(p,modulus), ap, cp, i;

  if(verbose)
    {
      if(vp>0) cout<<"q"; else cout<<"p";
      cout<<" = "<<p<<"\t";
    }
  for (i=0; i<n1ds; i++)
    {
      ap = apv[i];
      cp = (vp==0?ap:(vp==1?-ap:0));
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
  if(last>nquadprimes)
    {
      last=nquadprimes;
      cout<<"Cannot compute more than "<<nquadprimes
          <<" ap since we only have that many primes precomputed"<<endl;
    }
  if(last<=nap)
    {
      cout<<"Already have "<<nquadprimes <<" ap so no need to compute more"<<endl;
    }
  // now nap < last <= nquadprimes
  vector<Quad>::const_iterator pr = quadprimes.begin()+first-1;
  while((pr!=quadprimes.end()) && (nap<last))
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
  for(f=nflist.begin(); f!=nflist.end(); f++)
    {
      out<<setw(5)<<(f->sfe)<<" ";
      if(echo) cout<<setw(5)<<(f->sfe)<<" ";
    }
  out<<endl;  if(echo) cout<<endl;
  // Line 4: pdots
  for(f=nflist.begin(); f!=nflist.end(); f++)
    {
      out<<setw(5)<<(f->pdot)<<" ";
      if(echo) cout<<setw(5)<<(f->pdot)<<" ";
    }
  out<<endl;  if(echo) cout<<endl;
  // Line 5: dp0s
  for(f=nflist.begin(); f!=nflist.end(); f++)
    {
      out<<setw(5)<<(f->dp0)<<" ";
      if(echo) cout<<setw(5)<<(f->dp0)<<" ";
    }
  out<<endl;  if(echo) cout<<endl;
  // Line 6: cuspidal factors
  for(f=nflist.begin(); f!=nflist.end(); f++)
    {
      out<<setw(5)<<(f->cuspidalfactor)<<" ";
      if(echo) cout<<setw(5)<<(f->cuspidalfactor)<<" ";
    }
  out<<endl;  if(echo) cout<<endl;
  // Line 7: lambdas
  for(f=nflist.begin(); f!=nflist.end(); f++)
    {
      Quad lambda = f->lambda;
      out<<setw(5)<<real(lambda)<<" "<<imag(lambda)<<" ";
      if(echo) cout<<setw(5)<<real(lambda)<<" "<<imag(lambda)<<" ";
    }
  out<<endl;  if(echo) cout<<endl;
  // Line 8: lambdadots
  for(f=nflist.begin(); f!=nflist.end(); f++)
    {
      out<<setw(5)<<(f->lambdadot)<<" ";
      if(echo) cout<<setw(5)<<(f->lambdadot)<<" ";
    }
  out<<endl;  if(echo) cout<<endl;
  // Lines 9,10,11,12: a,b,c,d:
  Quad a;
  for(f=nflist.begin(); f!=nflist.end(); f++)
    {
      a = f->a;
      out<<setw(5)<<real(a)<<" "<<imag(a)<<" ";
      if(echo) cout<<setw(5)<<real(a)<<" "<<imag(a)<<" ";
    }
  out<<endl;  if(echo) cout<<endl;
  for(f=nflist.begin(); f!=nflist.end(); f++)
    {
      a = f->b;
      out<<setw(5)<<real(a)<<" "<<imag(a)<<" ";
      if(echo) cout<<setw(5)<<real(a)<<" "<<imag(a)<<" ";
    }
  out<<endl;  if(echo) cout<<endl;
  for(f=nflist.begin(); f!=nflist.end(); f++)
    {
      a = f->c;
      out<<setw(5)<<real(a)<<" "<<imag(a)<<" ";
      if(echo) cout<<setw(5)<<real(a)<<" "<<imag(a)<<" ";
    }
  out<<endl;  if(echo) cout<<endl;
  for(f=nflist.begin(); f!=nflist.end(); f++)
    {
      a = f->d;
      out<<setw(5)<<real(a)<<" "<<imag(a)<<" ";
      if(echo) cout<<setw(5)<<real(a)<<" "<<imag(a)<<" ";
    }
  out<<endl;  if(echo) cout<<endl;
  // Line 13: matdots
  for(f=nflist.begin(); f!=nflist.end(); f++)
    {
      out<<setw(5)<<(f->matdot)<<" ";
      if(echo) cout<<setw(5)<<(f->matdot)<<" ";
    }
  out<<endl;  if(echo) cout<<endl;
  out<<endl;  if(echo) cout<<endl;
  for(int i=0; i<nwq; i++)
    {
      for(f=nflist.begin(); f!=nflist.end(); f++)
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
      for(f=nflist.begin(); f!=nflist.end(); f++)
	{
	  out<<setw(5)<<(f->aplist)[i]<<" ";
	  if(echo) cout<<setw(5)<<(f->aplist)[i]<<" ";
	}
      out<<endl;      if(echo) cout<<endl;
    }
  out.close();
}

void newforms::find_jlist()
{
  int i, j, ok=0; j0=0;
  for(j=1; (!ok)&&(j<=h1->h1dim()); j++)
    {
      ok=1;
      for (i=0; (i<n1ds)&&ok; i++)
        ok=(nflist[i].basis[j]!=0);
      if(ok) j0=j;
    }
  if(ok)
    {
      if(verbose)  cout<<"j0="<<j0<<endl;
      jlist.insert(j0);
      for (i=0; i<n1ds; i++)
	{
	  nflist[i].j0 = j0;
          nflist[i].fac = nflist[i].basis[j0];
          if(nfhmod) nflist[i].facinv=invmod(nflist[i].fac,nfhmod);
	}
    }
  else
    {
      if(verbose)
	cout<<"Failed to find j0 such that nflist[i].basis[j0]!=0 for all i"
	    <<endl;
      // Find out which pivots we'll be using:
      for (i=0; i<n1ds; i++)
	{
	  vec& bas = nflist[i].basis;
	  j=1; while(bas[j]==0) j++;
	  jlist.insert(j);
	  nflist[i].j0 = j;
	  nflist[i].fac = nflist[i].basis[j];
	}
      if(verbose)  cout<<"jlist="<<jlist<<endl;
    }
}

//#define DEBUG_APVEC

void update(const mat& pcd, vec& imagej, long ind, long hmod)
{
  vec part;
  if(ind>0)
    part = pcd.row(ind);
  else
    part = -pcd.row(-ind);
#ifdef DEBUG_APVEC
  cout<<"--adding "<<part<<" (ind="<<ind<<") to imagej, ";
#endif
  if(hmod)
    {
      imagej.addmodp(part,hmod);
      imagej=reduce_modp(imagej);
    }
  else
    imagej+=part;
#ifdef DEBUG_APVEC
  cout<<"updated imagej is "<<imagej<<endl;
#endif
}

vector<long> newforms::apvec(const Quad& p)  // computes a[p] for each newform
{
#ifdef DEBUG_APVEC
  cout<<"In apvec with p = "<<p<<endl;
#endif
  vector<long> apv(n1ds);
  vec v;
  long i,j,ap,aq,normp=quadnorm(p);

  int vp = val(p,modulus);

  if (vp>0) // bad prime, we already know the eigenvalues
    {
      for (i=0; i<n1ds; i++)
        {
          int ip = find(plist.begin(),plist.end(),p)-plist.begin();
          aq = nflist[i].aqlist[ip];
          if(!((aq==1)||(aq==-1)))
            {
              cout<<"Error: Atkin-Lehner eigenvalue "<<aq<<" for q="<<p
                  <<" for form # "<<(i+1)<<" is neither +1 nor -1"<<endl;
              cout<<"------------------------"<<endl;
              nflist[i].display();
              cout<<"------------------------"<<endl;
            }
          apv[i] = aq;
          // if(vp==1)
          //   apv[i] = -aq;
          // else
          //   apv[i] = 0;
        }
      return apv;
    }

  // now p is a good prime

  long maxap=(long)(2*sqrt((double)normp)); // for validity check

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

  for(std::set<long>::const_iterator jj=jlist.begin(); jj!=jlist.end(); jj++)
    {
      vec imagej=vec(n1ds); // initialised to 0
      j=*jj;
      symb s = h1->symbol(h1->freegens[j-1]);
#ifdef DEBUG_APVEC
      cout<<"Computing image under T("<<p<<") of "<<j<<"'th M-symbol"
          <<" = "<<s<<"..."<<flush;
#endif
      Quad u=s.cee(),v=s.dee();
      mat& pcd = h1->projcoord;
      //cout<<"projcoord = "<<pcd;
// Matrix [1,0;0,p]
      long ind = h1->coordindex[h1->index2(u,p*v)];
#ifdef DEBUG_APVEC
      cout<<"u1="<<u<<", u2="<<p*v<<", ind="<<ind<<endl;
#endif
      if(ind) update(pcd,imagej,ind,nfhmod);
// Matrix [p,0;0,1]
      ind = h1->coordindex[h1->index2(p*u,v)];
#ifdef DEBUG_APVEC
      cout<<"u1="<<p*u<<", u2="<<v<<", ind="<<ind<<endl;
#endif
      if(ind) update(pcd,imagej,ind,nfhmod);
// Other matrices
      vector<Quad> resmodp=residues(p);
      vector<Quad>::const_iterator res=resmodp.begin();
      while(res!=resmodp.end())
        {
          b = *res++;
          if(b==0) continue; // handled above as special case
          a = -p;
          u1=u*p; u2=v-u*b;
          ind = h1->coordindex[h1->index2(u1,u2)];
#ifdef DEBUG_APVEC
          cout<<"Residue class "<<b<<": ";
          cout<<"a="<<a<<", b="<<b<<", u1="<<u1<<", u2="<<u2<<", ind="<<ind<<endl;
#endif
          if(ind) update(pcd,imagej,ind,nfhmod);
          while(b!=0)
            {
              q=a/b; c=a-b*q; u3=q*u2-u1;
              a=-b; b=c; u1=u2; u2=u3;
              ind = h1->coordindex[h1->index2(u1,u2)];
#ifdef DEBUG_APVEC
              cout<<"a="<<a<<", b="<<b<<", u1="<<u1<<", u2="<<u2<<", ind="<<ind<<endl;
#endif
              if(ind) update(pcd,imagej,ind,nfhmod);
            }
#ifdef DEBUG_APVEC
          cout<<" image after term is "<<imagej<<endl;
#endif
        }
      images[j]=imagej/(h1->h1denom());
#ifdef DEBUG_APVEC
      cout<<" image after scaling is "<<images[j]<<endl;
#endif
    }

  for (i=0; i<n1ds; i++)
    {
// recover eigenvalue:
#ifdef DEBUG_APVEC
      cout << "numer = " << images[nflist[i].j0][i+1] << endl;
      cout << "denom = " << nflist[i].fac << endl;
#endif
      if(nfhmod)
        ap=xmodmul(images[nflist[i].j0][i+1],nflist[i].facinv,nfhmod);
      else
        ap=images[nflist[i].j0][i+1]/nflist[i].fac;
      apv[i]=ap;
#ifdef DEBUG_APVEC
      cout << "ap = " << ap << endl;
#endif
// check it is in range:
      if((ap>maxap)||(-ap>maxap))
	{
	  cout<<"Error:  eigenvalue "<<ap<<" for p="<<p
	      <<" for form # "<<(i+1)<<" is outside valid range "
	      <<-maxap<<"..."<<maxap<<endl;
          //          exit(1);
	}
    }
  return apv;
}

//end of newforms.cc
