//  newforms.cc

#include <iostream>
#include <sstream>
#include "oldforms.h"
#include "newforms.h"

newform::newform(const newforms* nfs, const vec& v, const vector<long>& ap)
 :basis(v),aplist(ap) 
{
  dp0    =  1 + quadnorm(nfs->p0) - aplist[nfs->npdivs];
  pdot   =  (nfs->mvp)*basis;
  loverp = rational(abs(pdot),dp0*(Quad::nunits));  //No division by h1denom()
  int iq;
  for(iq=0, sfe=-1; iq<nfs->npdivs; iq++) sfe*=aplist[iq];
}

void newforms::makeh1plus(void)
{
  if(!h1) h1 = new homspace(modulus,1,0);
}

newforms::newforms(const Quad& n, int useolddata, int disp)
 :level(n), verbose(disp)
{
  nap=level::nap;
  nwq=level::npdivs;
  ntp=nap-nwq;
  h1=0; 
  makeh1plus();   // necessary until we implement newforms file data
  hmod = h1->h1hmod();
  of=0;
//
// get smallest "good" prime and manin-vector:
//
  vector<Quad>::const_iterator pr=quadprimes.begin();
  p0 = *pr; 
  while (div(p0,modulus)) p0=*pr++;     // First "good" prime
  mvp=h1->maninvector(p0);
  if(verbose) 
    {
      cout<<"Constructing newforms\n";
      cout<<"p0 = "<<p0<<", mvp = "<<mvp<<"\n";
    }

  if(useolddata)createfromeigs(); 
  else createfromscratch();

// Finally find the twisting primes for each newform (more efficient to do 
// this here instead of within the newform constructors, as one lambda might
// work for more than one newform).
//
  if(n1ds==0) return; // no work to do

//#define DEBUG_LAMBDA

  int* gotlambda = new int[n1ds];
  int i, nfound=0;
  // first check for "easy" forms where no lambda needed:
#ifdef DEBUG_LAMBDA
  if(verbose)
  {
   cout<<"Newform data so far:\n";
   for(int k=0; k<n1ds; k++) nflist[k].display();
   cout<<"Now looking for twisting primes.\n";
  }
#endif
  for (i=0; i<n1ds; i++) 
    if(nflist[i].pdot!=0) 
      {
#ifdef DEBUG_LAMBDA
        if(verbose)cout<<"Newform "<<i<<": lambda=1 will do.\n";
#endif
        nflist[i].lambda=Quad(1);
        gotlambda[i]=1;
        nfound++;
      }
    else gotlambda[i]=0;
#ifdef DEBUG_LAMBDA
  if(verbose)cout<<nfound<<" easy cases out of "<<n1ds<<endl;
#endif
  pr=quadprimes.begin();
  for(; pr!=quadprimes.end() && (nfound<n1ds); pr++)
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
                      int dot = abs(mvtw*nfj.basis);
                      if(dot&&((chimod*nfj.sfe)==+1))
                        {
#ifdef DEBUG_LAMBDA
                          if(verbose)cout<<"Success! ";
#endif
                          nfj.loverp = rational(dot,(Quad::nunits));
                          nfj.lambda = lam;
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
  of = new oldforms(this,verbose>1);
  if(verbose)of->display();
  maxdepth = nap;
  long mindepth = npdivs;
  dimsplit = n1ds = 0;
  long olddimall = (of->olddimall);
  if(verbose>1) cout<<"olddimall = "<<olddimall<<endl;
  nnflist = upperbound = (h1->h1cuspdim()) - olddimall;
  if(verbose>1) cout<<"upperbound = "<<upperbound<<endl;
  if(upperbound>0)  // Else no newforms certainly so do no work!
    {
       form_finder ff(this,1,maxdepth,mindepth,1,0,verbose);
       ff.find();
     }
  if(verbose>1) cout<<"n1ds = "<<n1ds<<endl;
  n2ds=upperbound-n1ds; // dimension of new, non-rational forms
  if(verbose>1) cout<<"n2ds = "<<n2ds<<endl;
  if(verbose)
    {cout << "Total dimension " << h1->h1dim() << " made up as follows:\n";
     cout << "dim(newforms) = " << n1ds+n2ds << " of which " << n1ds << " is rational; \n";
     cout << "dim(oldforms) = " << olddimall << " of which " << (of->olddim1) << " is rational; \n";
     cout<<endl;
   }
  delete of;
}

void newforms::use(const vec& b1, const vec& b2, const vector<long> aplist)
{
  if (n1ds<upperbound)
    {
      nflist.push_back(newform(this,b1,aplist));
      n1ds++;
    }
  else 
    {  
      cerr << "Error in splitting eigenspaces: apparently found more ";
      cerr << "1D newforms ("<< n1ds+1 <<") than the total new-dimension ("
           <<upperbound<<").\n";
    }
}

void newforms::display(void) const
{
 if (n1ds==0) {cout<<"No newforms."<<endl; return;}
 cout << "\n"<<n1ds<<" newform(s) at level " << modulus << ":" << endl;
 for(int i=0; i<n1ds; i++)
   {cout<<i+1<<":\t";
    nflist[i].display();
  }
}

void newform::display(void) const
{
 cout << "basis = " << basis << ";\taplist = " << aplist << endl;
 cout << "Sign of F.E. = " << sfe << endl;
 cout << "Twisting prime lambda = " << lambda << endl;
 cout << "L/P ratio    = " << loverp << endl;
}

vec newforms::proj(const vec& v)  //returns vec of components in each eig-space
{
  vec ans(n1ds);
  for (int i=0; i<n1ds; i++) ans[i+1]=(nflist[i].basis)*v;
  return ans;
}

void newforms::allproj() //Replaces "coord" member of homspace with projections
                      //onto eigenspaces, to save time
{
  int ncoord = nrows(h1->coord);
  h1->projcoord.init(ncoord,n1ds);
  for (int i=1; i<=ncoord; i++)
    { 
      //      vec coordi = h1->kernelpart(h1->coord.row(i));
      vec coordi = h1->coord.row(i);
      for (int j=1; j<=n1ds; j++)
        { 
	  long pcij = coordi * (nflist[j-1].basis);
	  if (hmod) pcij = mod(pcij,hmod);
          h1->projcoord.set(i,j, pcij);
        }
    }
}

// Second constructor, in which data is read from file in directory newforms/ 
// If not, recreates it from eigs

void newforms::createfromeigs()
{
  if(verbose) cout << "Retrieving newform data for N = " << modulus << endl;
//
//STEP ONE: Read number of newforms and their eigenvalues from file
//
  eigdata filedata(this,modulus,-1,0);  // neigs=-1 means get ALL from file
  nnflist=n1ds=filedata.nforms;
  n2ds=filedata.nforms2;
  if(n1ds==0) return;
//
//STEP 2: Find bases
//We now have the aplists for the newforms (in filedata.eigs), and must 
//find the corresponding basis vectors. 
//

  form_finder splitspace(this, 1, nap, 0, 1, verbose);
  upperbound = n1ds; n1ds=0;  
  splitspace.recover(filedata.eigs);
  if(n1ds!=upperbound)
    {
      cout << "Error: out of " << upperbound << " newform(s) expected, "
           << "only successfully reconstructed " << n1ds << endl;
    }			 

}

//end of newforms.cc
