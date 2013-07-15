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
  of = new oldforms(verbose>1);
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
#ifdef USE_XSPLIT
       form_finder ff(this,1,maxdepth,mindepth,1,0,verbose);
       ff.find();
#else
       tpmats = new mat[ntp];           //Won't be computed until needed
       tpknown= new int[ntp];
       int i=ntp;
       while(i--) tpknown[i]=0;
       wmats = new mat[npdivs]; i=0;          //All computed at start
       vector<Quad>::const_iterator pr;
       for(pr=plist.begin(); pr!=plist.end(); pr++, i++)
         {wmats[i] = transpose(h1->wop(*pr,0));}
       subspace startingspace(h1->dimension);            //The full space
       vector<long> v0;                               //Null list, length 0
       if(verbose)cout << "W matrices computed; calling wsplit." << endl;
       wsplit( startingspace, 0, 0, v0);
       delete[] tpmats; delete[] wmats;
       delete[] tpknown;
#endif
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
          h1->projcoord.set(i,j, coordi * (nflist[j-1].basis));
        }
    }
}

// Second constructor, in which data is read from file in directory newforms/ 
// If not, recreates it from eigs

void newforms::createfromeigs()
{
#ifndef USE_XSPLIT
  cout << "newforms::createfromeigs() not implemented without USE_XSPLIT"<<endl;
  return;
#endif

  if(verbose) cout << "Retrieving newform data for N = " << modulus << endl;
//
//STEP ONE: Read number of newforms and their eigenvalues from file
//
  eigdata filedata(modulus,-1,0);  // neigs=-1 means get ALL from file
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

#ifndef USE_XSPLIT

void newforms::usespace(const SUBSP& s, const vector<long>& aplist)
{
  if (verbose) 
  { 
    cout << "Found a 1D eigenspace! (#" << n1ds+1 << ")" << endl;
    cout << "List of ap: " << aplist << endl;
  }
  vec bas =  getbasis1(&s);
  use(bas,bas,aplist);  // second parameter is dummy
}

void newforms::tsplit(const SUBSP& s,int ip,int depth,const vector<long>& aplist)
{
   if (verbose) 
      cout << "In split, depth = " << depth << ", aplist = " << aplist << "\n";
   int dimsofar = dim(s), dimold = of->dimoldpart(aplist);
   if (verbose) cout<<"dimsofar="<<dimsofar
                    <<", dimold="<<dimold
                    <<", dimnew="<<dimsofar-dimold<<endl;
   if (dimold==dimsofar)
     {dimsplit+=dimsofar;
      if (verbose)
        {cout<<"Abandoning a common eigenspace of dimension "<<dimsofar;
         cout<<" which is a sum of oldclasses."<<endl;
         cout<<"List of ap:"<<aplist<<endl;
       }
      return;   // This branch of the recursion ends: all is old
    }
   if ((dimsofar==1) && (depth>npdivs)) //Want at least one good p!
   { usespace(s,aplist); 
     dimsplit ++;
     return;
   }
   if (depth==maxdepth)
   { if (1) // we want to see THIS message whatever the verbosity level!
     { cout << "\nFound a "<<dimsofar<<"D common eigenspace\n";
       cout << "List of ap: "<<aplist<<"\n";
       cout << "Abandoning, even though oldforms only make up ";
       cout << dimold << "D of this." << endl;
      }
     dimsplit += dimsofar;
     return;
   }
// The recursive part:
   int inewp=ip+1; Quad newp = primelist[npdivs+ip];
   if (! tpknown[ip])
   { tpmats[ip] = transpose(h1->heckeop(newp,0));
     tpknown[ip] = 1;
   }
   int aplim=0;
   while (aplim*aplim<=4*quadnorm(newp)) aplim++; aplim--;
   if (verbose) cout << "Using p = " << newp << ", |ap| up to "<<aplim<<"\n";
   MAT t = RESTRICT(tpmats[ip] , s);

   for (int a=1; (a<= 2*aplim+1) ; a++)
   { long ap = (odd(a) ? (1-a)/2 : a/2);
     SCALAR lambda = h1->h1denom()*ap*denom(s);
     SUBSP temp = EIGENSPACE(t,lambda);
     SUBSP newspace = COMBINE(s,temp);
     int newdim = dim(newspace);
     if (newdim>0)
     { int newdepth = depth+1;
       vector<long> newaplist(newdepth);
       int i;
       for (i=0; i<depth; i++) newaplist[i]=aplist(i);
       newaplist[depth] = ap;
       tsplit(newspace,inewp,newdepth,newaplist);
     }
   }
}

void newforms::wsplit(const SUBSP& s,int iq,int depth,const vector<long>& aplist)
{
   if (verbose) 
     cout << "In wsplit, depth = " << depth << ", aplist = " << aplist <<endl;
   int dimsofar = dim(s), dimold = of->dimoldpart(aplist);
   if (verbose) cout<<"dimsofar="<<dimsofar
                    <<", dimold="<<dimold
                    <<", dimnew="<<dimsofar-dimold<<endl;
   if (dimold==dimsofar)
     {dimsplit+=dimsofar;
      if (verbose)
        {cout<<"Abandoning a common eigenspace of dimension "<<dimsofar;
         cout<<" which is a sum of oldclasses."<<endl;
         cout<<"List of ap:"<<aplist<<endl;
       }
      return;  //This branch of the recursion ends
    }
   if (depth==maxdepth)  // Can't happen as maxdepth>nplist always!
   { if (verbose)
     { cout << "\nFound a "<<dimsofar<<"D common eigenspace\n";
       cout << "List of ap: "<<aplist<<"\n";
       cout << "Abandoning, even though oldforms only make up ";
       cout << dimold << "D of this." << endl;
      }
     dimsplit += dimsofar;
     return;
   }

// The recursive part:

   int inewq=iq+1; Quad newq=primelist[iq];
   if (verbose) cout << "Using q = " << newq << endl;
   MAT t = RESTRICT(wmats[iq] , s);

   for (int a=1; (a>-2) ; a-=2)
   { int aq = a;
     SCALAR lambda = h1->h1denom()*aq*denom(s);
     SUBSP temp = EIGENSPACE(t,lambda);
     SUBSP newspace = COMBINE(s,temp);
     int newdim = dim(newspace);
     if (newdim>0)
     { int newdepth = depth+1;
       vector<long> newaplist(newdepth);
       for (int i=0; i<depth; i++) newaplist[i]=aplist(i);
       newaplist[depth] = aq;
       if (inewq==npdivs)
          tsplit(newspace,    0,newdepth,newaplist);
       else
          wsplit(newspace,inewq,newdepth,newaplist);
     }
   }
}

#endif // #ifndef USE_XSPLIT
 
//end of newforms.cc
