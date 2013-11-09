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

newform::newform(const newforms* nfs, const vec& v, const vector<long>& ap)
 :basis(v)
{
  //cout<<"Constructing newform with aplist = "<<ap<<endl;
  copy(ap.begin(), ap.begin()+(nfs->npdivs), back_inserter(aqlist));
  dp0    =  1 + quadnorm(nfs->p0) - ap[nfs->npdivs];
  pdot   =  (nfs->mvp)*basis;
  if(nfs->hmod) pdot=mod(pdot,nfs->hmod);
  //cout<<"dp0 = "<<dp0<<", pdot="<<pdot<<endl;
  loverp = rational(abs(pdot),dp0*(Quad::nunits));  //No division by h1denom()
  sfe = std::accumulate(aqlist.begin(),aqlist.end(),-1,std::multiplies<long>());
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

  if(n1ds==0) return; // else no work to do

// Find the twisting primes for each newform (more efficient to do
// this here instead of within the newform constructors, as one lambda
// might work for more than one newform). NB If the level is square
// SFE=-1 then no such lambda will exist.

  get_lambda();

  makeh1plus(); // In case it wasn't done in the newforms constructor
  allproj();    // Compute homspace::projcoord, so projcycle can be used
  find_jlist();
}

void newforms::get_lambda()
{
//#define DEBUG_LAMBDA
  int* gotlambda = new int[n1ds];
  int i, nfound=0;
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
    else 
      {
        nflist[i].lambda=Quad(0); // indicates no lambda exists
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
                      int dot = abs(mvtw*nfj.basis);
                      if(hmod) dot=mod(dot,hmod);
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
  if(upperbound<0) // check for error condition
    {
      cout<<"Error:  total old dimension = "<<olddimall<<" as computed is greater than total cuspidal dimension "<<(h1->h1cuspdim())<<" -- aborting"<<endl;
      exit(1);
    }
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
      cout << "Error in splitting eigenspaces: apparently found more ";
      cout << "1D newforms ("<< n1ds+1 <<") than the total new-dimension ("
           <<upperbound<<").\n";
      exit(1);
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
 cout << "basis = " << basis 
      << ";\taqlist = " << aqlist 
      << ";\taplist = " << aplist << endl;
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

scalar dotmodp(const vec& v1, const vec& v2, scalar pr)
{
  scalar ans=0;
  for(long i=1; i<=dim(v1); i++) ans=xmod(ans+xmodmul(v1[i],v2[i],pr),pr);
  return mod(ans,pr);
}

void newforms::allproj() //Replaces "coord" member of homspace with projections
                      //onto eigenspaces, to save time
{
  int ncoord = nrows(h1->coord); long pcij;
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

void newforms::getoneap(const Quad& p, int verbose)
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
      nflist[i].aplist.push_back(cp);
    }
  if(verbose)
    cout << endl;
}

void newforms::getap(int first, int last, int verbose)
{
  nap=0;
  if(n1ds>0)
    {
      vector<Quad>::const_iterator pr;
      for(pr=quadprimes.begin(); pr!=quadprimes.end(); pr++)
	{
	  long index = pr-quadprimes.begin()+1;
	  if ( (index>=first) && (index<=last) )
	    getoneap(*pr,verbose);
	}
      nap = nflist[0].aplist.size();
    }
}

void newforms::output_to_file(string eigfile) const
{
  int echo=0;
  ofstream out;
  out.open(eigfile.c_str());
  out<<n1ds<<" "<<n2ds<<endl;
  if(echo) cout<<n1ds<<" "<<n2ds<<endl;
  out<<nap<<endl;
  if(echo) cout<<nap<<endl;
  vector<newform>::const_iterator f;
  vector<long>::const_iterator ap;
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
  out<<endl;
  if(echo) cout<<endl;
  for(int i=0; i<nap; i++)
    {
      for(f=nflist.begin(); f!=nflist.end(); f++)
	{
	  out<<setw(5)<<(f->aplist)[i]<<" ";
	  if(echo) cout<<setw(5)<<(f->aplist)[i]<<" ";
	}
      out<<endl;
      if(echo) cout<<endl;
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

  map<long,vec> images; // [j,v] stores image of j'th M-symbol in v
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
      if(ind) update(pcd,imagej,ind,h1->hmod);
// Matrix [p,0;0,1]
      ind = h1->coordindex[h1->index2(p*u,v)];
#ifdef DEBUG_APVEC
      cout<<"u1="<<p*u<<", u2="<<v<<", ind="<<ind<<endl;
#endif
      if(ind) update(pcd,imagej,ind,h1->hmod);
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
          if(ind) update(pcd,imagej,ind,h1->hmod);
          while(b!=0)
            {
              q=a/b; c=a-b*q; u3=q*u2-u1;
              a=-b; b=c; u1=u2; u2=u3;
              ind = h1->coordindex[h1->index2(u1,u2)];
#ifdef DEBUG_APVEC
              cout<<"a="<<a<<", b="<<b<<", u1="<<u1<<", u2="<<u2<<", ind="<<ind<<endl;
#endif
              if(ind) update(pcd,imagej,ind,h1->hmod);
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
      apv[i]=ap=images[nflist[i].j0][i+1]/nflist[i].fac;
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

void newforms::addap(long last) // adds ap for primes up to the last'th prime
{
}

//end of newforms.cc
