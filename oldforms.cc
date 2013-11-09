#include <iostream>
#include <sstream>
#include "oldforms.h"
#include "newforms.h"

inline int testbit(long a, long i) {return (a& (1<<i));}

// Implementation of eigdata constructor -- reads data from file
eigdata::eigdata(const level *iN, const Quad& m, int neigs, int verbose) :sublevel(m)
{
  N = iN;
  if(verbose) cout << "Getting eigdata for " << m << endl;
  string eigfilename = eigfile(m);
  ifstream data(eigfilename.c_str());
  if (!data)
    {
      if(verbose)
        {
          cout << "No data file for m = " << m;
          cout << "  so creating newforms at that level..." << endl;
        }
      newforms olddata(m,0,verbose);
      olddata.getap(1,iN->nap,0);
      olddata.output_to_file(eigfilename);
      if(verbose)
        {
          cout << "  finished creating newforms at level " << m << endl;
          olddata.display();
        }
      data.open(eigfilename.c_str());
    }
  int i,j,neigsonfile;
  nforms = nforms2 = neigsonfile = 0;
  data >> nforms >> nforms2 >> neigsonfile;
  if(verbose)
    {
      cout<<"neigs="<<neigs<<", neigsonfile = "<<neigsonfile<<endl;
      cout<<nforms<<" newforms found at level "<<m<<", total new dimension = "<<(nforms+nforms2)<<endl;
    }
  if(nforms>0){
    if(neigs<0)nap=neigs=neigsonfile;
    else nap= (neigsonfile<neigs?neigsonfile:neigs);
    int nwq=N->npdivs;
    int ntp=nap-nwq;
    vector<Quad> qlist = pdivs(m);
    int nq = qlist.size();
    vector<vector<long> > aps, aqs;
    eigs.resize(nforms);
    aps.resize(nforms);
    aqs.resize(nforms);
    for(i=0; i<nforms; i++)
      {
	eigs[i].resize(nap);
	aqs[i].resize(nq);
	aps[i].resize(nap);
      }

    // First read the W-eigenvalues at level M into aqs:
    vector<vector<long> >::iterator f; long eig;
    for(i=0; i<nq; i++)
      for(f=aqs.begin(); f!=aqs.end(); f++)
	{
	  data>>eig;
	  (*f)[i]=eig;
	}

    // Next read the coefficients at level M into aps:
    for(i=0; i<neigs; i++)
      for(f=aps.begin(); f!=aps.end(); f++)
	{
	  data>>eig;
	  (*f)[i]=eig;
	}

    data.close();

    if(verbose)
      {
	cout << "Finished reading eigdata." << endl;
	cout << "aqs = " << endl;
	for(i=0; i<nforms; i++) cout<<i<<": "<<aqs[i]<<endl;
	cout << "aps = " << endl;
	for(i=0; i<nforms; i++) cout<<i<<": "<<aps[i]<<endl;
      }

    // Now construct the eigenvalue sequence, first Wq eigenvalues for
    // bad primes then Tp-eigenvalues for good primes

    int countp=0, countq=0, pindex;
    vector<Quad>::const_iterator pr;
    for (pr=quadprimes.begin();
	 (pr-quadprimes.begin())<=neigsonfile && ((countp<ntp) || (countq<nwq));
	 pr++)
      {
	pindex = pr-quadprimes.begin();
	if (div(*pr,N->modulus))
          {
	    if(verbose)
	      cout<<"p="<<(*pr)<<" = bad prime # "<<countq<<" [";
	    // if p also divides m we can pick up the W-eigenvalue,
	    // otherwise the value is not needed.  Note that N may
	    // have more prime factors than m, so the index in the aqs
	    // may be different from countq.
	    if (div(*pr,m)) // pr divides m (and N)
	      {
		j = find(qlist.begin(),qlist.end(),*pr)-qlist.begin();
		for(i=0; i<nforms; i++)
		  {
		    eigs[i][countq] = aqs[i][j];
		    if(verbose) cout<<" "<<eigs[i][countq];
		  }
	      }
	    else // pr divides N but not m
	      {
		for(i=0; i<nforms; i++)
		  {
		    eigs[i][countq] = 999; // dummy value, will be overwritten
		    if(verbose) cout<<" *";//<<eigs[i][countq];
		  }
	      }
	    countq++;
	    if(verbose) cout<<" ]"<<endl;
          }
	else if (countp<ntp)
	  {
	    if(verbose)
	      cout<<"p="<<(*pr)<<" = good prime # "<<countp<<" [";
	    for(i=0; i<nforms; i++)
	      {
		eigs[i][nwq+countp] = aps[i][pindex];
		if(verbose) cout<<" "<<eigs[i][nwq+countp];
	      }
	    countp++;
	    if(verbose) cout<<" ]"<<endl;
	  }
      }
    if(countp<ntp)
      {
	cout<<"Error: not enough T_p eigs in file "
	    <<eigfilename<<endl;
      }
    if(countq<nwq)
      {
	cout<<"Error: not enough W_q eigs in file "
	    <<eigfilename<<endl;
      }
    if(verbose)
      {
	cout << "eigs = " << endl;
	for(i=0; i<nforms; i++) cout<<i<<": "<<eigs[i]<<endl;
      }
  }
}

// Implementation of oldform member functions

//This must include newforms in the minus space for tests to pass!
static long min_newform_level_norm[12] = {0,65,
                                          32,
                                          49,
                                          0,0,25,
                                          0,0,0,0,9};

oldforms::oldforms(const level* iN, int verbose)
{
   N = iN;
   nap = N->nap;
   ntp = nap-N->npdivs;
   noldclasses=olddim1=olddim2=0;
   vector<Quad>::const_iterator d=(N->dlist).begin();
   long min_norm = min_newform_level_norm[Quad::d];
   while(d!=(N->dlist).end())
     {
       if (quadnorm(*d)<min_norm)
         {
           if(verbose)
             cout<<"Skipping oldforms from sublevel "<<(*d)<<" of norm "<<quadnorm(*d)<<" which is less than "<< min_norm <<endl;
         }
       else
         {
           getoldclasses(*d,verbose);
         }
       d++;
     }
   for (int i=0; i<noldclasses; i++) olddim1+=oldclassdims[i];
   olddimall = olddim1+olddim2;
   if(verbose)
     {
       cout<<"Leaving oldform constructor with olddim1 = "<<olddim1;
       cout<<", olddim2 = "<<olddim2<<", olddimall="<<olddimall<<endl;
     }
}

//really a subroutine of the constructor
void oldforms::getoldclasses(const Quad& d, int verbose)
{
  long normd = quadnorm(d);
  if ((normd>1) && (N->normod>normd))
    {
      if(verbose) cout << "Getting oldclasses for divisor " << d << endl;
      eigdata olddata(N,d,nap,verbose);
      int nforms=olddata.nforms;
      Quad m = N->modulus/d;
      int k=0, oldmult=1, xmult, mult, j, beta;
      vector<long> betalist;
      vector<Quad>::const_iterator p=(N->plist).begin();
      while(p!=(N->plist).end())
        {
	  beta=val(*p++,m);
	  oldmult*=1+beta;
	  if(beta>0) k++;
	  betalist.push_back(beta);
	}
      if(verbose) cout<<"betas="<<betalist<<", each oldspace dimension is "<<oldmult<<endl;
      olddim2+=oldmult*olddata.nforms2;
      if(verbose) cout << "Computing W multiplicities." << endl;
      vector<long> nextoldformap(nap);
      for(int iform=0; iform<nforms; iform++)
        { for (int c=0; c<(1<<k); c++)
            {
               if(verbose) cout << "c = " << c << endl;
               mult=1; j=0;
	       vector<Quad>::const_iterator q;
               for (q=(N->plist).begin();
		    q!=(N->plist).end()&&(mult>0); q++)
                 {  int i = q-(N->plist).begin();
		    beta=betalist[i];
                    if (beta>0)
                      { int bit = testbit(c,j); j++;
			nextoldformap[i] = bit?1:-1;
                        if (odd(beta)) xmult =  (beta+1)/2;
                        else if(div(*q,d) && (olddata.eigs[iform][i]==-1)) xmult=beta/2;
                        else xmult=(beta/2)+1;
                        if (!bit) xmult=1+beta-xmult;
                        mult*=xmult;
                      }
		    else nextoldformap[i] = olddata.eigs[iform][i];
                  }
               if(verbose) cout << "Multiplicity = " << mult << endl;
               if (mult>0)
                 {
                   for(int i=N->npdivs; i<nap; i++)
		     nextoldformap[i] = olddata.eigs[iform][i];
		   oldformap.push_back(nextoldformap);
		   oldclassdims.push_back(mult);
		   oldlevels.push_back(d);
                   noldclasses++;
                 }
             }
        }
    }
}

long oldforms::dimoldpart(vector<long> aplist)
{ int ans = 0;
  if (aplist.size()==0) return 0;   // all lists "start with" a null list!
  //  cout<<"dimoldpart: aplist="<<aplist<<endl;
  for (int i=0; i<noldclasses; i++)
    {
      //      cout<<"           oldformap["<<i<<"]="<<oldformap[i]<<endl;
      if (startswith(oldformap[i] , aplist, aplist.size()))
        ans += oldclassdims[i];
    }
  //  cout<<"dimoldpart: ans="<<ans<<endl;
  return ans;
}

void oldforms::display(void) const
{
  if (noldclasses>0)
  {
    cout << "\nOld classes\n~~~~~~~~~~~\n";
    cout << "Level   Dimension " << N->primelist << endl;
    for (int i=0; i<noldclasses; i++)
    { cout << oldlevels[i] << "       " << oldclassdims[i] << "       ";
      cout << oldformap[i] << endl;
    }
  }
 cout<<"Total number of (rational) oldclasses = "<<noldclasses<<endl;
 cout<<"Total dimension of (rational) oldclasses = "<<olddim1<<endl;
 cout<<"Total dimension of all oldclasses = "<<olddimall<<endl;
}

