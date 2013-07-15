#include <iostream>
#include <sstream>
#include "oldforms.h"

inline int testbit(long a, long i) {return (a& (1<<i));}

// Implementation of eigdata constructor -- reads data from file
eigdata::eigdata(const Quad& m, int neigs, int verbose) :sublevel(m)
{
  if(verbose) cout << "Getting eigdata for " << m << endl;
  string eigfilename = eigfile(m);
  ifstream data(eigfilename.c_str());
  if (!data)
    {
      if(verbose)
        {cout << "No data file for m = " << m; 
         cout << "  so assuming no newforms at that level." << endl;
       }
      nforms=nforms2=nap=0;
      return;
    }
  int dump,i,neigsonfile;
  nforms = nforms2 = neigsonfile = 0;
  data >> nforms >> nforms2 >> neigsonfile;
  neigsonfile++;  //NB this is because the no. on file is the last index,
                  //which is 1 less than the number of eigs!
  if(verbose)
    cout<<nforms<<" newforms found at level "<<m<<", total new dimension = "<<(nforms+nforms2)<<endl;
  if(nforms>0){
    if(neigs<0)nap=neigs=neigsonfile;
    else nap= (neigsonfile<neigs?neigsonfile:neigs);
    int nwq=level::npdivs;
    int ntp=nap-nwq;
    eigs.resize(nforms);
    for(i=0; i<nforms; i++) eigs[i].resize(nap);
    int countp=0, countq=0;
    vector<Quad>::const_iterator pr;
    for (pr=quadprimes.begin(); 
	 (pr-quadprimes.begin())<=neigsonfile && ((countp<ntp) || (countq<nwq)); pr++)
      { 
	if(verbose) cout<<"p="<<(*pr)<<endl;
	if (div(*pr,level::modulus))
          {
            for(i=0; i<nforms; i++) data>>eigs[i][countq];
            countq++; 
          }
      else if (countp<ntp)
        {
          for(i=0; i<nforms; i++)  data>>eigs[i][nwq+countp];
          countp++;
        }
      else for(i=0; i<nforms; i++) data >> dump;
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
    if(verbose) cout << "Finished reading eigdata." << endl;
    if(verbose) cout << "eigs = " << endl;
    if(verbose) for(i=0; i<nforms; i++) cout<<i<<": "<<eigs[i]<<endl;
    data.close();
  }
}

// Implementation of oldform member functions

oldforms::oldforms(int verbose)
{  
   nap = level::nap;
   ntp = nap-level::npdivs;
   noldclasses=olddim1=olddim2=0;
   vector<Quad>::const_iterator d=(level::dlist).begin();
   while(d!=(level::dlist).end()) getoldclasses(*d++,verbose);
   for (int i=0; i<noldclasses; i++) olddim1+=oldclassdims[i];
   olddimall = olddim1+olddim2;
   if(verbose)
     {
       cout<<"Leaving oldform constructor with olddim1 = "<<olddim1;
       cout<<", olddim2 = "<<olddim2<<", olddimall="<<olddimall<<endl;
     }
}

void oldforms::getoldclasses(const Quad& d, int verbose) //really a subroutine of the
{                                                 //constructor
  Quad N = level::modulus;
  long normd = quadnorm(d);
  if ((normd>1) && (level::normod>normd))
    {
      if(verbose) cout << "Getting oldclasses for divisor " << d << endl;
      eigdata olddata(d,nap,verbose);
      int nforms=olddata.nforms;
      Quad m = N/d;
      int k=0, oldmult=1, xmult, mult, j, beta; 
      vector<long> betalist;
      vector<Quad>::const_iterator p=(level::plist).begin();
      while(p!=(level::plist).end()) 
        {
	  beta=val(*p++,m);
	  oldmult*=1+beta;
	  if(beta>0) k++;
	  betalist.push_back(beta);
	}
      olddim2+=oldmult*olddata.nforms2;
      if(verbose) cout << "Computing W multiplicities." << endl;
      vector<long> nextoldformap(nap);
      for(int iform=0; iform<nforms; iform++)
        { for (int c=0; c<(1<<k); c++)
            {  
               if(verbose) cout << "c = " << c << endl;
               mult=1; j=0;
	       vector<Quad>::const_iterator q;
               for (q=(level::plist).begin(); 
		    q!=(level::plist).end()&&(mult>0); q++)
                 {  int i = q-(level::plist).begin(); 
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
                   for(int i=level::npdivs; i<nap; i++)
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
  for (int i=0; i<noldclasses; i++)
    if (startswith(oldformap[i] , aplist, aplist.size())) 
      ans += oldclassdims[i];
  return ans;
}

void oldforms::display(void) const
{
  if (noldclasses>0)
  {
    cout << "\nOld classes\n~~~~~~~~~~~\n";
    cout << "Level   Dimension " << level::primelist << endl;
    for (int i=0; i<noldclasses; i++)
    { cout << oldlevels[i] << "       " << oldclassdims[i] << "       ";
      cout << oldformap[i] << endl;
    }
  }
 cout<<"Total number of (rational) oldclasses = "<<noldclasses<<endl;
 cout<<"Total dimension of (rational) oldclasses = "<<olddim1<<endl;
 cout<<"Total dimension of all oldclasses = "<<olddimall<<endl;
}


string ideal_code(const Quad& d) // string code for a (principal)  ideal
{
  stringstream s;
  long r=real(d), i=imag(d);
  if(r<0)    s << "m";
  s << abs(r);
  s << "i";
  if(i<0)    s << "m";
  s << abs(i);
  s << char(0);
  return s.str();
}

string eigfile(const Quad& d)    //returns filename for eigs at level d
{
  stringstream s;
  s << getenv("NF_DIR");
  if (s.str().empty()) {s.clear(); s<<"./newforms";}
  s << "/Qsqrt-" << Quad::d << "/e";
  s << ideal_code(d);
  return s.str();
}

