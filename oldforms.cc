#include <iostream>
#include <sstream>
#include "newforms.h"

inline int testbit(long a, long i) {return (a& (1<<i));}

string ideal_code(const Quad& N); // string code for a (principal)  ideal

string eigfile(const Quad& N)    //returns filename for eigs at level N
{
  stringstream s;
  s << getenv("NF_DIR");
  if (s.str().empty()) {s.clear(); s<<"./newforms";}
  s << "/2.0." << (Quad::disc) << ".1/";
  s << ideal_code(N);
  return s.str();
}

string eigfile(Qideal& N)    //returns filename for eigs at level N
{
  return eigfile(N.gen()); // temporary for principal ideals
}

// Implementation of eigdata constructor -- reads data from file
eigdata::eigdata(Qideal& iN, const Qideal& iM, int neigs, int verbose)
  : N(iN), M(iM)
{
  if(verbose) cout << "Getting eigdata for " << M << endl;
  string eigfilename = eigfile(M.gen());
  ifstream data(eigfilename.c_str());
  if (!data)
    {
      if(verbose)
        {
          cout << "No data file for M = " << M;
          cout << "  so creating newforms at that level..." << endl;
        }
      newforms olddata(M,verbose);
      olddata.createfromscratch();
      olddata.getap(1,20,0);
      olddata.output_to_file(eigfilename);
      if(verbose)
        {
          cout << "  finished creating newforms at level " << M << endl;
          olddata.display();
        }
      data.open(eigfilename.c_str());
    }
  int neigsonfile;
  nforms = nforms2 = neigsonfile = 0;
  data >> nforms >> nforms2 >> neigsonfile;
  if(verbose)
    {
      cout<<"neigs="<<neigs<<", neigsonfile = "<<neigsonfile<<endl;
      cout<<nforms<<" newforms found at level "<<M<<", total new dimension = "<<(nforms+nforms2)<<endl;
    }
  if(nforms>0){
    if(neigs<0)nap=neigs=neigsonfile;
    else nap= (neigsonfile<neigs?neigsonfile:neigs);
    int nwq=N.factorization().size();
    int ntp=nap-nwq;
    vector<Quadprime> qlist = pdivs(M);
    int nq = qlist.size();
    eigs.resize(nforms);
    aps.resize(nforms);
    aqs.resize(nforms);
    intdata.resize(nforms);
    Quaddata.resize(nforms);
    int i;
    for(i=0; i<nforms; i++)
      {
	eigs[i].resize(nap);
	aqs[i].resize(nq);
	aps[i].resize(nap);
        intdata[i].resize(6);
        Quaddata[i].resize(5);
      }

    // Read the auxiliary data:
    for (i=0; i<nforms; i++) data>>intdata[i][0];  // sfe
    for (i=0; i<nforms; i++) data>>intdata[i][1];  // pdot
    for (i=0; i<nforms; i++) data>>intdata[i][2];  // dp0
    for (i=0; i<nforms; i++) data>>intdata[i][3];  // cuspidalfactor
    for (i=0; i<nforms; i++) data>>Quaddata[i][0]; // lambda
    for (i=0; i<nforms; i++) data>>intdata[i][4];  // lambdadot
    for (i=0; i<nforms; i++) data>>Quaddata[i][1]; // a
    for (i=0; i<nforms; i++) data>>Quaddata[i][2]; // b
    for (i=0; i<nforms; i++) data>>Quaddata[i][3]; // c
    for (i=0; i<nforms; i++) data>>Quaddata[i][4]; // d
    for (i=0; i<nforms; i++) data>>intdata[i][5];  // matdot

    //  Read the W-eigenvalues at level M into aqs:
    vector<vector<long> >::iterator f; long eig;
    for(i=0; i<nq; i++)
      for(f=aqs.begin(); f!=aqs.end(); ++f)
	{
	  data>>eig;
	  (*f)[i]=eig;
	}

    // Next read the coefficients at level M into aps:
    for(i=0; i<neigs; i++)
      for(f=aps.begin(); f!=aps.end(); ++f)
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
	cout << "intdata = " << endl;
	for(i=0; i<nforms; i++) cout<<i<<": "<<intdata[i]<<endl;
	cout << "Quaddata = " << endl;
	for(i=0; i<nforms; i++) cout<<i<<": "<<Quaddata[i]<<endl;
      }

    // Now construct the eigenvalue sequence, first Wq eigenvalues for
    // bad primes then Tp-eigenvalues for good primes

    int countp=0, countq=0;
    for (vector<Quadprime>::const_iterator Pi=Quadprimes::list.begin();
	 ((countp<ntp) || (countq<nwq));
	 ++Pi)
      {
	int pindex = Pi - Quadprimes::list.begin();
        Quadprime P = *Pi;
	if (P.divides(N))
          {
	    if(verbose)
	      cout<<"P="<<P<<" = bad prime # "<<countq<<" [";
	    // if P also divides M we can pick up the W-eigenvalue,
	    // otherwise the value is not needed.  Note that N may
	    // have more prime factors than M, so the index in the aqs
	    // may be different from countq.
	    if (P.divides(M)) // P divides M (and N)
	      {
                // find the index j of P in the list of prime divisors of M:
		int j = find(qlist.begin(),qlist.end(),P)-qlist.begin();
		for(i=0; i<nforms; i++)
		  {
		    eigs[i][countq] = aqs[i][j];
		    if(verbose) cout<<" "<<eigs[i][countq];
		  }
	      }
	    else // P divides N but not M
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
	      cout<<"P="<<P<<" = good prime # "<<countp<<" [";
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
    if(verbose)
      {
	cout << "eigs = " << endl;
	for(i=0; i<nforms; i++) cout<<i<<": "<<eigs[i]<<endl;
      }
  }
}

// Implementation of oldform member functions

//This must include newforms in the minus space for tests to pass!
static long min_newform_level_norm[20] = {0,65, // d=1
                                          32,   // d=2
                                          49,   // d=3
                                          0,0,25,           // d=7
                                          0,0,0,0,9,        // d=11
                                          0,0,0,0,0,0,0,4}; // d=19

oldforms::oldforms(Qideal& iN, const vector<Quadprime>& pr, int verbose)
  : N(iN), plist(pr)
{
   nap = plist.size();
   noldclasses = olddim1 = olddim2 = 0; // will be incremented in getoldclasses()
   vector<Qideal> DD = alldivs(N);
   for(vector<Qideal>::iterator Di = DD.begin(); Di!=DD.end(); ++Di)
     {
       getoldclasses(*Di,verbose); // will skip D==N
     }
   for (int i=0; i<noldclasses; i++)
     olddim1 += oldclassdims[i];
   olddimall = olddim1 + olddim2;
   if(verbose)
     {
       cout<<"Leaving oldform constructor with olddim1 = "<<olddim1;
       cout<<", olddim2 = "<<olddim2<<", olddimall="<<olddimall<<endl;
     }
}

//really a subroutine of the constructor
void oldforms::getoldclasses(Qideal& D, int verbose)
{
  if (D==N)
    return;
  long min_norm = (Quad::d <=19? min_newform_level_norm[Quad::d]: 1);
  if (D.norm() < min_norm)
    {
      if(verbose)
        cout<<"Skipping oldforms from sublevel "<<D<<" of norm "<<D.norm()
            <<" which is less than "<< min_norm <<endl;
      return;
    }
  if(verbose)
    cout << "Getting oldclasses for divisor " << D << endl;
  eigdata olddata(N,D,nap,verbose);
  int nforms=olddata.nforms;
  Qideal M = N/D;
  int k=0, oldmultiplicity=1, xmultiplicity, multiplicity, j, beta;
  vector<long> betalist;
  for (vector<Quadprime>::const_iterator Pi = plist.begin(); Pi!=plist.end(); ++Pi)
    {
      beta = val(*Pi, M);
      oldmultiplicity *= 1+beta;
      if(beta>0) k++;
      betalist.push_back(beta);
    }
  if(verbose) cout<<"betas="<<betalist<<", each oldspace dimension is "<<oldmultiplicity<<endl;
  olddim2+=oldmultiplicity*olddata.nforms2;
  if(verbose) cout << "Computing W multiplicities." << endl;
  vector<long> nextoldformap(nap);
  for(int iform=0; iform<nforms; iform++)
    { for (int c=0; c<(1<<k); c++)
        {
          if(verbose) cout << "c = " << c << endl;
          multiplicity=1; j=0;
          ;
          for (vector<Quadprime>::const_iterator Qi = plist.begin(); Qi != plist.end()&&(multiplicity>0); ++Qi)
            {
              int i = Qi - plist.begin();
              beta=betalist[i];
              if (beta>0)
                { int bit = testbit(c,j); j++;
                  nextoldformap[i] = bit?1:-1;
                  if (odd(beta)) xmultiplicity =  (beta+1)/2;
                  else if((*Qi).divides(D) && (olddata.eigs[iform][i]==-1)) xmultiplicity=beta/2;
                  else xmultiplicity=(beta/2)+1;
                  if (!bit) xmultiplicity=1+beta-xmultiplicity;
                  multiplicity*=xmultiplicity;
                }
              else nextoldformap[i] = olddata.eigs[iform][i];
            }
          if(verbose) cout << "Multiplicity = " << multiplicity << endl;
          if (multiplicity>0)
            {
              for(int i=plist.size(); i<nap; i++)
                nextoldformap[i] = olddata.eigs[iform][i];
              oldformap.push_back(nextoldformap);
              oldclassdims.push_back(multiplicity);
              oldlevels.push_back(D);
              noldclasses++;
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
    cout << "\nOld classes for level "<<N<<"\n~~~~~~~~~~~\n";
    cout << "Level   Dimension " << plist << endl;
    for (int i=0; i<noldclasses; i++)
    { cout << oldlevels[i] << "       " << oldclassdims[i] << "       ";
      cout << oldformap[i] << endl;
    }
  }
 cout<<"Total number of (rational) oldclasses = "<<noldclasses<<endl;
 cout<<"Total dimension of (rational) oldclasses = "<<olddim1<<endl;
 cout<<"Total dimension of all oldclasses = "<<olddimall<<endl;
}

