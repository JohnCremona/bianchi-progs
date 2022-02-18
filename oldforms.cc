#include <iostream>
#include <sstream>
#include "newforms.h"

string ideal_code(const Quad& N); // string code for a (principal)  ideal

string eigfile(const Quad& N, long p)    //returns filename for eigs at level N, characteristic p
{
  stringstream s;
  s << getenv("NF_DIR");
  if (s.str().empty()) {s.clear(); s<<"./newforms";}
  s << "/" << field_label() << "/";
  s << ideal_code(N);
  if (p)
    s << "_mod_"<<p;
  return s.str();
}

string eigfile(Qideal& N, long p)    //returns filename for eigs at level N, characteristic p
{
  stringstream s;
  s << getenv("NF_DIR");
  if (s.str().empty()) {s.clear(); s<<"./newforms";}
  s << "/" << field_label() << "/";
  if (Quad::class_number==1) // for backwards compatibility of data file names
    s << ideal_code(N.gen());
  else
    s << ideal_label(N);
  if (p)
    s << "_mod_"<<p;
  return s.str();
}

// Implementation of oldform member functions

//This must include newforms in the minus space for tests to pass!
static long min_newform_level_norm[20] = {0,65, // d=1
                                          32,   // d=2
                                          49,   // d=3
                                          0,0,25,           // d=7
                                          0,0,0,0,9,        // d=11
                                          0,0,0,0,0,0,0,4}; // d=19

oldforms::oldforms(const newforms* nfs)
  : nf(nfs)
{
  plist = nf->goodprimes;
  plist.insert(plist.end(), nf->goodprimes.begin(), nf->goodprimes.end());
  nap = plist.size();
  noldclasses = olddimall = olddim1 = olddim2 = 0; // will be incremented in getoldclasses()
  if (nf->characteristic>0)
    return;
  N = nf->N;
  vector<Qideal> DD = alldivs(N);
  if (nf->verbose>1)
    {
      cout<<"List of all divisors of N: ";
      for(vector<Qideal>::iterator Di = DD.begin(); Di!=DD.end(); ++Di)
        {
          if (Di!=DD.begin())
            cout << ", ";
          cout << ideal_label(*Di) << " = " << (*Di);
        }
      cout << endl;
    }
  for(vector<Qideal>::iterator Di = DD.begin(); Di!=DD.end(); ++Di)
    {
      getoldclasses(*Di); // will skip D==N
    }
  for (int i=0; i<noldclasses; i++)
    olddim1 += oldclassdims[i];
  olddimall = olddim1 + olddim2;
  if(nf->verbose)
    {
      cout<<"Leaving oldform constructor with olddim1 = "<<olddim1;
      cout<<", olddim2 = "<<olddim2<<", olddimall="<<olddimall<<endl;
    }
}

//really a subroutine of the constructor
void oldforms::getoldclasses(Qideal& D)
{
  if (D==N)
    return;
  long min_norm = 1;
  if (Quad::d <=19 && (nf->characteristic==0))
    min_norm = min_newform_level_norm[Quad::d];
  if (D.norm() < min_norm)
    {
      if(nf->verbose)
        cout<<"Skipping oldforms from sublevel "<<ideal_label(D)<<" whose norm "<<D.norm()
            <<" is less than "<< min_norm <<endl;
      return;
    }
  if(nf->verbose)
    cout << "\nGetting oldclasses for divisor " << ideal_label(D) << endl;

  newforms olddata(D, nf->verbose, nf->characteristic);
  olddata.read_from_file_or_find();
  int nforms=olddata.n1ds;

  // Compute the oldform multiplicity as the number of divisors of N/D
  int oldmult;
  if (nf->characteristic>0)
    oldmult=0;
  else
    {
      vector<Quadprime>::const_iterator Pi;
      for (Pi = nf->allbadprimes.begin(), oldmult=1; Pi!=nf->allbadprimes.end(); ++Pi)
        {
          int fac = (1 + val(*Pi, N) - val(*Pi, D));
          //cout<<" prime "<<(*Pi)<<" contributes a factor of "<<fac<<" to the oldspace multiplicity"<<endl;
          oldmult *= fac;
        }
    }

  // As we are not using W's, all oldform multiplicities are equal to
  // oldmult, we do not need to subdivide this into 2^k pieces.

  if(nf->verbose) cout<<" each oldspace dimension is "<<oldmult<<endl;
  olddim2+=oldmult*olddata.n2ds;

  for(int iform=0; iform<nforms; iform++)
    {
      oldformap.push_back(olddata.nflist[iform].oldform_eigs(N));
      oldclassdims.push_back(oldmult);
      oldlevels.push_back(D);
      noldclasses++;
    }
}
      // for (int c=0; c<(1<<k); c++)
      //   {
      //     if(nf->verbose) cout << "c = " << c << endl;
      //     multiplicity=1; j=0;

      //     for (vector<Quadprime>::const_iterator Qi = plist.begin(); Qi != plist.end()&&(multiplicity>0); ++Qi)
      //       {
      //         int i = Qi - plist.begin();
      //         beta=betalist[i];
      //         if (beta>0)
      //           { int bit = testbit(c,j); j++;
      //             nextoldformap[i] = bit?1:-1;
      //             if (odd(beta)) xmultiplicity =  (beta+1)/2;
      //             else if((*Qi).divides(D) && (oldeigs[iform][i]==-1)) xmultiplicity=beta/2;
      //             else xmultiplicity=(beta/2)+1;
      //             if (!bit) xmultiplicity=1+beta-xmultiplicity;
      //             multiplicity*=xmultiplicity;
      //           }
      //         else nextoldformap[i] = oldeigs[iform][i];
      //       }
      //     if(nf->verbose) cout << "Multiplicity = " << multiplicity << endl;
      //     if (multiplicity>0)
      //       {
      //         for(int i=plist.size(); i<nap; i++)
      //           nextoldformap[i] = oldeigs[iform][i];
      //         oldformap.push_back(nextoldformap);
      //         oldclassdims.push_back(multiplicity);
      //         oldlevels.push_back(D);
      //         noldclasses++;
      //       }
      //   }
      //}

long oldforms::dimoldpart(vector<long> aplist)
{ int ans = 0;
  int debug=0;
  if (noldclasses==0) return 0;          // no oldforms
  if (nf->characteristic!=0) return 0;   // until we work out how to compute this
  if (aplist.size()==0) return 0;        // all lists "start with" a null list!
  if (debug) cout<<"In dimoldpart with aplist="<<aplist<<endl;
  for (int i=0; i<noldclasses; i++)
    {
      if (debug) cout<<"           oldformap["<<i<<"]="<<oldformap[i]<<endl;
      if (startswith(oldformap[i] , aplist, aplist.size()))
        ans += oldclassdims[i];
    }
  if (debug) cout<<"dimoldpart: total="<<ans<<endl;
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

