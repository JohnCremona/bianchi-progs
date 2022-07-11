#include <iostream>
#include <sstream>
#include "newforms.h"

// Implementation of oldform member functions

oldforms::oldforms(const newforms* nfs)
  : n2r(Quad::class_group_2_rank), nf(nfs)
{
  nchi = 1<<n2r;
  old1dims.resize(nchi); // will sum to olddim1
  old2dims.resize(nchi); // will sum to olddim2
  olddims.resize(nchi);  // will sum to olddimall
  noldclasses = 0;       // will be incremented in getoldclasses()
  olddim1 = olddim2 = olddimall = 0;

  N = nf->N;
  vector<Qideal> DD = alldivs(N);
  for(auto Di = DD.begin(); Di!=DD.end(); ++Di)
    if (*Di!=N)
      getoldclasses(*Di);

  assert(olddim1   == std::accumulate(old1dims.begin(), old1dims.end(), 0, std::plus<int>()));
  assert(olddim2   == std::accumulate(old2dims.begin(), old2dims.end(), 0, std::plus<int>()));
  assert(olddimall == std::accumulate(olddims.begin(), olddims.end(), 0, std::plus<int>()));

  if (nf->verbose >1)
  {
    cout<<"noldclasses = "<<noldclasses<<endl;
    cout<<"olddim1 = "  <<olddim1;
    if (nchi>1) cout << " (by character: "<<old1dims<<")";
    cout<<endl;
    cout<<"olddim2 = "  <<olddim2;
    if (nchi>1) cout << " (by character: "<<old2dims<<")";
    cout<<endl;
    cout<<"olddimall = "<<olddimall;
    if (nchi>1) cout << " (by character: "<<olddims<<")";
    cout<<endl;
  }
  assert (olddimall == olddim1 + olddim2);

  if(nf->verbose)
    {
      cout<<"Leaving oldform constructor with olddim1 = "<<olddim1
          <<", olddim2 = "<<olddim2
          <<", olddimall="<<olddimall<<endl;
    }
}

//really a subroutine of the constructor
void oldforms::getoldclasses(Qideal& D)
{
  if(nf->verbose)
    cout << "\nGetting oldclasses at level "<<ideal_label(N)<<" from divisor " << ideal_label(D) << endl;

  newforms olddata(D, nf->verbose, nf->characteristic);
  olddata.read_from_file_or_find();
  int old1ds=olddata.n1ds, old2ds=olddata.n2ds;
  noldclasses += old1ds;

  if (old1ds==0 && old2ds==0)
    return;

  // Compute the oldform multiplicities.

  // Note that in even class number, the multiplicity may be different
  // for different newforms coming from the same level, if some of
  // them are self-twist.

  Qideal M = N/D;
  vector<Qideal> divisors = alldivs(M);
  int oldmult = divisors.size(); // the usual multiplicity

  vector<int>oldmults(old1ds, oldmult);   // list of multiplicities of each newform
  for(int i=0; i<old1ds; i++)
    {
      QUINT CMD = olddata.nflist[i].CMD;
      if (!is_zero(CMD))
        oldmults[i] = old_multiplicity(CMD, divisors);
    }

  if(nf->verbose)
    cout<<" oldspace dimensions for "<<old1ds<<" rational forms are "<<oldmults<<endl;

  for(int iform=0; iform<old1ds; iform++)
    {
      oldformap.push_back(olddata.nflist[iform].oldform_eigs(N));
      oldclassdims.push_back(oldmults[iform]);
      oldlevels.push_back(D);
    }

  vector<int> old1dimsD = old_multiplicities(olddata.new1dims, divisors);
  vector<int> old2dimsD = old_multiplicities(olddata.new2dims, divisors);
  int this_olddim1=0, this_olddim2=0, this_olddimall=0;
  for (int i=0; i<nchi; i++)
    {
      old1dims[i] += old1dimsD[i];
      old2dims[i] += old2dimsD[i];
      this_olddim1 += old1dimsD[i];
      this_olddim2 += old2dimsD[i];
      olddims[i] += (old1dimsD[i] + old2dimsD[i]);
      this_olddimall += (old1dimsD[i] + old2dimsD[i]);
    }
  olddim1 += this_olddim1;
  olddim2 += this_olddim2;
  olddimall += this_olddimall;
  if(nf->verbose)
    {
      cout<<" total oldspace dimension from divisor "<<D<<" is "
          <<this_olddim1<<"+"<<this_olddim2<<"="<<this_olddimall<<endl;
      cout<<" cumulative total oldspace dimension from divisors so far is "
          <<olddim1<<"+"<<olddim2<<"="<<olddimall<<endl;
    }
}

long oldforms::dimoldpart(vector<long> aplist)
{
  if (noldclasses==0) return 0;          // no oldforms
  if (nf->characteristic!=0) return 0;   // until we work out how to compute this
  if (aplist.size()==0) return 0;        // all lists "start with" a null list!

  int debug=0;
  if (debug) cout<<"In dimoldpart with aplist="<<aplist<<endl;
  int ans = 0;
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
    cout << "Level   Dimension ";
    int r = Quad::class_group_2_rank;
    if (r>0)
      {
        cout<<"[";
        while(r--) cout<<" nu";
        cout<<" ] ";
      }
    cout << nf->goodprimes << endl;
    for (int i=0; i<noldclasses; i++)
    { cout << oldlevels[i] << "       " << oldclassdims[i] << "       ";
      cout << oldformap[i] << endl;
    }
  }
 cout<<"Total number of (rational) oldclasses = "<<noldclasses<<endl;
 cout<<"Total dimension of (rational) oldclasses = "<<olddim1<<endl;
 cout<<"Total dimension of all oldclasses = "<<olddimall<<endl;
}


// Usually, the oldform multiplicity (i.e. the dimension of the
// oldspace) at level N from a newform at level D is the number of
// divisors of N/D, but if the class number is even and the form is
// self-twist by an unramified quadratic character chi, we only count
// the divisors in the kernel of chi.


// Return the oldspace dimension at level N of a (rational) newform at
// level D which is self-twist by discriminant d
int old_multiplicity(Qideal D, QUINT d, Qideal N)
{
  Qideal M = N/D;
  vector<Qideal> divisors = alldivs(M);
  return old_multiplicity(d, divisors);
}

int old_multiplicity(QUINT d, vector<Qideal>& divisors)
{
  int mult = 0;
  for_each(divisors.begin(), divisors.end(),
           [d, &mult](Qideal D)
           {mult += int(D.genus_character(d)==+1);}
           );
  return mult;
}

// Given the new dimensions at level D and a multiple N of D, return
// the oldspace dimensions at level N
vector<int> old_multiplicities(Qideal D, vector<int> newdimsD, Qideal N)
{
  Qideal M = N/D;
  vector<Qideal> divisors = alldivs(M);
  return old_multiplicities(newdimsD, divisors);
}

// The same with the list of divisors of N/D given
vector<int> old_multiplicities(vector<int> newdimsD, vector<Qideal>& divisors)
{
  vector<int> ans(newdimsD.size());
  std::transform(newdimsD.begin(), newdimsD.end(), Quad::all_disc_factors.begin(),
                 ans.begin(),
                 [&divisors](int d, QUINT D) {return d*old_multiplicity(D, divisors);}
                 );
  return ans;
}
