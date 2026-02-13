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
  for( auto& D : DD)
    if (D!=N)
      getoldclasses(D);

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
    cout << "\nGetting oldclasses at level "<<label(N)<<" from divisor " << label(D) << endl;

  newforms olddata(D, nf->modulus, nf->verbose, nf->characteristic);
  olddata.read_from_file_or_find();
  int old1ds=olddata.n1ds, old2ds=olddata.n2ds;
  noldclasses += old1ds;

  if (old1ds==0 && old2ds==0)
    {
      if(nf->verbose)
        cout<<" oldspace dimension is 0"<<endl;
      return;
    }

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
      INT CMD = olddata.nflist[i].CMD;
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
      cout<<" total oldspace dimension from divisor "<<label(D)<<" is "
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

// Usually, in the principal homology, the old multiplicity (i.e. the
// dimension of the oldspace) at level N from a new eigensystem at
// level D is the number of divisors of N/D, but if the class number
// is even and the eigensystem is self-twist by an unramified
// quadratic character chi, the multiplicity is the number of divisors
// in the kernel of chi.

// NB The preceding observation only applies to old multiplicities in
// the (principal) homology space.  In the Bianchi newform spaces the
// usual formula always holds, i.e. the oldform multiplcity is the
// total number of divisors, whether or not the form is self-twist.

// Return the oldspace dimension at level N of a new eigensystem at
// level D which is self-twist by genus character with discriminant d
int old_multiplicity(const Qideal& D, INT d, const Qideal& N)
{
  Qideal M = N/D;
  vector<Qideal> divisors = alldivs(M);
  return old_multiplicity(d, divisors);
}

// The same with the list of divisors of N/D given
int old_multiplicity(INT d, vector<Qideal>& divisors)
{
  int m = 0;
  for_each(divisors.begin(), divisors.end(),
           [d, &m](Qideal D)
           {m += int(D.genus_character(d)==+1);}
           );
  return m;
}

// Given a list of the new homology dimensions at level D (indexed by
// self-twist genus character), and a multiple N of D, return the old
// homology dimensions (similarly indexed) at level N.
vector<int> old_multiplicities(const Qideal& D, vector<int> newdimsD, const Qideal& N)
{
  Qideal M = N/D;
  vector<Qideal> divisors = alldivs(M);
  return old_multiplicities(newdimsD, divisors);
}

// The same with the list of divisors of N/D given
vector<int> old_multiplicities(vector<int> newdimsD, vector<Qideal>& divisors)
{
  vector<int> ans(newdimsD.size());
  // Here (d,D) runs over pairs where d is a dimension in the list
  // newdimsD and D is a discriminant divisor.  Each d is multiplied
  // by the appropriate old multiplicity depending on the self-twist
  // genus character D: the number of Ds is the same as the size of
  // newdimsD.
  std::transform(newdimsD.begin(), newdimsD.end(), Quad::all_disc_factors.begin(),
                 ans.begin(),
                 [&divisors](int d, INT D) {return d*old_multiplicity(D, divisors);}
                 );
  return ans;
}
