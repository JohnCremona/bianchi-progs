#include <iostream>
#include <sstream>
#include "newforms.h"

// Implementation of oldform member functions

// List for small fields of least norm with positive dimension.  This
// must include newforms in the minus space for tests to pass!  We
// could omit this, it just saves reading a lot of almost trivial data
// files for small levels which have no newforms.

static long min_newform_level_norm[24] = {0,65,32,49, // d=-,1,2,3
                                          0,8,3,25,   // d=-,5,6,7
                                          0,0,1,9,    // d=-,-,10,11
                                          0,1,1,12,   // d=-,13,14,15
                                          0,1,0,4,    // d=-,17,-,19
                                          0,1,1,6};   // d=-,21,22,23

oldforms::oldforms(const newforms* nfs)
  : nf(nfs)
{
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
  if (Quad::d <=23 && nf->characteristic==0)
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
  int old1ds=olddata.n1ds, old2ds=olddata.n2ds;
  if (old1ds==0 && old2ds==0)
    return;

  // Compute the oldform multiplicity.

  // Usually, this is the number of divisors of N/D, but if the class
  // number is even and the form is self-twist by an unramified
  // quadratic character chi, we only count the divisors in the kernel
  // of chi.  The effect of this is that if there is at least one
  // prime P dividing NN/D with chi(P)=-1, then the multiplicity is
  // halved (rounded up).  For example, if N/D=P^e then usually the
  // multiplicity is (e+1) but if chi(P)=-1 it is ceiling((e+1)/2),
  // i.e. for e=1,2,3,4,5,... instead of the multiplicity being
  // 2,3,4,5,... respectively it is 1,2,2,3,...

  // Note that in the special case, the multiplicity may be different
  // for different newforms at the same level.

  int oldmult = 1, any_half_mults=0;
  vector <int> old_mult(old1ds,0), half_mult(old1ds,0);
  if (nf->characteristic==0)
    {
      for (auto Pi = nf->badprimes.begin(); Pi!=nf->badprimes.end(); ++Pi)
        {
          Quadprime P = *Pi;
          int fac = (1 + val(P, N) - val(P, D));
          oldmult *= fac;
          if ((fac>1) && (Quad::class_group_2_rank>0))
            {
              for(int iform=0; iform<old1ds; iform++)
                {
                  QUINT CMD = olddata.nflist[iform].CMD;
                  if (!is_zero(CMD))
                    {
                      if (nf->verbose)
                        cout<<"Form "<<(iform+1)<<" has unramified self-twist by "<<CMD<<endl;
                      if (P.genus_character(CMD)==-1)
                        {
                          if (nf->verbose)
                            cout<<"chi("<<P<<") = -1, so halving oldform multiplicity"<<endl;
                          half_mult[iform]=1;
                          any_half_mults=1;
                        }
                      else
                        {
                          if (nf->verbose)
                            cout<<"chi("<<P<<") = +1"<<endl;
                        }
                    }
                }
            }
        }
      for(int iform=0; iform<old1ds; iform++)
        {
          old_mult[iform] = (half_mult[iform]? (1+oldmult)/2: oldmult);
        }
    }

  // As we are not using W's, all oldform multiplicities are equal to
  // oldmult, we do not need to subdivide this into 2^k pieces.

  if(nf->verbose)
    {
      if (any_half_mults)
        cout<<" oldspace dimensions are "<<old_mult<<endl;
      else
        cout<<" each oldspace dimension is "<<oldmult<<endl;
    }
  olddim2+=oldmult*old2ds;

  for(int iform=0; iform<old1ds; iform++)
    {
      oldformap.push_back(olddata.nflist[iform].oldform_eigs(N));
      oldclassdims.push_back(old_mult[iform]);
      oldlevels.push_back(D);
      noldclasses++;
    }
}

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

