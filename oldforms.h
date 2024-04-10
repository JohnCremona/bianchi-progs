// File OLDFORMS.H

#if     !defined(_OLDFORMS_H)
#define _OLDFORMS_H      1       //flags that this file has been included

#include "primes.h"

class newforms;

class oldforms {
 public:
  Qideal N;
  int noldclasses;
  int olddim1, olddim2, olddimall;         // total dim of rational (resp. non-rational, all) oldforms
  vector<int> old1dims, old2dims, olddims; // partition of previous by unramified quadratic character
  int n2r, nchi;
private:
  const newforms* nf;
  vector< vector<long> >  oldformap;
  vector<long> oldclassdims;
  vector<Qideal> oldlevels;
  void getoldclasses(Qideal& D);
 public:
  explicit oldforms(const newforms* nfs);
  long dimoldpart(const vector<long> aplist);
  void display(void) const;
};

// Given the new dimensions at level D and a multiple N of D, return
// the oldspace dimensions at level N
vector<int> old_multiplicities(const Qideal& D, vector<int> newdimsD, const Qideal& N);

// The same with the list of divisors of N/D given
vector<int> old_multiplicities(vector<int> newdimsD, vector<Qideal>& divisors);

// Return the oldspace dimension at level N of a (rational) newform at
// level D which is self-twist by discriminant d
int old_multiplicity(const Qideal& D, INT d, const Qideal& N);

// The same with the list of divisors of N/D given
int old_multiplicity(INT d, vector<Qideal>& divisors);

#endif
