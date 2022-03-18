// File OLDFORMS.H

#if     !defined(_OLDFORMS_H)
#define _OLDFORMS_H      1       //flags that this file has been included

#include "primes.h"

class newforms;

class oldforms {
 public:
  Qideal N;
  int noldclasses;
  int olddim1, olddim2, olddimall;  // total dim of rational oldforms, resp. non-rational, all oldform
 private:
  const newforms* nf;
  vector< vector<long> >  oldformap;
  vector<long> oldclassdims;
  vector<Qideal> oldlevels;
  void getoldclasses(Qideal& D);
 public:
  oldforms(const newforms* nfs);
  long dimoldpart(const vector<long> aplist);
  void display(void) const;
};

// Given the new dimensions at level D and a multiple N of D, return
// the oldspace dimensions at level N
vector<int> old_multiplicities(Qideal D, vector<int> newdimsD, Qideal N);

#endif
