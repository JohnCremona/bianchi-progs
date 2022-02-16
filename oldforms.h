// File OLDFORMS.H

#if     !defined(_OLDFORMS_H)
#define _OLDFORMS_H      1       //flags that this file has been included

#include "primes.h"

string eigfile(const Quad& N, long p=0);    //returns filename for eigs at level N, characteristic p
string eigfile(Qideal& N, long p=0);        //returns filename for eigs at level N, characteristic p

class newforms;

class oldforms {
 public:
  int noldclasses, nap;
  int olddim1, olddim2, olddimall;  // total dim of rational oldforms, resp. non-rational, all oldform
 private:
  const newforms* nf;
  Qideal N;
  vector<Quadprime> plist;
  vector< vector<long> >  oldformap;
  vector<long> oldclassdims;
  vector<Qideal> oldlevels;
  void getoldclasses(Qideal& D);
 public:
  oldforms(const newforms* nfs);
  long dimoldpart(const vector<long> aplist);
  void display(void) const;
};

#endif
