// File OLDFORMS.H

#if     !defined(_OLDFORMS_H)
#define _OLDFORMS_H      1       //flags that this file has been included

#include "moddata.h"

class eigdata {  
// Contains eigs info about a sub-level
// Constructor can be used to get eigs from a lower level to construct oldforms
//  and also to retrieve data for this level.
public:
  Quad sublevel;
  int nforms,nforms2;
  int nap;
  vector<vector<long> > eigs;
  eigdata(const Quad& m, int neigs=-1, int verbose=0);
};


class oldforms {
 public:
  int noldclasses, nap, ntp;
  int olddim1, olddim2, olddimall;  // total dim of rational oldforms, resp. non-rational, all oldform
 private:
  vector< vector<long> >  oldformap;
  vector<long> oldclassdims; 
  vector<Quad> oldlevels;
  void getoldclasses(const Quad& d, int verbose);
 public:
  oldforms(int verbose=0);   
  long dimoldpart(const vector<long> aplist);
  void display(void) const;
};

string eigfile(const Quad& d);    //returns filename for eigs at level d

#endif
