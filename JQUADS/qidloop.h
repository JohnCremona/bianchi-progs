// FILE QIDLOOP.H

//////////////////////////////////////////////////////////////////////
// Usage:
//     Qideal a;
//     for (Qidealooper avar(firstn,lastn,both); avar.ok(); avar++)
//         { a=(Qideal)avar;
//           ...
//         }
// for looping through Qideals of norm from firstn to lastn.
// If both=1 (default=0) both conjugates are returned (if different)
//////////////////////////////////////////////////////////////////////

#include "arith.h"
//#include "quadarith.h"
//#include "qideal.h"
#include "quads.h"

class Qidealooper {

public: 
  Qidealooper(long nn, long nmax, int bothval=0);   // default only 1 conj
  ~Qidealooper() {;}
  void operator++(int) {do {bstep();} while(!found());}
  int ok() const {return n<=nlim;}
  operator Qideal() const {return Qideal(a,b,c);}

private: 
  int both;
  long n,nlim,a,b,blim,c,f,df;
  vector<long> clist;
  vector<long>::iterator cvar;
  int found() const {return (f%a==0);}
  void bstep();
  void initbloop();
  void nstep();
};

// END OF FILE QIDLOOP.H
