// FILE QIDLOOP.H

/////////////////////////////////////////////////////////////////////////////
//
// Usage:
//
//     Qidealooper looper(firstn,lastn,both,sorted);
//     while (looper.not_finished())
//         { Qideal A = looper.next();
//           ...
//         }
// for looping through Qideals of norm from firstn to lastn.
// If both=1 (default=0) both conjugates are returned (if different).
// If sorted=1 (default 0) the ideals are returned in standard sorted order.
/////////////////////////////////////////////////////////////////////////////

#if     !defined(_QIDLOOP_H)
#define _QIDLOOP_H      1       //flags that this file has been included

#include "qideal.h"
#include <list>

class Qidealooper {

public:
  Qidealooper(QUINT nmin, QUINT nmax,
              int both_conjugates=0,   // default only 1 conj
              int sorted=0);           // default unsorted within each norm
  Qideal next();
  int not_finished() const {return not I_norm_n.empty();}

private:
  QUINT n, maxn; // current norm, norm bound
  int both;     // flag for whether or not to include conjugates
  int sort;     // flag for whether or not to yield ideals of fixed norm into standard order
  std::list<Qideal> I_norm_n; // list of ideals of norm n (not a vector as we use pop_front())
  void advance(); // advance if necessary, return true if not finished
};

#endif
// END OF FILE QIDLOOP.H
