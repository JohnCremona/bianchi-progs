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

#include <eclib/arith.h>
#include "qideal.h"
#include <list>

class Qidealooper {

public:
  Qidealooper(long nmin, long nmax,
              int both_conjugates=0,   // default only 1 conj
              int sorted=0);           // default unsorted within each norm
  Qideal next();
  int not_finished() const {return not I_norm_n.empty();}

private:
  long n, maxn; // current norm, norm bound
  int both;     // flag for whether or not to include conjugates
  int sort;     // flag for whether or not to yield ideals of fixed norm into standard order
  std::list<Qideal> I_norm_n; // list of ideals of norm n (not a vector as we use pop_front())
  void advance(); // advance if necessary, return true if not finished
};

// END OF FILE QIDLOOP.H
