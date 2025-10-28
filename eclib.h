// FILE ARITH_EXTRAS.H: miscellaneous integer (int/long) functions

#if     !defined(_ARITH_EXTRAS_H)
#define _ARITH_EXTRAS_H      1       //flags that this file has been included

// for convenience we put all the includes from eclib here, since all
// other files include this one directly or indirectly:

#include <eclib/xsplit.h>   // which includes linalg.h
#include <eclib/polys.h>
#include <eclib/bigrat.h>
#include <eclib/unimod.h>
#include <eclib/bitspace.h>
#include <eclib/timer.h>
#include <eclib/curvesort.h> // for letter codes
#include <eclib/int.h>  // for INT wrapping fmpz
#include <eclib/frat.h> // for RAT wrapping fmpq

#endif
