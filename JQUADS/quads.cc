// FILE QUADS.CC
// This file only exists for backwards-compatibility.  It is supposed to be
// compiled to  quads.o .  Linking in quads.o is supposed to be equivalent
// to linking in  quadarith.o ,  qideal.o (if needed) and primes.o .

#ifndef MAX_CLASSNUM
#define MAX_CLASSNUM 1
#endif

#include "primes.h"
#include "field.cc"
#include "quadarith.cc"
#if MAX_CLASSNUM>1
#include "qideal.cc"
#endif
#include "primes.cc"

// END OF FILE QUADS.CC
