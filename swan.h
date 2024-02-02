// FILE SWAN.H: declaration of Swan's algorithm functions

#if     !defined(_SWAN_H)
#define _SWAN_H      1       //flags that this file has been included

#include <iostream>
#include <set>

#include "swan_utils.h"

// return  a saturated irredundant list of alphas, and list of sigmas, in the fundamental rectangle
pair<CuspList,CuspList> find_alphas_and_sigmas(int debug=0, int verbose=0);

#endif
