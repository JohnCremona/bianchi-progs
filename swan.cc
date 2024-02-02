// FILE SWAN.CC: implementation of Swan's algorithm

#include "swan.h"
#include "swan_sigmas.h"
#include "swan_alphas.h"

// return  a saturated irredundant list of alphas, and list of sigmas, in the fundamental rectangle
pair<CuspList,CuspList> find_alphas_and_sigmas(int debug, int verbose)
{
  auto sigmas = singular_points();
  auto alphas = find_alphas(sigmas, debug, verbose);
  return {alphas, sigmas};
}

