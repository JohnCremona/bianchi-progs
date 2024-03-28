// FILE SWAN.H: declaration of Swan's algorithm functions

#if     !defined(_SWAN_H)
#define _SWAN_H      1       //flags that this file has been included

#include <iostream>
#include <set>

#include "swan_utils.h"

class SwanData {
public:
  CuspList alist;
  CuspList slist;
  H3pointList corners;

  SwanData();
  void make_sigmas();
  CuspList get_sigmas() {
    make_sigmas();
    return slist;
  }
  void make_alphas(int verbose=0);
  CuspList get_alphas(int verbose=0) {
    make_alphas(verbose);
    return alist;
  }
  H3pointList get_corners() const {
    return corners;
  }
  vector<Quad> shifts; // x+yw with x,y in {-1,0,1}
private:
  Quadlooper denom_looper; // default init: norms from 1 to oo, both conjugates
  CuspList alistx; // list of alphas + 8 integer translates
  CuspList alistF4; // sublist of those in quarter rectangle
  INT maxn; // max denom norm of alphas considered systematically
  CuspList alist_ok, alist_open; // partition of current alphas (ok=surrounded, open=not yet)
  map<RatQuad, CuspList, RatQuad_comparison> nbrs, nbrs_ok, nbrs_open;
  CuspList slistx; // list of sigmas + 8 integer translates
  H3pointList cornersx;

  // add one alpha; use covered=1 after fiding covering alohas and
  // saturating with more
  int add_one_alpha(const RatQuad& a, int covered=0, int verbose=0);
  // add next batch of alphas from denom_looper, return number added
  int add_new_alphas(int verbose=0);

  void find_covering_alphas(int verbose=0);

  // list of singular corners [s,0] on S_a (s in slist or a translate)
  H3pointList singular_corners(const RatQuad& a);

  // Find and fill corners list, replacing alist/alistF4 with sublist of alphas (/in F4) on >=3 corners
  void find_corners(int verbose=0);
  void saturate_alphas(int verbose=0);

  int is_sigma_surrounded(const RatQuad& s, int verbose=0);
  int are_sigmas_surrounded(int verbose=0);
  // test if a is singular by reducing to recatngle and comparing with slist
  int is_singular(const RatQuad& a);

  int is_alpha_surrounded(const RatQuad& a, int verbose=0);
  int are_intersection_points_covered(const RatQuad& a, const RatQuad& b, int verbose=0);
  int are_alphas_surrounded(int verbose=0);
};

#endif
