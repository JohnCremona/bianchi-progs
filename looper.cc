#include "looper.h"

void Quadlooper::setblims()
{
  INT zero(0);
  switch (d)
    {
    case 1:
      bmin = zero;
      bmax = include_conjugates? (n-1).isqrt() : (n/2).isqrt();
      break;
    case 3:
      bmin = zero;
      bmax = include_conjugates? (n-1).isqrt() : (n/3).isqrt();
      break;
    default:
      INT m = (d%4==3? 4*n: n);
      bmax = (m/d).isqrt();
      bmin = include_conjugates? -bmax : zero;
      if (d*bmin*bmin==m) bmin+=1;
    }
}

void Quadlooper::bstep()
{
  b+=1;
  if(b>bmax) {nstep(); b=bmin;}
}

void Quadlooper::nstep()
{
  n+=1;
  if (n>nlim) return;
  while(kronecker(disc,n)==-1)
    {
      n+=1;
      if (n>nlim) return;
    }
  setblims();
}

void Quadlooper::operator++()
{
  bstep();
  while(!finda()) bstep();
}

int Quadlooper::finda()
{
  if (d%4==3)
    {
      INT c, csq = 4*n-d*b*b;
      if(csq.is_square(c)) {a=(c-b)/2; return 1;}
      return 0;
    }
  return (n-d*b*b).is_square(a);
}
