#include "looper.h"

INT floorsqrt(INT asq)
{
  INT r(0);
  if(asq>0)
    Iasb(r, sqrt(to_bigfloat(asq)+0.1));
  return r;
}

int issquare(INT asq, INT& a)
{
  if(asq<0) return 0;
  a = floorsqrt(asq);
  return (asq==a*a);
}

void Quadlooper::setblims()
{
  INT zero(0);
  switch (d)
    {
    case 1:
      bmin = zero;
      bmax = include_conjugates? floorsqrt(n-1) : floorsqrt(n/2);
      break;
    case 3:
      bmin = zero;
      bmax = include_conjugates? floorsqrt(n-1) : floorsqrt(n/3);
      break;
    default:
      INT m = (d%4==3? 4*n: n);
      bmax = floorsqrt(m/d);
      bmin = include_conjugates? -bmax : zero;
      if (d*bmin*bmin==m) bmin++;
    }
}

void Quadlooper::bstep()
{
  b++;
  if(b>bmax) {nstep(); b=bmin;}
}

void Quadlooper::nstep()
{
  n++;
  if (n>nlim) return;
  while(kronecker(disc,n)==-1)
    {
      n++;
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
      if(issquare(csq,c)) {a=(c-b)/2; return 1;}
      return 0;
    }
  return issquare(n-d*b*b,a);
}
