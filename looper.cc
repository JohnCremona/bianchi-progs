#include "looper.h"

QUINT floorsqrt(QUINT asq)
{
  if(asq<0) return 0;
  return (QUINT)(sqrt(double(asq)+0.1));
}

int issquare(QUINT asq, QUINT& a)
{
  if(asq<0) return 0;
  a = floorsqrt(asq);
  return (asq==a*a);
}

void Quadlooper::setblims()
{
  switch (d)
    {
    case 1:
      bmin = 0;
      bmax = include_conjugates? floorsqrt(n-1) : floorsqrt(n/2);
      break;
    case 3:
      bmin = 0;
      bmax = include_conjugates? floorsqrt(n-1) : floorsqrt(n/3);
      break;
    default:
      QUINT m = (d%4==3? 4*n: n);
      bmax = floorsqrt(m/d);
      bmin = include_conjugates? -bmax : 0;
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
  while(kronecker(-d,n)==-1)
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
      QUINT c, csq = 4*n-d*b*b;
      if(issquare(csq,c)) {a=(c-b)/2; return 1;}
      return 0;
    }
  return issquare(n-d*b*b,a);
}
