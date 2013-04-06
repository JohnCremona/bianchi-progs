#include <eclib/arith.h>
#include "looper.h"

int floorsqrt(int asq)
{
  if(asq<0) return 0;
  return (int)(sqrt(double(asq)+0.1));
}

int issquare(int asq, int& a)
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
    case 2:
      bmax = floorsqrt(n/2);
      bmin = include_conjugates? -bmax : 0;
      if (2*bmin*bmin==n) bmin++;
      break;
    case 3:
      bmin = 0;
      bmax = include_conjugates? floorsqrt(n-1) : floorsqrt(n/3);
      break;
    default: // cases 7, 11
      bmax = floorsqrt(4*n/d);
      bmin = include_conjugates? -bmax : 0;
      if (d*bmin*bmin==4*n) bmin++;
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
  if(d<3)
    {
      return issquare(n-d*b*b,a);
    }
  else
    {
      int c, csq = 4*n-d*b*b;
      if(issquare(csq,c)) {a=(c-b)/2; return 1;}
      return 0;
    }
}
