#include "arith.h"
#include "looper.h"

int issquare(int asq, int& a)
{
  if(asq<0) return 0;
  a=(int)(sqrt(double(asq))+0.1);
  return (asq==a*a);
}

void Quadlooper::setbsqlim()
{
  switch (d) 
    {
    case 1: 
    case 2:         bsqlim=n/2; break;
    case 3:         bsqlim=n/3; break;
    default:        bsqlim=(4*n)/d;
	  }		    
}

void Quadlooper::bstep()
{
  b++;
  if(b*b>bsqlim) {nstep(); b=0;}
}

void Quadlooper::nstep()
{
  n++;
  while(kronecker(-d,n)==-1) n++;
  setbsqlim();
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
