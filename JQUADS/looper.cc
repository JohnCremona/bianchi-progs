// FILE looper.cc
// modified by JSB to allow underlying Field to be non-Euclidean,
// and to support options on conjugates and associates

//#include <builtin.h>
#include "arith.h"
#include "looper.h"

int issquare(long asq, long& a)
{
  if(asq<0) return 0;
  a=(long)(sqrt(asq)+0.1);
  return (asq==a*a);
}

Quadlooper::operator Quad() const
{
  Quad ans = Quad(a,b);
  if (conjflag) { ans = ans.conj(); }
  if (nretd>0)
    { return ans*Quadunits::u(nretd); }
  else
    { return ans; }
}

void Quadlooper::initbloop()
{
  while(kronecker(n,Field::d)==-1) n++;  // no chance of sol if n not sq mod d
  if (Field::t)
    { if (Field::d==3) bsqlim=n/3; else bsqlim=(4*n)/Field::d; }
  else
    { if (Field::d==1) bsqlim=n/2; else bsqlim=n/Field::d; }
  b=0;
}

void Quadlooper::bstep()
{
  b++;
  if(b*b>bsqlim) { n++; initbloop(); }
}

void Quadlooper::operator++(int)
{
  if ( (conjflag!=both)      // supposed to try conjugate next ...
      && (b!=0)              // ... and conjugate *is* different ...
      && ( (all==0) ||       // ... and conjugate won't arise as associate
	   ( (all==1) && (Quad(a,b) != Field::fundunit*(Quad(a,b).conj()) ) )
	  )
      )
    { conjflag = 1; }
  else
    {
      conjflag = 0;
      if ( (all==1) && (n>0) && ((++nretd)<Field::nunits) )
	{ ; }  // increment nretd - but this is incorporated in the if-test
      else
	{
	  nretd = 0;
	  do {bstep();} while (!finda());
	}
    }
}

int Quadlooper::finda()
{
  if(Field::t)
    {
      long c, csq = 4*n-Field::d*b*b;
      if(issquare(csq,c)) {a=(c-b)/2; return 1;}
      return 0;
    }
  else
    {
      return issquare(n-Field::d*b*b,a);
    }
}

// END OF FILE looper.cc
