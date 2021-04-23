// FILE QIDLOOP.CC

//#include <builtin.h>
#include "arith.h"
#include "qidloop.h"

Qidealooper::Qidealooper(long nn, long nmax, int bothval)
                         : n(nn), nlim(nmax), both(bothval)
{ 
  clist = sqdivs(n);
  cvar = clist.begin();
  initbloop();
  while(!found()) bstep();
}

void Qidealooper::bstep()
{
  if ( (++b) < blim )
    {f=(f+df)%a; df+=2;}  // without %a 1st ideal lost by overflow at a=47743
  else
    {
      cvar++;
      if (cvar==clist.end()) nstep(); // after nstep, new *cvar *is* ok ( = 1 )
      initbloop();
    }
  //  cout<<"After bstep, (a,b,c)=("<<a<<","<<b<<","<<c<<")"<<endl;
}

void Qidealooper::initbloop()
{
  c=*cvar;   
  a=n/(c*c);
  b=0;
  f=Field::n;
  df=1+(Field::t);
  if (both)
    { blim=a; }
  else
    { blim = a - Field::t;
      blim = ( blim - (blim%2) ) / 2 ;
      blim+=1;
    }
}

void Qidealooper::nstep()
{
  n++;
  clist = sqdivs(n);
  //  cout<<"n="<<n<<", clist = "<<clist<<endl;
  cvar=clist.begin();
}

// END OF FILE QIDLOOP.CC
