// FILE CUSP.CC

#include "cusp.h"

//#define DEBUG_CUSP_EQ

int cusplist::cuspeq(const RatQuad& c1, const RatQuad& c2) const
{
#ifdef DEBUG_CUSP_EQ
  cout<<"Testing equivalence of cusps "<<c1<<" and "<<c2;
  cout<<" (N="<<N->modulus<<")"<<endl;
#endif
  Quad p1 = num(c1), p2 = num(c2), q1 = den(c1), q2 = den(c2);
  Quad s1,r1,s2,r2,temp;
  temp=quadbezout(p1,q1,s1,r1);  s1*=q2;
  temp=quadbezout(p2,q2,s2,r2);  s2*=q1;
  Quad q3 = quadgcd(q1*q2,N->modulus);
#ifdef DEBUG_CUSP_EQ
  cout<<"s1 =  "<<s1<<", s2 = " << s2 << ", q3 = "<<q3<<endl;
#endif
  int equiv=0; Quad u=1;
  for(int i=0; (!equiv)&&(i<(Quad::nunits)); i++,  u*=fundunit)
    {equiv = div(q3,(s1-u*s2));
     if(!(N->plusflag)) {i++; u*=fundunit;}
   }
#ifdef DEBUG_CUSP_EQ
  cout<<"Returning "<<equiv<<endl;
#endif
  return equiv;
}

int cusplist::index(const RatQuad& c)
{
  // adds c to list if not there already, and return index
  int ans=-1;
  for (int i=0; (i<number) && (ans<0); i++) if (cuspeq(c,list[i]))  ans=i;
  if (ans==-1) {list[number]=c; ans=number; number++;}
  return ans;
}

