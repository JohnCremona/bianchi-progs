// FILE CUSP.CC

#include "cusp.h"

//#define DEBUG_CUSP_EQ

int cusp::eq(const cusp& c) const
{
#ifdef DEBUG_CUSP_EQ
  cout<<"Testing equivalence of cusps "<<(*this)<<" and "<<c;
  cout<<" (N="<<level::modulus<<")"<<endl;
#endif
  Quad p1 = num(*this), p2 = num(c), q1 = den(*this), q2 = den(c);
  Quad s1,r1,s2,r2,temp;
  temp=quadbezout(p1,q1,s1,r1);  s1*=q2;
  temp=quadbezout(p2,q2,s2,r2);  s2*=q1;
  Quad q3 = quadgcd(q1*q2,level::modulus);
#ifdef DEBUG_CUSP_EQ
  cout<<"s1 =  "<<s1<<", s2 = " << s2 << ", q3 = "<<q3<<endl;
#endif
  int equiv=0; Quad u=1;
  for(int i=0; (!equiv)&&(i<(Quad::nunits)); i++,  u*=fundunit)
    {equiv = div(q3,(s1-u*s2)); 
     if(!(level::plusflag)) {i++; u*=fundunit;}
   }
#ifdef DEBUG_CUSP_EQ
  cout<<"Returning "<<equiv<<endl;   
#endif
  return equiv;
}
 
int cusplist::index(const cusp& c)
{  
  // adds c to list if not there already, and return index
  int ans=-1;
  for (int i=0; (i<num) && (ans<0); i++) if (c.eq(list[i]))  ans=i;
  if (ans==-1) {list[num]=c; ans=num; num++;}
  return ans;
}

