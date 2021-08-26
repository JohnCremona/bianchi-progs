// FILE CUSP.CC

#include "cusp.h"

void cusplist::display() const
{
  vector<RatQuad>::const_iterator ci;
  int i;
  for(ci = cusps.begin(), i=0; ci !=cusps.end(); ++ci, ++i)
    cout<<i<<"\t"<< *ci <<endl;
}

int cusplist::index(const RatQuad& c)
{
  // add c to list if not there already, and return index
  Quad n = N.gen();
  for(vector<RatQuad>::iterator ci = cusps.begin(); ci != cusps.end(); ++ci)
    {
      int t1 = cuspeq(*ci, c, n, plusflag);
      int t2 = cuspeq(*ci, c, N, plusflag);
      if (t1!=t2)
        {
          cout<<"*********************************************************************************************************"<<endl;
          cout<<"cuspeq("<<(*ci)<<","<<c<<") should be "<<t1<<" but ideal version gives "<<t2<<endl;
          cout<<"*********************************************************************************************************"<<endl;
        }
      if (t1)
        return ci-cusps.begin();
    }
  // not found:
  cusps.push_back(c);
  return cusps.size()-1;
}

