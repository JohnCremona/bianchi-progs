// FILE CUSP.CC

#include "cusp.h"

void cusplist::display() const
{
  vector<RatQuad>::const_iterator ci;
  int i;
  for(ci = cusps.begin(), i=0; ci !=cusps.end(); ++ci, ++i)
    cout<<i<<"\t"<< *ci <<endl;
}

//#define DEBUG_CUSPS

int cusplist::index(const RatQuad& c)
{
  // add c to list if not there already, and return index
  Quad n = N.gen();
  for(vector<RatQuad>::iterator ci = cusps.begin(); ci != cusps.end(); ++ci)
    {
      int t = cuspeq(*ci, c, N, plusflag);
      if (Quad::class_number==1)
        {
          int t1 = cuspeq(*ci, c, n, plusflag);
          if (t1!=t)
            {
              cout<<"*********************************************************************************************************"<<endl;
              cout<<"cuspeq("<<(*ci)<<","<<c<<") should be "<<t1<<" but ideal version gives "<<t<<endl;
              cout<<"*********************************************************************************************************"<<endl;
            }
        }
      if (t)
        {
          int nc = ci-cusps.begin();
#ifdef DEBUG_CUSPS
          cout << "cusp "<<c<<" equivalent to cusp #"<<nc<<" ("<<cusps[nc]<<")"<<endl;
#endif
          return nc;
        }
    }
  // not found:
  cusps.push_back(c);
  int nc = cusps.size()-1;
#ifdef DEBUG_CUSPS
  cout << "cusp "<<c<<" is new #"<<nc<<endl;
#endif
  return nc;
}

