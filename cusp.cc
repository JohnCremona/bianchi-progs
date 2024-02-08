// FILE CUSP.CC

#include "cusp.h"

void cusplist::display() const
{
  int i=0;
  for(const auto& c : cusps)
    cout<<i++<<"\t"<< c <<endl;
}

//#define DEBUG_CUSPS

int cusplist::index(const RatQuad& c)
{
#ifdef DEBUG_CUSPS
  cout<<"Testing cusp "<<c<<" against existing list"<<endl;
#endif
  // add c to list if not there already, and return index
  Quad n = N.gen();
  int ic = 0;
  for(const auto& ci : cusps)
    {
#ifdef DEBUG_CUSPS
      cout<<" comparing cusp "<<c<<" with cusp "<<ci<<endl;
      int t = cuspeq_conj(ci, c, N, plusflag);
#else
      int t = cuspeq(ci, c, N, plusflag);
#endif
      if (Quad::class_number==1)
        {
          int t1 = cuspeq(ci, c, n, plusflag);
          if (t1!=t)
            {
              cout<<"*********************************************************************************************************"<<endl;
              cout<<"cuspeq("<<ci<<","<<c<<") should be "<<t1<<" but ideal version gives "<<t<<endl;
              cout<<"*********************************************************************************************************"<<endl;
            }
        }
      if (t)
        {
#ifdef DEBUG_CUSPS
          cout << "cusp "<<c<<" equivalent to cusp #"<<ic<<" ("<<cusps[ic]<<")"<<endl;
#endif
          return ic;
        }
      ic++;
    }
  // not found:
  cusps.push_back(c);
  ic = cusps.size()-1;
#ifdef DEBUG_CUSPS
  cout << "cusp "<<c<<" is new #"<<ic<<endl;
#endif
  return ic;
}

