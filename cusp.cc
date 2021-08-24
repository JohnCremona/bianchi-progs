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
  for(vector<RatQuad>::iterator ci = cusps.begin(); ci != cusps.end(); ++ci)
    {
      if (cuspeq(*ci, c, N, plusflag))
        return ci-cusps.begin();
    }
  // not found:
  cusps.push_back(c);
  return cusps.size()-1;
}

