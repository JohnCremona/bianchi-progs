// FILE CUSP.CC

#include "cusp.h"

void cusplist::display() const
{
  vector<RatQuad>::const_iterator ci;
  int i;
  for(ci = cusps.begin(), i=0; ci !=cusps.end(); ci++, i++)
    cout<<i<<"\t"<< *ci <<endl;
}

using namespace std::placeholders;

int cusplist::index(const RatQuad& c)
{
  // adds c to list if not there already, and return index
  vector<RatQuad>::iterator ci = std::find_if(cusps.begin(), cusps.end(), std::bind(::cuspeq,_1, c, N->modulus, N->plusflag));
  if (ci==cusps.end())
    {
      cusps.push_back(c);
      return cusps.size()-1;
    }
  else
    return ci-cusps.begin();
}

RatQuad mat22::operator()(const RatQuad& q)const
{
  Quad r = num(q), s = den(q);
  apply_left(r, s);
  return RatQuad(r,s);
}

