#include <iostream>

#include <eclib/arith.h>
#include "intprocs.h"

void roundtest(long alim=10, long blim=8)
{
  long a, b, q, r;
  for(b=2; b<=blim; b++)
    {
      for(a=-alim; a<=alim; a++)
        {
          cout << "("<<a<<","<<b<<"):\t";
          q = rounded_division(a,b);
          r = a-b*q;
          cout<<"--> (q,r) = ("<<q<<","<<r<<")\n";
        }
      cout<<endl;
    }
}

void lctest(int dim=3, int bound=2)
{
  vector<vector<int>> lclist = all_linear_combinations(dim, bound);
  cout << "All "<< lclist.size()
       <<" primitive integer vectors of length " << dim
       << ", max entry " << bound << " (up to sign)" << endl;
  for(auto v: lclist)
    cout << v << endl;
}

int main ()
{
  roundtest();
  int dim=4, bnd=1;
  for(int d=1; d<=dim; d++)
    lctest(d,bnd);
}
