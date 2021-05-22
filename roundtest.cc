#include <iostream>

#include <eclib/arith.h>
#include "quads.h"

int main ()
{
  long a, b, q, r;
  for(b=2; b<=8; b++)
    {
      for(a=-10; a<=10; a++)
        {
          cout << "("<<a<<","<<b<<"):\t";
          q = roundover(a,b);
          r = a-b*q;
          cout<<"--> (q,r) = ("<<q<<","<<r<<")\n";
        }
      cout<<endl;
    }
}
