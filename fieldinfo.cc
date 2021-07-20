#include <iostream>

#include <eclib/arith.h>
#include "primes.h"

int main ()
{
  int d,max;
  cerr << "Enter field: " << flush;  cin >> d;
  cerr << "Enter max. norm for primes: " << flush;  cin >> max;
  cerr << endl;
  Quad::field(d,max);
  Quad::displayfield(cout);
  if (Quad::class_number>1)
    {
      Quadprimes::display(cout, max);
      exit(0);
    }
  //  Quadprimes::display(cout, 0);
  Quad p; long np, ip, e, f;
  vector<Quad>::iterator pr = quadprimes.begin();
  while((pr!=quadprimes.end()))
    {
      p = *pr++;
      vector<long> hnf = HNF(p);
      if (hnf[2]==1)
        {
          ip = np = hnf[0];
          f = 1;
        }
      else
        {
          ip = hnf[2];
          np = ip*ip;
          f = 2;
        }
      e = 1;      if (((d==1)&&(ip==2)) || (ip==d)) e=2;
      cout << np << " "
           << ideal_label(p)
           << " (" << p << ") "
           << ip << " "
           << e << " "
           << f << endl;
    }
}
