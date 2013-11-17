#include <iostream>

#include <eclib/arith.h>
#include "quads.h"

int main ()
{
  int d,max;
  cerr << "Enter field: " << flush;  cin >> d;
  cerr << "Enter max. norm for primes: " << flush;  cin >> max;
  cerr << endl;
  Quad::field(d,max);
  cout<<"K = Q(sqrt("<<-d<<")) = Q("<<Quad::name<<"), disc(K) = -"<<Quad::disc;
  cout<<", min poly("<<Quad::name<<") = x^2";
  if(Quad::t) cout<<"-x";
  cout<<"+"<<Quad::n<<".\n";

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
      e = 1;      if (((d==1)&&(p==2)) || (p==d)) e=2;
      cout << np << " "
           << ideal_label(p) << " "
           << ip << " "
           << e << " "
           << f << endl;
    }
}
