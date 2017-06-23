#include <iostream>
#include "looper.h"

int main(void)
{
  int d,max=1000, conj;
  cout << "Enter field: " << flush;  cin >> d;
  Quad::field(d,max);
  long firstn, lastn; Quad n;
  
  cout << "Include conjugates? " << flush;  cin >> conj;
  cout<<"Enter first and last norm for Quad loop: ";
  cin >> firstn >> lastn;
  
  for(Quadlooper alpha(d,firstn,lastn,conj); alpha.ok(); ++alpha)
    {
      Quad a = (Quad)alpha;  long norma = quadnorm(a);
      cout << "Quad = " << a << "\twith norm " << norma << endl;
      //cout << field_label() << " " << ideal_label(makepos(alpha)) << endl;
    }
}

