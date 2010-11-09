#include <iostream>
#include "looper.h"

int main(void)
{
  int d,max=1000;
  cout << "Enter field: " << flush;  cin >> d;
  Quad::field(d,max);
  long firstn, lastn; Quad n;
  cout<<"Enter first and last norm for Quad loop: ";
  cin >> firstn >> lastn;
  
  for(Quadlooper alpha(d,firstn,lastn); alpha.ok(); ++alpha)
    {
      Quad a = (Quad)alpha;  long norma = quadnorm(a);
      cout << "Quad = " << a << "\twith norm " << norma << endl;
    }
}

