#include "looper.h"

int main(void)
{
  long d, max(1000);
  int conj;
  cout << "Enter field: " << flush;  cin >> d;
  Quad::field(d,max);
  long firstn, lastn; Quad n;

  cout << "Include conjugates? " << flush;  cin >> conj;
  cout<<"Enter first and last norm for Quad loop: ";
  cin >> firstn >> lastn;

  for(Quadlooper alpha(firstn,lastn,conj); alpha.ok(); ++alpha)
    {
      Quad a = (Quad)alpha;
      cout << "Quad = " << a << "\twith norm " << quadnorm(a) << endl;
      //cout << field_label() << " " << ideal_label(makepos(alpha)) << endl;
    }
}

