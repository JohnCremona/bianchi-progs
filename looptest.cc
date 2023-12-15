#include "looper.h"

int main(void)
{
  long d, max(1000);
  int conj;
  cerr << "Enter field: " << flush;  cin >> d;
  Quad::field(d,max);
  long firstn, lastn; Quad n;

  cerr << "Include conjugates? " << flush;  cin >> conj;
  cerr <<"Enter first and last norm for Quad loop: ";
  cin >> firstn >> lastn;
  cerr <<endl;

  cout << "Quads in order of norm:" << endl;
  for(Quadlooper alpha(firstn,lastn,conj); alpha.ok(); ++alpha)
    {
      Quad a = (Quad)alpha;
      cout << "Quad = " << a << "\twith norm " << quadnorm(a) << endl;
      //cout << field_label() << " " << ideal_label(makepos(alpha)) << endl;
    }
  cout << "-------------------------------------------\n\n";

  cout << "Quads with each possible norm from "<<firstn<<" to "<<lastn<<":" << endl;
  Quadlooper alpha(firstn,lastn,conj);
  while(alpha.ok())
    {
      vector<Quad> values = alpha.values_with_current_norm();
      INT n = values[0].norm();
      cout << "Norm " << n << ":\t" << values << endl;
    }
  cout << "-------------------------------------------\n\n";
}
