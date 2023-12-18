#include "looper.h"

int main(void)
{
  long d, max(1000);
  int conj;
  cerr << "Enter field: " << flush;  cin >> d;
  Quad::field(d,max);
  long firstn, lastn; Quad a, n;

  cerr << "Include conjugates? " << flush;  cin >> conj;
  cerr <<"Enter first and last norm for Quad loop: ";
  cin >> firstn >> lastn;
  cerr <<endl;

  cout << "Quads in order of norm:" << endl;
  for(Quadlooper looper(firstn,lastn,conj); looper.ok(); ++looper)
    {
      a = (Quad)looper;
      cout << "Quad = " << a << "\twith norm " << quadnorm(a) << endl;
    }
  cout << "-------------------------------------------\n\n";

  cout << "Quads with each possible norm from "<<firstn<<" to "<<lastn<<":" << endl;
  Quadlooper looper(firstn,lastn,conj);
  while(looper.ok())
    {
      vector<Quad> values = looper.values_with_current_norm();
      cout << "Norm " << values[0].norm() << ":\t" << values << endl;
    }
  cout << "-------------------------------------------\n\n";

  INT n1=25, n2=50;
  looper = Quadlooper(1,0,conj); // 0 means no upper bound set
  cout << "Quads with norms up to "<<n1<<" in one list: " << endl;
  cout << looper.values_with_norm_up_to(n1) << endl;

  cout << "Quads with norms from "<<n1+1<<" to "<<n2<<" in separate lists: " << endl;
  while(n1<n2)
    {
      vector<Quad> values = looper.values_with_current_norm();
      n1 = values.front().norm();
      cout << "Norm " << n1 << ":\t" << values << endl;
    }
  cout << "-------------------------------------------\n\n";
}
