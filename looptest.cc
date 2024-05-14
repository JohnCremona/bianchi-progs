#include "looper.h"

const int show_both = 0; // set to 1 to display output from standalone functions as well as Quadlooper

int main(void)
{
  long d, maxpnorm(1000);
  cerr << "Enter field: " << flush;  cin >> d;
  Quad::field(d,maxpnorm);
  Quad::displayfield(cout);
  cout << endl;

  long firstn, lastn; Quad a;
  vector<Quad> values, values2;

  cerr <<"Enter first and last norm for Quad loop: ";
  cin >> firstn >> lastn;
  cerr <<endl;
  for ( int conj : {1, 0})
    {
      values.clear();
      values2.clear();
      if(conj)
        cout << "all including conjugates, up to units" <<endl;
      else
        cout << "all excluding conjugates, up to units" <<endl;
      cout << "-------------------------------------" <<endl<<endl;

#if(1)
      values = quads_of_norm_between(INT(firstn), INT(lastn), conj, 1);
      cout << "Quads in order of norm (using Quadlooper):" << endl;

      for(Quadlooper looper2(firstn,lastn,conj); looper2.ok(); ++looper2)
        {
          a = (Quad)looper2;
          cout << "Quad = " << a << "\twith norm " << a.norm() << endl;
          values2.push_back(a);
        }
      std::sort(values2.begin(), values2.end());
      assert (values==values2);
      cout << "-------------------------------------------\n\n";
#endif
#if(1)

  cout << "Quads with each possible norm from "<<firstn<<" to "<<lastn<<":" << endl;
  Quadlooper looper(firstn,lastn,conj);
  while(looper.ok())
    {
      values = looper.values_with_current_norm(1); // sorted
      std::sort(values.begin(), values.end());
      INT n = values.front().norm();
      cout << "Norm " << n << ":\t" << values;
      if (show_both) cout << " using Quadlooper";
      cout << endl;
      values2 = quads_of_norm(n, conj, 1); // sorted
      if (show_both)
        cout << "Norm " << n << ":\t" << values2 << " using standalone function" << endl;
      assert (values==values2);
    }
  cout << "-------------------------------------------\n\n";

#endif
#if(1)
  INT n1(25), n2(50);
  looper = Quadlooper(1,0,conj); // 0 means no upper bound set
  cout << "Quads with norms up to "<<n1<<" in one list: " << endl;
  values = looper.values_with_norm_up_to(n1, 1); // sorted
  if (show_both) cout << "Using Quadlooper: ";
  cout << values << endl;
  values2 = quads_of_norm_up_to(n1, conj, 1); // sorted
  if (show_both)
    cout << "Using standalone: " << values2 << endl;
  assert (values == values2);

  cout << "Quads with norms from "<<n1+1<<" to "<<n2<<" in separate lists: " << endl;
  while(n1<n2)
    {
      values = looper.values_with_current_norm(1); //sorted
      n1 = values.front().norm();
      cout << "Norm " << n1 << ":\t" << values;
      if (show_both) cout << " using Quadlooper";
      cout << endl;
      values2 = quads_of_norm(n1, conj, 1); // sorted
      if (show_both)
        cout << "Norm " << n1 << ":\t" << values2 << " using standalone" << endl;
      assert (values == values2);
    }
  cout << "-------------------------------------------\n\n";
#endif
    }
}
