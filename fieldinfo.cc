#include <iostream>

#include "primes.h"

#define MAX_DISC 100

int main ()
{
  long f, max;
  int one_field = 0;
  vector<long> fields = valid_field_discs();
  cerr << "Enter field (0 for all): " << flush;  cin >> f;
  cerr << "Enter max. norm for primes: " << flush;  cin >> max;
  cerr << endl;
  if (f)
    {
      fields = {f};
      one_field = 1;
    }
  for (auto di = fields.begin(); di!=fields.end(); ++di)
    {
      long D = *di;
      if ((!one_field) && D>MAX_DISC)
        break;
      long d = (D%4==0? D/4: D);
      Quad::field(d,max);
      Quad::displayfield(cout, 1); // 1 means also show info on 2-part of class group
      Quadprimes::display(cout, max, 1);
      cout << "------------------------------------------------------------" << endl;
    }
}
