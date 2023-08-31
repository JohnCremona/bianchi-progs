#include <iostream>

#include "primes.h"

int main ()
{
  long f, max;
  vector<long> fields = valid_fields;
  cerr << "Enter field (0 for all): " << flush;  cin >> f;
  cerr << "Enter max. norm for primes: " << flush;  cin >> max;
  cerr << endl;
  if (f)
    fields = {f};
  for (auto di = fields.begin(); di!=fields.end(); ++di)
    {
      long d = *di;
      Quad::field(d,max);
      Quad::displayfield(cout, 1); // 1 means also show info on 2-part of class group
      Quadprimes::display(cout, max, 1);
    }
}
