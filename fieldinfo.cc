#include <iostream>

#include "primes.h"

#define MAX_DISC 100

int main ()
{
  long f, maxpnorm;
  vector<long> fields = valid_field_discs();
  cerr << "Enter field (0 for all up to "<<MAX_DISC<<", -1 for all): " << flush;  cin >> f;
  cerr << "Enter max. norm for primes: " << flush;  cin >> maxpnorm;
  cerr << endl;
  if (f>0)
    fields = {f};
  for (auto di = fields.begin(); di!=fields.end(); ++di)
    {
      long D = *di;
      if ((f==0) && (D>MAX_DISC))
        break;
      long d = (D%4==0? D/4: D);
      Quad::field(d,maxpnorm);
      Quad::displayfield(cout, 1); // 1 means also show info on 2-part of class group
      Quadprimes::display(cout, maxpnorm, 1);
      cout << "------------------------------------------------------------" << endl;
    }
}
