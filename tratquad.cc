#include <iostream>
#include "eclib/arith.h"
#include "ratquads.h"

int main ()
{
 int d,max;
 cout << "Enter field: " << flush;  cin >> d;
 cout << "Enter max. norm for primes: " << flush;  cin >> max;
 Quad::field(d,max);
 cout << "The field is "; Quad::displayfield(cout); cout << endl;
 Quad w(0,1);
 Quad a,b,c;

 cout << "Enter Quads, a and b: ";
 cin >> a >> b;
 cout << "a = " << a << endl;
 cout << "b = " << b << endl;
 RatQuad q(a,b);
 cout << "q = a/b = " << q << endl;
 cout << "b*q==a? " << (b*q==a) << endl;
 cout << "round(q) = " << round(q) << endl;
 cout << "reduce(q) = " << reduce(q) << endl;
 cout << "q==round(q)+reduce(q)? " << (q==round(q)+reduce(q)) << endl;
 cout << "recip(q) = " << recip(q) << endl;
}
