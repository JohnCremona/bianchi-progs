#include "ratquads.h"
#include "primes.h"

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
 RatQuad q(a,b, check_field(d));
 cout << "q = a/b = " << q << endl;
 cout << "b*q==a? " << (b*q==a) << endl;
 cout << "round(q) = " << q.round() << endl;
 cout << "translation_reduce(q) = " << q.translation_reduce() << endl;
 cout << "q==round(q)+translation_reduce(q)? " << (q==q.round() + q.translation_reduce()) << endl;
 cout << "recip(q) = " << q.recip() << endl;

 if (Quad::class_number==1)
   exit(0);

 if (!check_field(d))
   {
     cout << "Skipping gcd and bezout tests as this field not yet fully implemented"<<endl;
     exit(0);
   }


 Qideal I = q.ideal();
 cout << "ideal(q) = "<<I<<" = "<<I.factorization()<<endl;
 cout << "q = " << q;
 cout << (q.is_principal()? " is": " is not") << " principal"<<endl;
 Qideal J = q.denominator_ideal();
 cout << "denominator_ideal(q) = "<<J<<" = "<<J.factorization()<<endl;
 cout << "representation of q with ideal coprime to 2*3*5: ";
 q.reduce(30);
 I = q.ideal();
 cout << q << ", which has ideal "<<I<<" = "<<I.factorization();
 assert (J == q.denominator_ideal());
 cout << " (and the same denominator ideal)"<<endl;
}
