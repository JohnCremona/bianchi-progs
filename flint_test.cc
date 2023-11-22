#include "arith_extras.h"
#include "flint.h"

int main ()
{
  INT a(-13), b(3), c;
  cout << "a = "<<a<<" with sign "<<sign(a)<<endl;
  cout << "b = "<<b<<" with sign "<<b.sign()<<endl;
  cout << "c = "<<c<<" with sign "<<sign(c)<<endl;
  cout << "-a = "<<-a<<" with sign "<<sign(-a)<<endl;
  cout << "a.abs() = "<<a.abs()<<endl;
  cout << "abs(b) = "<<abs(b)<<endl;
  c = a+b;
  cout << "a+b = "<<c<<endl;
  cout << "a-b = "<<a-b<<endl;
  cout << "a*b = "<<a*b<<endl;
  cout << "2+a = "<<2+a<<endl;
  cout << "a+2 = "<<a+2<<endl;
  cout << "2-a = "<<2-a<<endl;
  cout << "a-2 = "<<a-2<<endl;
  cout << "2*a = "<<2*a<<endl;
  cout << "a*2 = "<<a*2<<endl;
  cout << "a%b = "<<a%b<<endl;
  cout << "mod(a,b) = "<<mod(a,b)<<endl;
  cout << "posmod(a,b) = "<<posmod(a,b)<<endl;
  cout << "gcd(a,b) = "<<gcd(a,b)<<endl;
  INT x, y;
  INT g = bezout(a,b,x,y);
  cout << "bezout(a,b,x,y) = "<<g<<" with x="<<x<<", y="<<y<<endl;
  int e=9;
  INT p = a^e;
  cout << "a^"<<e<<" = "<<p<<endl;
  long f = 3;
  cout << "b^"<<f<<" = "<<(b^f)<<endl;
  cout << "a==b? " << (a==b) << endl;
  cout << "a!=b? " << (a!=b) << endl;
  cout << "a==-13? " << (a==-13) << endl;
  cout << "a!=-13? " << (a!=-13) << endl;
  cout << "a==0? " << (a==0) << endl;
  cout << "a!=0? " << (a!=0) << endl;
  INT n = a*(a+1)*((a^2)+1);
  cout<<"Prime factors of "<<n<<" : "<<pdivs(n) << endl;
  cout<<"Divisors whose square divides "<<n<<" : "<<sqdivs(n) << endl;
  INT pr(7);
  for (long i=0; i<pr; i++)
    cout<<"("<<i<<"|"<<pr<<") = "<<legendre(INT(i),pr)<<endl;
  cout<<"Is "<<a<<" square? "<<a.is_square()<<endl;
  INT a2 = a*a;
  cout<<"Is "<<a2<<" square? "<<a2.is_square()<<endl;
  INT s = a2.isqrt();
  cout<<"(a square root is "<<s<<")"<<endl;
  cout<<"Enter an integer: "<<flush;
  cin >> a;
  cout<<"   value entered: "<<a<<endl;
  if (a.is_long())
    {
      cout<<" as a long int: "<<I2long(a)<<endl;
    }
  else
    {
      cout<<" does not fit in a long int!"<<endl;
    }
}
