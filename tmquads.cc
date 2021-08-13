#include <iostream>

#include "eclib/marith.h"
#include "mquads.h"

int main ()
{
 int d,max;
 cout << "Enter field: " << flush;  cin >> d;
 cout << "Enter max. norm for primes: " << flush;  cin >> max;
 mQuad::field(d,max);
 cout << "The field is "; mQuad::displayfield(cout); cout << endl;
 mQuad w(0,1);
 cout << "w   = " << w << endl;
 cout << "w*w = " << (w*w) << endl;
 mQuad a,b,c;

 cout << "Enter a mQuad, a: ";
 cin >> a;
 cout << "a = " << a << endl;
 cout << "real(a) = " << real(a) << endl;
 cout << "imag(a) = " << imag(a) << endl;
 cout << "mquadconj(a) = " << mquadconj(a) << endl;
 cout << "mquadnorm(a) = " << mquadnorm(a) << endl;
//Testing output to a string:
 char* s = to_string(a);
 cout << "As a string, a = " << s << endl;
 cout << "Here are the first 20 primes as strings (one per line):\n";
 for(mQuadvar p(mquadprimes); p.index<21; ++p)
   cout << "string: " << to_string(p) << "; normal: " << p << endl;

 b=a*a;
 cout << "b=a*a=" << b << endl;
 cout << "val(a,b) = " << val(a,b) << endl;
 cout << "div(a,b) = " << div(a,b) << endl;
 cout << "div(b,a) = " << div(b,a) << endl;
 c=b/a;
 cout << "b/a (=a?) = " << c;  
 if (c==a) {cout<<" \t--OK!";} 
 else      {cout<<" \t--NO!";} 
 cout<<endl;
 c=a/3;
 cout << "a/3 = " << c << endl;
 c=a+b;
 cout << "c = a+b = " << c << endl;
 c=c-b;
 cout << "c-b (=a?)  = " << c;
 if (c==a) {cout<<" \t--OK!";} 
 else      {cout<<" \t--NO!";} 
 cout<<endl;
 cout << "pos(a) = " << pos(a) << endl;

//Test of residues
 mQuadlist resmoda = residues(a);
 cout << "Residues modulo a: "<<resmoda<<endl;

//Test of primes and divisor functions
 cout<<"Testing primes divisor functions.\nEnter a: "; cin>>a;
 mQuad pda = primdiv(a);
 cout << "A prime divisor of a: " << pda << endl;
 mQuadlist plist = pdivs(a);
 cout << "The list of all prime divisors of a: " << plist << endl;
 mQuadlist elist(plist.length);
 for (mQuadvar pr(plist); pr.ok(); ++pr) elist[pr.index]=val(pr,a);
 cout << "Exponents: " << elist << endl;
 mQuadlist dlist = posdivs(a);
 cout << "The list of "<<dlist.length<<" divisors of a (up to units): \n";
 cout << dlist << endl;
 dlist = alldivs(a);
 cout << "The list of all "<<dlist.length<<" divisors of a: \n";
 cout << dlist << endl;
 dlist = sqdivs(a);
 cout << "The list of "<<dlist.length<<" divisors of a whose square divides: \n";
 cout << dlist << endl;
 dlist = sqfreedivs(a);
 cout << "The list of "<<dlist.length<<" square-free divisors of a: \n";
 cout << dlist << endl;


//Test of gcd and bezout
 mQuad g,x,y;
 cout << "Testing gcd and bezout."<<endl;
 while(cout<<"Enter mQuads a and b: ", cin >> a >> b, a!=0)
   {
     cout << "gcd("<<a<<","<<b<<") = "<<mquadgcd(a,b)<<endl;
     g=mquadbezout(a,b,x,y);
     cout << "bezout returns x="<<x<<", y="<<y<<", g="<<g<<".  ";
     if(g==a*x+b*y)cout<<"OK!"; else cout<<"WRONG!";
     cout<<endl;
   }
 cout << "Systematic test of bezout?"<<flush;
 int ans; cin>>ans;
 if(ans)
 {
 for (mQuadvar ap(mquadprimes); ap.index<100; ++ap)
   for (mQuadvar bp(mquadprimes); bp.index<ap.index; ++bp)
     {
       g=mquadbezout(ap,bp,x,y);
       if(g!=1)cout<<"mquadbezout("<<ap<<","<<bp<<") returns "<<g<<endl;
     }
 }
}
