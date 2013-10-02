#include <iostream>

#include <eclib/arith.h>
#include "quads.h"

int main ()
{
 int d,max;
 cout << "Enter field: " << flush;  cin >> d;
 cout << "Enter max. norm for primes: " << flush;  cin >> max;
 Quad::field(d,max);
 cout << "The field is "; Quad::displayfield(cout); cout << endl;
 Quad w(0,1);
 cout << "w   = " << w << endl;
 cout << "w*w = " << (w*w) << endl;
 Quad a,b,c,p;

 cout << "Enter a Quad, a: ";
 cin >> a;
 cout << "a = " << a << endl;
 cout << "real(a) = " << real(a) << endl;
 cout << "imag(a) = " << imag(a) << endl;
 cout << "quadconj(a) = " << quadconj(a) << endl;
 cout << "quadnorm(a) = " << quadnorm(a) << endl;
 cout << "HNF(a) = " << HNF(a) << endl;
 cout << "ideal_label(a) = " << ideal_label(a) << endl;

 cout << "Here are the first 20 primes out of "<<quadprimes.size()
      << " (one per line):\n";
 vector<Quad>::iterator pr = quadprimes.begin();
 while((pr-quadprimes.begin() < 21) && (pr!=quadprimes.end()))
   {
     p = *pr++;
     cout << p << " has norm "<<quadnorm(p)<<" and label "<<ideal_label(p)<<endl;
   }
   
 /*
 cout << "Here are all the primes:\n";
 pr = quadprimes.begin();
 while(pr!=quadprimes.end())
   cout << *pr++ << endl;
 */

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
 vector<Quad> resmoda = residues(a);  long norma=quadnorm(a);
 cout << norma << " residues modulo a: "<<resmoda<<endl;

//Test of primes and divisor functions
 cout<<"Testing primes divisor functions.\nEnter a: "; cin>>a;
 Quad pda = primdiv(a);
 cout << "A prime divisor of a: " << pda << endl;
 vector<Quad> plist = pdivs(a);
 cout << "The list of all prime divisors of a: " << plist << endl;
 vector<Quad> elist;
 pr=plist.begin();
 while(pr!=plist.end()) elist.push_back(val(*pr++,a));
 cout << "Exponents: " << elist << endl;
 vector<Quad> dlist = posdivs(a);
 cout << "The list of "<<dlist.size()<<" divisors of a (up to units): \n";
 cout << dlist << endl;
 dlist = alldivs(a);
 cout << "The list of all "<<dlist.size()<<" divisors of a: \n";
 cout << dlist << endl;
 dlist = sqdivs(a);
 cout << "The list of "<<dlist.size()<<" divisors of a whose square divides: \n";
 cout << dlist << endl;
 dlist = sqfreedivs(a);
 cout << "The list of "<<dlist.size()<<" square-free divisors of a: \n";
 cout << dlist << endl;


//Test of gcd and bezout
 Quad g,x,y;
 cout << "Testing gcd and bezout."<<endl;
 while(cout<<"Enter Quads a and b (a=0 to stop): ", cin >> a >> b, a!=0)
   {
     cout << "gcd("<<a<<","<<b<<") = "<<quadgcd(a,b)<<endl;
     g=quadbezout(a,b,x,y);
     cout << "bezout returns x="<<x<<", y="<<y<<", g="<<g<<".  ";
     if(g==a*x+b*y)cout<<"OK!"; else cout<<"WRONG!";
     cout<<endl;
   }
 cout << "Systematic testing of bezout..."<<flush;
 vector<Quad>::const_iterator ap,bp;
 for (ap=quadprimes.begin(); ap-quadprimes.begin()<10; ap++)
   for (bp=quadprimes.begin(); bp-quadprimes.begin()<10; bp++)
     {
       cout<<"quadbezout("<<*ap<<","<<*bp<<") = " <<quadbezout(*ap,*bp,x,y);
       cout<<endl;
     }
 cout<<"done.\n";

}
