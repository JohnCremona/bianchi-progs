// FILE tquads.cc   -- test program for quadarith, quadprimes, etc.

#include <iostream>
//#include <builtin.h>

#include "arith.h"
#include "quads.h"

void init()
{
  long d;
  cout << "Enter field: " << flush;
  cin >> d;
  Field::init(d); 
  cout << "The field is "; Field::display(); cout << endl;

  long max;
  cout << "Enter max. norm for primes: " << flush;
  cin >> max;
  Quadprimes::init(max);
  cout << "The Quadprimes summary is "; Quadprimes::display(cout); cout<<endl;
}

int confirm_test()
{
  int ans=0;
  cout << "Proceed ? (0=No, 1=Yes) ";
  do
    {
      cin >> ans;
      if ((ans==0)||(ans==1)) return ans;
      cout << "Retry: ";
    }
  while (1);
}

void test_quad_basics(Quad& a)
{Quad w(0,1);
 cout << "w   = " << w << endl;
 cout << "w*w = " << (w*w) << endl;
 
 Quad b,c;
 cout << "Enter a Quad, a: ";
 cin >> a;
 cout << "a = " << a << endl;
 cout << "real(a) = " << real(a) << endl;
 cout << "imag(a) = " << imag(a) << endl;
 cout << "a.conj() = " << a.conj() << endl;
 cout << "a.norm() = " << a.norm() << endl;
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
 cout << "a.is_pos() = " << a.is_pos() << endl;
}

//Test of residues
void test_residues(Quad& a)
{
  Quadlist resmoda = residues(a);
  long norma=a.norm();
  cout<<"Now test the function for returning the residues modulo a Quad"<<endl;
  cout << "Residues modulo a: "<<resmoda<<endl;
}

void getqideal(Qideal& f)
{ cout << "Please enter x,y,z for the Z-basis   z[x,y+a] : ";
  cin >> f;
  cout << endl;  
}

void primes_as_strings()
{
  cout<<"List the first (up to) 20 primes as strings (one per line).\n";
  if (confirm_test())
    {
      char* s; 

// The next loop uses "while" rather than "for", because with
//    for(Primevar pvar(Quadprimes::list); pvar.ok()&&(pvar.index<21); pvar++)
// the symbol "pvar" in the for loop would be unavailable to the debugger

      Primevar pvar(Quadprimes::list);
      while (pvar.ok()&&(pvar.index<21))
	{
	  Quadprime p = (Quadprime)pvar;  
	  s = to_string(p);
	  cout << "string: " << s << "; normal: " << p << endl;
	  delete s;
	  pvar++;
	}
    }
  cout << endl;
}

void test_divisor_functions()
{
  cout<<"Testing primes divisor functions.\n";

  PRIMETYPE a;
#if MAX_CLASSNUM<2
  cout<<"Please enter a Quad: ";
  cin >>a;
#else
  getqideal(a);
#endif
  cout << endl << "You entered a = " << a << endl;
  
  Quadprime pda = primdiv(a);
  cout << "The first prime divisor of a, (found by primdiv): " << pda << endl;

  Primelist plist = pdivs(a);
  cout <<"The list of all prime divisors of a (pdivs): "<< plist <<endl;

  Tlist<long> elist(plist.getlength());
  for (Primevar pr(plist); pr.ok(); pr++)
    elist[pr.index]=val((Quadprime)pr,a);
  cout << "Their exponents, found with val: " << elist << endl;

  prime_factn pp(a);
  cout << "Now the official prime factorisation (prime_factn):" << endl;
  pp.display();
  cout << endl;

  Tlist<PRIMETYPE> dlist = posdivs(a);
  cout << "The "<<dlist.getlength()<<" 'positive' divisors of a (posdivs):";
  cout << endl << dlist << endl;
  
  dlist = alldivs(a);
  cout << "The "<<dlist.getlength()<<" divisors of a (alldivs):";
  cout << endl << dlist << endl;
  
  dlist = sqdivs(a);
  cout <<"The "<<dlist.getlength()<<" divisors whose square divides (sqdivs):";
  cout << endl << dlist << endl;
  
  dlist = sqfreedivs(a);
  cout <<"The "<<dlist.getlength()<<" square-free divisors of a (sqfreedivs):";
  cout << endl << dlist << endl;
  cout <<endl;
}

//Test of gcd and bezout
void test_bezout()
{ Quad g,x,y,a,b;
 cout << "Testing gcd and bezout."<<endl;
 while(cout<<"Enter Quads a and b: ", cin >> a >> b, a!=0)
   {
     cout << "gcd("<<a<<","<<b<<") = "<<quadgcd(a,b)<<endl;
     g=quadbezout(a,b,x,y);
     cout << "bezout returns x="<<x<<", y="<<y<<", g="<<g<<".  ";
     if(g==a*x+b*y)cout<<"OK!"; else cout<<"WRONG!";
     cout<<endl;
   }
}

#if MAX_CLASSNUM<2
void test_bezout_on_primes()
{
  cout << "The next test tries bezout on all pairs of primes in the list.\n";
  if (confirm_test())
    {
      Quad g,x,y;
      for (Quadvar ap(Quadprimes::list); ap.index<100; ap++)
	for (Quadvar bp(Quadprimes::list); bp.index<ap.index; bp++)
	  {
	    g=quadbezout((Quad)ap,(Quad)bp,x,y);
	    if(g!=1)cout<<"quadbezout("<<(Quad)ap<<","<<(Quad)bp<<") returns "<<g<<endl;
	  }
    }
}

void classnum_one_tests()
{
  test_bezout_on_primes();
}

void classnum_two_tests() {;}

#else

void test_addition_on_primes()
{
  cout << "The next test tries addition of all pairs of primes in the list.\n";
  long l = Quadprimes::list.getlength();
  cout << "There are "<< l <<" stored primes, making "<<l*(l+1)/2<<" pairs.\n";
  if (confirm_test())
    {
      Qideal g,one=1;
      cout << "Systematic testing of addition of coprime ideals..."<<endl;
      for (Primevar ap(Quadprimes::list); ap.index<l; ap++)
	for (Primevar bp(Quadprimes::list); bp.index<ap.index; bp++)
	  {
	    g=(Quadprime)ap+(Quadprime)bp;
	    if(g!=one)cout<<"Qideal addition"<< (Quadprime)ap<<"+"<<(Quadprime)bp<<") returns "<<g<<endl;
	  }
    }
  cout << endl;
}

void test_qideal_basics()
{ Qideal fa;
  Qideal fb(Quad(1));
  Qideal fc(2);
  Qideal fe(1,0,1);
  Qideal ff;

  cout << endl;
  cout<<"The next test attempts some invalid Qideal constructor calls." <<endl;
  cout<<"A warning message should be generated on each invalid attempt."<<endl;
  if (confirm_test())
    {
      for (primevar pr; pr.ok()&&pr<10; pr++)
	{ for (long i=0; i<pr; i++) 
	    { ff = Qideal(pr,i,1);
	    }
	  cout <<endl;
	}
    }

  ff=Qideal(2,1,1);
  if (fa!=fb) {cout << "Error:  fa!=fb" << "\n";}
  if (fb!=fe) {cout << "Error:  fb!=fe" << "\n";}
  if (ff==fc) {cout << "Error:  ff==fc" << "\n";}
  Quad a(4);
  if (!ff.divides(a)) {cout << "Error:  " << ff << " !divides " << a << endl;} 

  getqideal(fa);
  cout << endl << "You entered: " << fa << endl;
  if (fa.divides(5)) {cout << fa << " divides 5." << endl;}
  if (fa.divides(Qideal(5))) {cout << fa << " divides " << Qideal(5)<< "." << endl;}
  ff=Qideal(Quad(0,1)); cout << ff << endl;
}

void test_qideal_addition()
{ Qideal fa(4), fb(6), fc(24), fd(5,0,3);
  cout << fa << " " << fb << " " << fc << " " << fd << endl;
  cout << fa+4 << " " << fa+Quad(0,2) <<" "<< fa+fb <<" "<< fa+fb+fc << endl;
  Qideal fe=fc+12;  cout << fe << endl;
  fe+=6;  cout << fe << endl;
  fe+=Quad(3,0); cout << fe << endl;
  fe+=Qideal(1,0,1); cout << fe << endl;
}

void list_qidealprimes()
{
  long l = Quadprimes::list.getlength();
  cout <<"List all the "<<l<<" stored primes."<<endl;
  if (confirm_test())
    {
      cout << Quadprimes::list << endl;;
      // for (Tvar<Qideal> f(Quadprimes::list); f.ok(); f++) {; //do something}
    }
  cout << endl;
}

void classnum_one_tests() {;}

void classnum_two_tests()
{
  test_qideal_basics();
  test_qideal_addition();
  list_qidealprimes();
  test_addition_on_primes();
}

#endif


/*
// for debugging Tlist - 4 mandatory member functions
class Myclass {
  int d;
public:
  Myclass(int dd=0) {cout << "ctor Myclass"<< endl; d=dd;}
  Myclass (const Myclass&m): d(m.d) { cout<<"copy ctor ";}  // optional
  ~Myclass() { cout <<"  dtor"<< endl; }
  int operator==(const Myclass&a) const {return (d==a.d);}
  friend ostream& operator<<(ostream&s, const Myclass&a)
    {s<<"("<<a.d<<")"; return s;}
};

Tlist<Myclass> fred(5);
*/

int main ()
{ init();
  Quad a;
  test_quad_basics(a);
  test_residues(a);
  test_bezout();
  //  primes_as_strings();
  
  cout << "Tests of the divisor functions." << endl;
  if (confirm_test()) {test_divisor_functions();}

//  classnum_one_tests();
  classnum_two_tests();
}

// END OF FILE tquads.cc
