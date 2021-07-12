// FILE qidltest.cc  -- for testing the Qidealooper

#include <iostream>
#include "qidloop.h"
#include "primes.h"

void test1(Qideal&a)
{
  long norma = a.norm();
  cout << "Ideal " << ideal_label(a) << " = " << a << " (norm " << norma << ")";
  if (a.is_principal())
    cout << " is principal with generator " << a.gen();
  else
    cout << " is not principal, generators " <<a.gens();
  prime_factn pf(a);
  cout << " with factorization: ";
  pf.display();
  Qideal test;
  for (long i=0; i<pf.num_primes(); i++)
    for (long e=0; e<pf.expo(i); e++)
      test*=pf.prime(i);
  if (test!=a)
    { cerr << "Error: multiplying up gave different answer!"<<endl;
      exit(1);
    }
  cout << endl;
}

void test2(Qideal&a)
{
  long norma = a.norm();
  Qideal b = a.conj();
  if (a.is_principal()&& !b.is_principal())
    {
      cout << "Qideal = " << a << " (norm " << norma <<")";
      cout << ", principal generator " << a.gen() << endl;
      cout << "Qideal = " << b << " (norm " << norma << ")"<<endl;
      cout << endl;
    }
}

void looptest()
{
  long firstn, lastn;
  cout<<"Enter first and last norm for Qideal loop: ";
  cin >> firstn >> lastn;
  if (lastn==0) exit(0);
  int both=1, sorted=1;
  Qidealooper loop(firstn, lastn, both, sorted);
  while( loop.not_finished() )
    {
      Qideal I = loop.next();
      cout<<flush;
      test1(I);
      cout<<flush;
    }
}

void labeltest()
{
  long firstn=1, lastn=100;
  cout<<"testing labels and label parsing for all ideal of norm from "<<firstn<<" to "<<lastn<<endl;
  int both=1, sorted=1;
  Qidealooper loop(firstn, lastn, both, sorted);
  while( loop.not_finished() )
    {
      Qideal I = loop.next();
      string s = ideal_label(I);
      Qideal J(s);
      assert (I==J);
      assert (ideal_label(J)==s);
    }
  cout<<"label test ok"<<endl;
}

void stringtest()
{
  string s;
  long N;
  while (1)
    {
      cout << "Enter an ideal label N.i : ";
      cin >> s;
      if (s[0]=='0') return;
      Qideal I(s);
      N = I.norm();
      cout << "ideal is "<<I<<" with norm " << N <<", index "<<I.get_index()<<", and label "<<ideal_label(I)<<endl;
    }
}

void show_primes()
{
  Quadprimes::display();
  for (vector<Quadprime>::iterator Pi = Quadprimes::list.begin(); Pi != Quadprimes::list.end(); Pi++)
    {
      if ((*Pi).norm()>200) break;
      cout << (*Pi) << " = " << ideal_label(*Pi) << " = " << (Qideal)(*Pi) <<endl;
    }
}

void init()
{
  long d;
  cout << "Enter field: " << flush;  cin >> d;
  Quad::field(d);
  Quad::displayfield(cout);
}

int main(void)
{
  cout << endl << "QIDEAL TEST PROGRAM" << endl;
  init();
  show_primes();
  looptest();
  stringtest();
  labeltest();
}

// END OF FILE qidltest.cc
