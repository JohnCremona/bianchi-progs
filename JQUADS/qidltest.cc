// FILE qidltest.cc  -- for testing the Qidealooper

#include <iostream>
#include "qidloop.h"

void test1(Qideal&a)
{
  long norma = a.norm();
  cout << "Qideal = " << a << "\twith norm " << norma;
  if (a.isprincipal()) cout << "\t principal generator " << a.gee0();
  cout << endl;
}

void test2(Qideal&a)
{
  long norma = a.norm();
  Qideal b = a.conj();
  if (a.isprincipal()&& !b.isprincipal())
    {
      cout << "Qideal = " << a << "\twith norm " << norma;
      cout << "\t principal generator " << a.gee0() << endl;
      cout << "Qideal = " << b << "\twith norm " << norma << endl;
      cout << endl;
    }
}

void test3(Qideal&a)
{
  long norma = a.norm();
  cout << "Qideal = " << a << "\twith norm " << norma <<endl;
  prime_factn pf(a);
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

void init()
{
  long d;
  cout << "Enter field: " << flush;  cin >> d;
  Field::init(d);
  Quadprimes::init(1000);
}

int main(void)
{
  cout << endl << "QIDEAL LOOP TEST PROGRAM" << endl;
  init();
  long firstn, lastn; int both;
  cout<<"Enter first and last norm for Qideal loop: ";
  cin >> firstn >> lastn;
  cout<<"Use both conjugates? (0=No, 1=Yes): ";
  cin >> both;
  long ch;
//  cout << "Verbose? "; cin >> verbose;
  cout << "OPTIONS:"<< endl;
  cout << "0) quit" << endl;
  cout << "1) List ideals of norm in given range, identifying any principal generators" << endl;
  cout << "2) Debug tool: Search for principal ideals with non-princ conjugates (!)" << endl;
  cout << "3) Print prime factn, and test by multiplying up" << endl;

  do { cout << endl << "Your choice: ";  cin >> ch; }
  while ((ch<0)||(ch>3));
  if (ch!=0)
  for(Qidealooper alpha(firstn,lastn,both); alpha.ok(); alpha++)
    {
      Qideal a = (Qideal)alpha;
      switch (ch) {
        case 1: test1(a); break;
        case 2: test2(a); break;
        case 3: test3(a); break;
      }
    }
}

// END OF FILE qidltest.cc
