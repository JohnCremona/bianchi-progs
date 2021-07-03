// FILE qidltest.cc  -- for testing the Qidealooper

#include <iostream>
#include "qidloop.h"
#include "primes.h"

void test1(Qideal&a)
{
  long norma = a.norm();
  cout << "Ideal " << a << " (norm " << norma << ")";
  if (a.isprincipal())
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
  if (a.isprincipal()&& !b.isprincipal())
    {
      cout << "Qideal = " << a << " (norm " << norma <<")";
      cout << ", principal generator " << a.gen() << endl;
      cout << "Qideal = " << b << " (norm " << norma << ")"<<endl;
      cout << endl;
    }
}

void init()
{
  long d;
  cout << "Enter field: " << flush;  cin >> d;
  Quad::field(d);
  cout << "Initialised field" <<endl;
  cout<<"K = Q(sqrt("<<-d<<")) = Q("<<Quad::name<<"), disc(K) = -"<<Quad::disc;
  cout<<", min poly("<<Quad::name<<") = x^2";
  if(Quad::t) cout<<"-x";
  cout<<"+"<<Quad::n<<".\n";
  cout<<"\nInitializing Quadprimes"<<endl;
  Quadprimes::init(1000);
  Quadprimes::display();
  // for (vector<Quadprime>::const_iterator Pi = Quadprimes::list.begin(); Pi != Quadprimes::list.end(); Pi++)
  //   cout << (*Pi) << " = " << (Qideal)(*Pi) <<endl;
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
  Qidealooper loop(firstn,lastn,both);
  while( loop.not_finished() )
    {
      Qideal I = loop.next();
      test1(I);
    }
}

// END OF FILE qidltest.cc
