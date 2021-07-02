// FILE qidltest.cc  -- for testing the Qidealooper

#include <iostream>
#include "qidloop.h"
#include "primes.h"

void test1(Qideal&a)
{
  long norma = a.norm();
  cout << "Qideal = " << a << "\twith norm " << norma;
  if (a.isprincipal())
    cout << "\t principal, generator " << a.gen();
  else
    cout << "\t non-principal, generators " <<a.gens();
  cout << endl;
}

void test2(Qideal&a)
{
  long norma = a.norm();
  Qideal b = a.conj();
  if (a.isprincipal()&& !b.isprincipal())
    {
      cout << "Qideal = " << a << "\twith norm " << norma;
      cout << "\t principal generator " << a.gen() << endl;
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

  for (int ch=1; ch<4; ch++)
    {
      switch(ch)
        {
        case 1:
          cout << "1) List ideals of norm in given range, identifying any principal generators" << endl;
          break;
        case 2:
          cout << "2) Debug tool: Search for principal ideals with non-princ conjugates (!)" << endl;
          break;
        case 3:
          cout << "3) Print prime factorization, and test by multiplying up" << endl;
        }
      Qidealooper Iloop(firstn,lastn,both);
      while( Iloop.not_finished() )
        {
          Qideal I = Iloop.next();
          switch (ch) {
          case 1: test1(I); break;
          case 2: test2(I); break;
          case 3: test3(I); break;
          }
        }
    }
}

// END OF FILE qidltest.cc
