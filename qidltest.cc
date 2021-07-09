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

void show_primes()
{
  Quadprimes::display();
  for (vector<Quadprime>::iterator Pi = Quadprimes::list.begin(); Pi != Quadprimes::list.end(); Pi++)
    {
      if ((*Pi).norm()>200) break;
      cout << (*Pi) << " = " << ideal_label(*Pi) << " = " << (Qideal)(*Pi) <<endl;
    }
}

void class_number()
{
  long MB = floor(2*sqrt(Quad::disc)/PI);
  cout << "Minkowski bound = " << MB << endl;
  vector<Qideal> class_reps; // prime ideals representing nontrivial ideal classes
  int nclasses = 1;
  Qidealooper loop(1, MB, 1, 1);
  while( loop.not_finished() )
    {
      Qideal I = loop.next();
      if (I.is_principal()) continue;
      int new_class=1;
      for (vector<Qideal>::iterator Ji = class_reps.begin(); (Ji != class_reps.end()) && new_class; Ji++)
        if (I.is_equivalent(*Ji)) new_class = 0;
      if (new_class)
        {
          nclasses++;
          cout << I << " is in a new ideal class (#" << nclasses << ")" << endl;
          class_reps.push_back(I);
        }
    }
  cout << "Class number = " << nclasses << " with representatives " << class_reps << endl;
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
}

int main(void)
{
  cout << endl << "QIDEAL LOOP TEST PROGRAM" << endl;
  init();
  show_primes();
  class_number();
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

// END OF FILE qidltest.cc
