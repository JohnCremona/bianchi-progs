// FILE qidltest.cc  -- for testing the Qidealooper

#include <iostream>
#include "qidloop.h"
#include "mat22.h"
#include "primes.h"

void test1(Qideal& I)
{
  cout << "Ideal " << label(I) << " = " << I << " (norm " << I.norm() << ")";
  if (I.is_principal())
    cout << " is principal, with generator " << I.gen();
  else
    cout << " is not principal, generators " <<I.gens();
  Factorization F = I.factorization();
  cout << " with factorization " << F;
  Qideal J(Quad::one);
  for (int i=0; i<F.size(); i++)
    {
      Qideal Q = F.prime_power(i);
      J*=Q;
    }
  if (J!=I)
    { cerr << "Error: multiplying up gave different answer!"<<endl;
      cerr << "(" << J << " instead of "<<I<<")"<<endl;
      exit(1);
    }
  cout << endl;
}

void test2(Qideal& I)
{
  Qideal J = I.conj();
  if (I.is_principal()&& !J.is_principal())
    {
      cout << "Qideal = " << I << " (norm " << I.norm() <<")";
      cout << ", principal generator " << I.gen() << endl;
      cout << "Qideal = " << J << " (norm " << J.norm() << ")"<<endl;
      cout << endl;
    }
}

void ABmatrixtest(Qideal& I)
{
  mat22 M = I.AB_matrix();
  cout<<"AB-matrix of "<<label(I)<<" is "<<M<<endl;
}

void looptest()
{
  long bound = 50; int n2r = Quad::class_group_2_rank;
  cout << "\nIdeals of norm up to "<<bound<<" (sorted, both conjugates";
  if (n2r>0)
    cout<<",with images of ideals under unramified quadratic class group characters" << endl;
  cout << "):" << endl;
  Qidealooper loop_both(1, bound, 1, 1);
  while( loop_both.not_finished() )
    {
      Qideal I = loop_both.next();
      string g = gens_string(I);
      cout << label(I) << " = " << I << " = " << g;
      if (n2r>0)
        {
          cout <<"\t[";
          for (int i=0; i<n2r; i++)
            {
              if (i>0) cout << ", ";
              //cout << "("<<Quad::ideal_class_mod_squares(I)<<")";
              cout << Quad::unramified_character(1<<i, I);
            }
          cout << "]";
        }
      cout << endl;
    }
  cout << "\nIdeals of norm up to "<<bound<<" (sorted, only one of each conjugate pair):" << endl;
  Qidealooper loop_one(1, bound, 0, 1);
  while( loop_one.not_finished() )
    {
      Qideal I = loop_one.next();
      cout << label(I) << " = " << I << " = " << gens_string(I) << endl;
    }

  long firstn, lastn;
  cerr<<"\nEnter first and last norm for Qideal loop: ";
  cin >> firstn >> lastn;
  if (is_zero(lastn)) exit(0);
  int both=1, sorted=1;
  Qidealooper loop(firstn, lastn, both, sorted);
  while( loop.not_finished() )
    {
      Qideal I = loop.next();
      test1(I);
      if (I.norm()<=100)
        {
          residuetest(I);
          ABmatrixtest(I);
        }
    }
}

void labeltest()
{
  long firstn(1), lastn(100);
  cout<<"testing labels and label parsing for all ideal of norm from "<<firstn<<" to "<<lastn<<endl;
  int both=1, sorted=1;
  Qidealooper loop(firstn, lastn, both, sorted);
  while( loop.not_finished() )
    {
      Qideal I = loop.next();
      string s = label(I);
      Qideal J(s);
      assert (I==J);
      assert (label(J)==s);
    }
  cout<<"label test ok"<<endl;
}

void stringtest()
{
  string s;
  while (1)
    {
      cerr << "Enter an ideal label N.i or a prime label Pp or Ppa or Ppb: ";
      Qideal I;
      Quadprime P;
      cin >> s;
      switch (s.front()) {
      case '0':
        return;
      case 'P':
        P = Quadprime(s);
        cout << "prime ideal is "<<P<<" with norm " << P.norm()
             << ", and ideal label "<<label(P)<<endl;
        break;
      default:
        I = Qideal(s);
        cout << "ideal is "<<I<<" with norm " << I.norm()
             << ", index "<<I.get_index()
             << ", and label "<<label(I)<<endl;
      };
    }
}

int main(void)
{
  cout << endl << "QIDEAL TEST PROGRAM" << endl;

  long d;
  cout << "Enter field: " << flush;  cin >> d;
  Quad::field(d);
  Quad::displayfield(cout);
  Quadprimes::display(cout, 200);

  looptest();
  stringtest();
  labeltest();
}

// END OF FILE qidltest.cc
