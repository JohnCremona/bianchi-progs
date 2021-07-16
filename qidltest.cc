// FILE qidltest.cc  -- for testing the Qidealooper

#include <iostream>
#include "qidloop.h"
#include "primes.h"

void test1(Qideal& I)
{
  long NI = I.norm();
  cout << "Ideal " << ideal_label(I) << " = " << I << " (norm " << NI << ")";
  if (I.is_principal())
    cout << " is principal, with generator " << I.gen();
  else
    cout << " is not principal, generators " <<I.gens();
  Factorization F(I);
  cout << " with factorization " << F;
  Qideal J;
  for (int i=0; i<F.size(); i++)
    J*=F.prime_power(i);
  if (J!=I)
    { cerr << "Error: multiplying up gave different answer!"<<endl;
      cerr << "(" << J << " instead of "<<I<<")"<<endl;
      exit(1);
    }
  cout << endl;
}

void test2(Qideal& I)
{
  long NI = I.norm();
  Qideal J = I.conj();
  if (I.is_principal()&& !J.is_principal())
    {
      cout << "Qideal = " << I << " (norm " << NI <<")";
      cout << ", principal generator " << I.gen() << endl;
      cout << "Qideal = " << J << " (norm " << J.norm() << ")"<<endl;
      cout << endl;
    }
}

void residuetest(Qideal& I)
{
  for (long i = 0; i<I.norm(); i++)
    {
      // assert (i==a.numres(a.resnum(i)));
      Quad r = I.resnum(i);
      long j = I.numres(r);
      if (i!=j) cout<<i<<" --> "<<r<<" --> "<<j<<" ************"<<endl;
    }
  vector<Quad> res = I.residues();
  assert ((long)res.size()==I.norm());
  cout << I.norm() << " residues mod "<<ideal_label(I)<<": "<<res<<endl;
  if (I.norm()==1) return;

  Factorization F(I);
  long phi=1;
  for (int i=0; i<F.size(); i++)
    {
      long np = F.prime(i).norm();
      phi *= (np-1);
      for (long e=1; e<F.exponent(i); e++)
        phi *= np;
    }

  pair<vector<Quad>, vector<Quad>> invres = I.invertible_residues();
  cout << phi << " invertible residues mod "<<ideal_label(I)<<":\n";
  cout<<invres.first<<endl;
  cout << " with inverses:\n";
  cout<<invres.second<<endl;
  assert ((long)invres.first.size()==phi);
}

void ABmatrixtest(Qideal& I)
{
  mat22 M = I.AB_matrix();
  cout<<"AB-matrix of "<<ideal_label(I)<<" is "<<M<<endl;
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
      Quadprime P = *Pi;
      if (P.norm()>200) break;
      vector<Quad> gg = P.gens();
      cout << P << " = " << ideal_label(P) << " = " << (Qideal)P << " = (" << gg[0] <<","<<gg[1] << ")";
      if (P.is_principal())
        cout << " = ("<<gg[0]<<") (principal)";
      else
        cout << " (not principal)";
      cout<<endl;
    }
}

void init()
{
  long d;
  cout << "Enter field: " << flush;  cin >> d;
  Quad::field(d);
  Quad::displayfield(cout);
  if (Quad::class_number==1)
    Quadprimes::init();
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
