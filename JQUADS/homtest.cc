// FILE homtest.cc  ---  for testing the class homspace

#include "homspace.h"   // which includes quads.h & moddata.h
#include "qidloop.h"

int verbose;
int plusflag=1;


void init()
{
  long d;
  cout << "Enter field: " << flush;  cin >> d;
  Field::init(d);
  geometry::init();
  Quadprimes::init(1000);
  cout << "Verbose? "; cin >> verbose;
  cout << endl;
}

void test1(Qideal&a)
{
  long norma = a.norm();
  cout << "Qideal = " << a << "\twith norm " << norma;
  if (a.isprincipal()) cout << "\t principal generator " << a.gee0();
  cout << endl;
}

void test2(Qideal a)
{
  long norma = a.norm();
  if (Field::is_princ())
    {
      if (!a.isprincipal())
	cerr << "Error: class number one but ideal not principal (sic!)"<<endl;
      else
	cout << ">>>> Level ("<< a.gee0() <<"), norm = "<<norma<<" <<<<\t";
    }
  else
    {
      cout << ">>>> Level ("<<a<<"), norm = "<<norma<<" <<<<\t";
    }

  homspace h(a,plusflag,verbose);  //level, plusflag, verbose
  cout << "Dimension = " << h.h1dim() << endl;
}

void singletest()
{
 Qideal alpha,zero=0;
 while(cout<<"Enter level (1 0 0 to finish): ", cin>>alpha, alpha!=zero)
   { test2(alpha); }
}

void loopedtest()
{
  long firstn, lastn; int both;
  cout<<"Enter first and last norm for Qideal loop: ";
  cin >> firstn >> lastn;
  cout<<"Use both conjugates? (0=No, 1=Yes): ";
  cin >> both;
  for(Qidealooper ql(firstn,lastn,both); ql.ok(); ql++)
   { test2((Qideal)ql); }
  cout << "END OF LOOP" << endl;
}

int main ()
{
  cout << endl << "TEST PROGRAM FOR QIDEAL VERSION OF HOMSPACE" << endl;
  init();
  if (verbose)
    { Field::display(cout);
      geometry::display(cout,verbose);
    }
  long ch;
  cout << "OPTIONS:"<< endl;
  cout << "0) quit" << endl;
  cout << "1) Input own levels" << endl;
  cout << "2) Loop through a range of levels" << endl;
  do { cout << endl << "Your choice: ";  cin >> ch; }
  while ((ch<0)||(ch>2));
  switch (ch) {
    case 1: singletest(); break;
    case 2: loopedtest(); 
  }
}

// END OF FILE homtest.cc
