// FILE SYMBTEST.CC --- for testing Qideal version of class symbdata

#include "symb.h"   // which includes quads.h & moddata.h
#include "qidloop.h"

int verbose;

void init()
{
  long d;
  cout << "Enter field: " << flush;  cin >> d;
  Field::init(d);
  Quadprimes::init(1000);
  cout << "Verbose? "; cin >> verbose;
}

/*
void test1(Qideal&a)
{
  long norma = a.norm();
  cout << "Qideal = " << a << "\twith norm " << norma;
  if (a.isprincipal()) cout << "\t principal generator " << a.gee0();
  cout << endl;
}
*/

void test3(const Qideal&n)
{
  long normn = n.norm();
  cout << "Qideal = " << n << "\twith norm " << normn;
  symbdata sd(n);
  if(verbose)
    { cout<<endl;
      sd.display(cout, verbose);
    }
  int ok = sd.check(verbose);
  if(ok)
    { cout<<" OK\n"; }
  else 
    { cout<<" !!!!!!CHECK FAILS!!!!!\n";
      cout<<" Halting program!"<<endl;
      exit(1);  // abort program
    }
  cout << endl;
}

void singletest()
{
 Qideal alpha,zero=0;
 while(cout<<"Enter level (1 0 0 to finish): ", cin>>alpha, alpha!=zero)
   { test3(alpha); }
}

void loopedtest()
{
  long firstn, lastn; int both;
  cout<<"Enter first and last norm for Qideal loop: ";
  cin >> firstn >> lastn;
  cout<<"Use both conjugates? (0=No, 1=Yes): ";
  cin >> both;
  cout << endl;
  for(Qidealooper fraka(firstn,lastn,both); fraka.ok(); fraka++)
   { test3((Qideal)fraka); }
}

int main ()
{
  cout << endl << "TEST PROGRAM FOR QIDEAL VERSION OF SYMBDATA" << endl;
  init();
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

// END OF FILE SYMBTEST.CC
