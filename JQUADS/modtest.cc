// FILE MODTEST.CC --- for testing the Qideal version of class moddata

#include "moddata.h"   // which includes quads.h
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

void test1(Qideal&a)
{
  long norma = a.norm();
  cout << "Qideal = " << a << "\twith norm " << norma;
  if (a.isprincipal()) cout << "\t principal generator " << a.gee0();
  cout << endl;
}

void test2(Qideal&a)    // construct moddata only
{
  long norma = a.norm();
  cout << ">>>> Level ("<<a<<"), norm = "<<norma<<" <<<<\t";
  moddata md(a);
  if(verbose) cout<<"\n",md.display();
/*
// obsolete bit
// relates to older version of moddata which had quotient ring as member
  int ok = (md.ring_modn)->check(verbose);
  if(ok)
    { cout<<" OK\n"; }
  else 
    { cout<<" !!!!!!CHECK FAILS!!!!!\n";
      cout<<" Halting loop!"<<endl;
      exit(1);  // abort program
    }
*/
  cout << endl;
}

void test3(const Qideal&a)    // construct moddata, quotient ring; test latter
{
  long norma = a.norm();
  cout << ">>>> Level ("<<a<<"), norm = "<<norma<<" <<<<\t";
  moddata md(a);
  if(verbose) cout<<"\n",md.display();
  quotient_ring qr(a,md.get_phi());
  int ok = qr.check(verbose);
  if(ok)
    { cout<<" OK\n"; }
  else 
    { cout<<" !!!!!!CHECK FAILS!!!!!\n";
      cout<<" Halting loop!"<<endl;
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
  for(Qidealooper ql(firstn,lastn,both); ql.ok(); ql++)
   { test3((Qideal)ql); }
}

int main ()
{
  cout << endl << "TEST PROGRAM FOR QIDEAL VERSION OF MODDATA" << endl;
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

int test_various_roundovers()
{
  long a,b;
  do {
    cout << "Input a, b : ";
    cin >> a >> b;
    cout << "a=" << a << "   b=" << b;
    cout <<"  JEC(a,b)=" << old_roundover_by_JEC(a,b);
    cout <<"  myJC(a,b)=" << roundover(a,b);
    cout <<"  JSB(a,b)=" << roundover_away_from_zero(a,b);
    cout <<"  gcc(a,b)=" << a/b << endl;
  } while (1);
}

int test_various_mods()
{
  long a,b;
  do {
    cout << "Input a, b : ";
    cin >> a >> b;
    cout << "a=" << a << "   b=" << b<<endl;
    cout <<"  mod(a,b)=" << mod(a,b)<<endl;;
    double c = fmod((double)a,(double)b);
    cout <<"  fmod(a,b)=" << c <<endl;
    double bd = ((double)b)/2;
    if (c>bd) c-=b; else if (c<=-bd) c+=b;
    cout <<"  mymod(a,b)=" << c <<endl;
  } while (1);
}

/*
int main() { test_various_mods(); }

int main() { test_various_roundovers(); }
*/

// END OF FILE MODTEST.CC
