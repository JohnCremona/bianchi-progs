// FILE looptest.cc  -- for tests using the Quadlooper

#include <iostream>
#include "looper.h"

int verbose;
int both,all;

void test1(Quad&a)
{
  long norma = a.norm();
  cout << "Quad = " << a << "\twith norm " << norma << endl;
}

void test2(Quad&a)
{
  int errflag=0;
  long norma = a.norm();
  Qideal ga(a);
  Qideal gb(ga.ay(), ga.bee(), ga.cee());
  if (ga!=gb)
    {cout << "Equality error !! "<<ga<< " not equal to "<<gb<<endl;}
  if (!gb.isprincipal())
    { errflag=1;
      cout << endl;
      cout << "Quad = " << a << "\t norm " << norma;
      cout << "\tQideal = " << ga << "\t norm " << ga.norm()<< endl;
      cout << "\t *** Error:  gb not returning principal !";
      cout << endl << "EXITING!" << endl; exit(1);
    }
  else
    {
      Quad posa=a.pos_assoc();
      Quad posb=gb.gee0().pos_assoc();
      if (posa!=posb)
	{ if (!errflag)
	    { cout << "Quad = " << a << "\t norm " << norma;
	      cout << "\tQideal = " << ga << "\t norm " << ga.norm();
	      cout << endl;
	    }
	  cout << "### Error: started with "<<a<<" (positive "<<posa;
	  cout << ") but obtained "<<gb.gee0()<<" (positive "<<posb<<")";
	}
    }
//  cout << endl;
}

void test3(Quad&alpha)
{
  Qideal a=Qideal(alpha);
  if (verbose) {cout << "alpha = "<<alpha<<" generating "<<a<<endl;}
  a/=alpha;
  if (a!=1)
    { cout << "Error:  alpha = "<<alpha<<" but <alpha>/=alpha = "<<a<<endl;}
}

void test4(Quad&alpha)
{
  Qideal a(alpha);
  if (verbose) {cout << "alpha = "<<alpha<<" generating "<<a<<endl;}
  for (Quadlooper ql(1,alpha.norm()); ql.ok(); ql++)
    {
      Quad beta=(Quad)ql;
      if (verbose) { cout << "beta="<<beta<<"  ";}
      a*=beta;
      a/=beta;
      if ((!a.isprincipal())||(alpha.pos_assoc()!=a.gee0().pos_assoc()))
	{ cout << "Error:  alpha = "<<alpha << endl;
	  cout << "         beta = "<<beta << endl;
	  cout << " but (<alpha>*=beta)/=alpha = "<<a<<endl;}
    }
  if (verbose) {cout << endl;}
}

void test5(Quad&alpha)
{
  Quad g,x,y;
  if (verbose) {cout << "alpha = "<<alpha<<endl;}
  for (Quadlooper ql(1,alpha.norm(),both,all); ql.ok(); ql++)
    {
      Quad beta=(Quad)ql;
      if (verbose) { cout << "beta="<<beta<<"  ";}
      g=quadbezout(alpha,beta,x,y);

      if ((Field::not_princ() || (g.div(alpha) && g.div(beta)) ) &&(g==x*alpha+y*beta))
	{ if (verbose) cout<<"g = "<<g<<" = ("<<x<<")"<<alpha<<" + ("<<y<<")"<<beta<<endl;}  //OK
      else
	{ exit(1); }
    }
  if (verbose) {cout << endl;}
}

void init()
{
  long d;
  cout << "Enter field: " << flush;  cin >> d;
  Field::init(d);
  Quadprimes::init(1000);
  cout << "Verbose ? (0=No, 1=Yes) ";
  cin >> verbose;
}

int main(void)
{
  cout << endl << "QUAD LOOP TEST PROGRAM" << endl;
  init();

  long ch;
  cout << "OPTIONS:"<< endl;
  cout << "0) quit" << endl;
  cout << "1) List Quads of norm in given range, with their norm" << endl;
  cout << "2) Debug tool: construct principal ideal, then try to recover generator" << endl;
  cout << "3) Debug tool: compute <alpha>/alpha" << endl;
  cout << "4) Debug tool: compute (<alpha>*beta)/beta for beta of smaller norm" << endl;
  cout << "5) Debug tool: compute quadbezout(alpha,beta) for beta of smaller norm" << endl;
  do { cout << endl << "Your choice: ";  cin >> ch; }
  while ((ch<0)||(ch>5));
  if (ch!=0)
    {
      long firstn, lastn;
      cout<<"Enter first and last norm for Quad loop: ";
      cin >> firstn >> lastn;
      cout<<"Use both conjugates (where different) ? (0=No, 1=Yes) ";
      cin >> both;
      cout<<"Use all associates (where different) ? (0=No, 1=Yes) ";
      cin >> all;
      long anzahl=0;
      for(Quadlooper alpha(firstn,lastn,both,all); alpha.ok(); alpha++)
	{
	  anzahl++;
	  Quad a = (Quad)alpha;
	  switch (ch) {
	    case 1: test1(a); break;
	    case 2: test2(a); break;
	    case 3: test3(a); break;
	    case 4: test4(a); break;
	    case 5: test5(a); break;
	   }
	}
      cout << "Anzahl = " << anzahl << endl;
    }
}

// END OF FILE looptest.cc
