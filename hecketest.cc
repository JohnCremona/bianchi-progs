// HECKETEST.CC  -- Test for Hecke operators

#include <eclib/subspace.h>
#include "homspace.h"
//#define LOOPER
#ifdef LOOPER
#include "looper.h"
#endif

#define MAXPRIME 10000

int main(void)
{
 int d,max=10000;
 int np,ip,jp,nq; 
 long firstn, lastn; Quad n; int mats, pols, plusflag;
 cout << "Enter field: " << flush;  cin >> d;
 if(!((d==1)||(d==2)||(d==3)||(d==7)||(d==11)))
   {
     cout<<"field must be one of: 1, 2, 3, 7, 11!\n";
     exit(1);
   }
 Quad::field(d,max);
 Quad::displayfield(cout);
 cout << "Plus space (0/1)? "; cin>>plusflag;
 cout << "See the hecke matrices (0/1)? "; cin >> mats;
 cout << "See the char polys (0/1)? "; cin >> pols;
#ifdef LOOPER
 cout<<"Enter first and last norm for Quad loop: ";
 cin >> firstn >> lastn;
 if(firstn<2) firstn=2;
 for(Quadlooper alpha(d,firstn,lastn); alpha.ok(); alpha++)
#else
 Quad alpha, p; 
 while(cout<<"Enter level: ", cin>>alpha, alpha!=0)
#endif
{
  n = makepos((Quad)alpha);  // makepos normalizes w.r.t. units
  long normn = quadnorm(n);
  cout << ">>>> Level " << ideal_label(n) <<" = ("<<n<<"), norm = "<<normn<<" <<<<" << endl;
  homspace h(n,plusflag,0);  //level, plusflag, verbose
  int d = h.h1dim();
  int den = h.h1denom();
  cout << "Dimension = " << d << endl;
  if(den!=1) cout << "denominator = " << den << endl;

  vector<Quad> badprimes = h.plist;
  vector<Quad>::const_iterator pr;
  nq = badprimes.size();
  if (d>0)
    {
      mat id = (den*den)*idmat(int(d));
      mat wq(d);
      vector<mat> wqlist;
      for (pr=badprimes.begin(); pr!=badprimes.end(); pr++)
	{
	  Quad q=*pr;
	  cout << "Computing W("<<q<<")...  " << flush;
	  wq=h.heckeop(q,0,mats);
	  cout << "done. " << flush;
          if (pols)
            cout << "char poly coeffs = " << charpoly(wq);
          cout << endl;
	  if (wq*wq==id) cout << "Involution!" << "\n";
	  else           cout << "NOT an involution...." << "\n";
	  wqlist.push_back(wq);
	}
      cout << "How many Hecke matrices T_p (max "<<nquadprimes<<")? "; 
      cin >> np;
      mat tp(d);
      vector<mat> tplist;
      ip=0;
      for (pr=quadprimes.begin(); 
	   pr!=quadprimes.end()&&((pr-quadprimes.begin())<np); 
	   pr++)
	{
	  while (n%(*pr)==0) {pr++; np++;}
	  p=*pr;
	  cout << "Computing T_p for p = " << p << "\n";
	  tp = h.heckeop(p,0,mats);
	  cout << "done. " << flush;
          if (pols)
            cout << "char poly coeffs = " << charpoly(tp);
          cout<<endl;
	  for (int kp=0; kp<nq; kp++)
	    {
	      if (wqlist[kp]*tp!=tp*wqlist[kp])
	      {
		cout << "Problem: T_p matrix for p = "<<p
		     <<" and W_q matrix "<<kp<<" do not commute!" << "\n";
	      }
	    }
	  for (jp=0; jp<ip; jp++)
	    {
	      if (tp*tplist[jp]!=tplist[jp]*tp)
		{
		  cout << "Problem: T_p matrix for p= "<<p
		       <<" does not commute!" << "\n";
		}
	    }
	  tplist.push_back(tp);
	}
    }      // end of if(d>0)

}       // end of while()
exit(0);
}       // end of main()
