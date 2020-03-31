// HECKETEST.CC  -- Test for Hecke operators

#include <eclib/subspace.h>
#include <eclib/mmatrix.h>
#include <NTL/mat_ZZ.h>
#include <NTL/mat_poly_ZZ.h>
#include <NTL/ZZXFactoring.h>
#include "homspace.h"
//#define LOOPER
#ifdef LOOPER
#include "looper.h"
#endif

#define MAXPRIME 10000

vector<bigint> char_poly(mat_m A, int show_factors=0); // using NTL

int main(void)
{
 int d,max=10000;
 int np,ip,jp,nq; 
 Quad n; int mats, pols, facs, plusflag;
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
 cout << "Factor the char polys (0/1)? "; cin >> facs;
#ifdef LOOPER
 long firstn, lastn;
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
  int cden = h.h1cdenom();
  long hmod = h.h1hmod();
  cout << "Dimension = " << d << endl;
  if(den!=1) cout << "denominator = " << den << endl;
  if(cden!=1) cout << "cuspidal denominator = " << cden << endl;
  if(hmod) cout << "modulus = " << hmod << endl;

  vector<Quad> badprimes = h.plist;
  vector<Quad>::const_iterator pr;
  nq = badprimes.size();
  vector<bigint> charpol;
  if (d>0)
    {
      mat_m id = (den*den)*idmat(int(d));
      mat_m wq(d), wq2;
      vector<mat_m> wqlist;
      for (pr=badprimes.begin(); pr!=badprimes.end(); pr++)
	{
	  Quad q=*pr;
	  cout << "Computing W("<<q<<")...  " << flush;
	  wq = reduce_modp(h.heckeop(q,0,mats),MODULUS);
	  cout << "done. " << flush;
          charpol = char_poly(wq, facs);
          if (pols)
            cout << "char poly coeffs = " << charpol;
          cout << endl;
          //wq2 = reduce_modp(matmulmodp(wq,wq,MODULUS),MODULUS);
          wq2 = wq*wq;
	  if (wq2==id) cout << "Involution!" << "\n";
	  else
            {
              if(d<20) cout << "wq^2 = " << wq2 << endl;
              cout << "NOT an involution...." << "\n";
            }
	  wqlist.push_back(wq);
	}
      cout << "How many Hecke matrices T_p (max "<<nquadprimes<<")? "; 
      cin >> np;
      mat_m tp(d), tpwq, wqtp;
      vector<mat_m> tplist;
      ip=0;
      for (pr=quadprimes.begin(); 
	   pr!=quadprimes.end()&&((pr-quadprimes.begin())<np); 
	   pr++)
	{
	  while (n%(*pr)==0) {pr++; np++;}
	  p=*pr;
	  cout << "Computing T_p for p = " << p << "\n";
	  tp = reduce_modp(h.heckeop(p,0,mats),MODULUS);
	  cout << "done. " << flush;
          charpol = char_poly(tp, facs);
          if (pols)
            cout << "char poly coeffs = " << charpol;
          cout<<endl;
	  for (int kp=0; kp<nq; kp++)
	    {
              //tpwq = reduce_modp(matmulmodp(tp,wqlist[kp],MODULUS),MODULUS);
              //wqtp = reduce_modp(matmulmodp(wqlist[kp],tp,MODULUS),MODULUS);
              tpwq = tp*wqlist[kp];
              wqtp = wqlist[kp]*tp;
	      if (wqtp!=tpwq)
	      {
		cout << "Problem: T_p matrix for p = "<<p
		     <<" and W_q matrix "<<kp<<" do not commute!" << "\n";
	      }
	    }
	  for (jp=0; jp<ip; jp++)
	    {
              //tpwq = reduce_modp(matmulmodp(tp,tplist[jp],MODULUS),MODULUS);
              //wqtp = reduce_modp(matmulmodp(tplist[jp],tp,MODULUS),MODULUS);
              tpwq = tp*tplist[jp];
              wqtp = tplist[jp]*tp;
	      if (tpwq!=wqtp)
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

vector<bigint> char_poly(mat_m A, int show_factors) // using NTL
{
  int i, j, d = A.nrows();
  mat_ZZ ntl_A;
  ntl_A.SetDims(d,d);
  for(i=1; i<=d; i++)
    for(j=1; j<=d; j++)
      ntl_A(i,j)=A(i,j);
  ZZX ntl_cp;
  CharPoly(ntl_cp, ntl_A);
  vector<bigint> cp(d+1);
  for(i=0; i<=d; i++)
    cp[i] = coeff(ntl_cp,i);
  if(show_factors)
    {
      vec_pair_ZZX_long factors;
      ZZ cont;
      factor(cont,factors,ntl_cp);
      cout<<"Factors are:"<<endl;
      long nf = factors.length();
      for(i=0; i<nf; i++)
        {
          cout<<(i+1)<<":\t"<<factors[i].a
              <<"\t(degree "<<deg(factors[i].a)<<")";
          cout<<"\t to power "<<factors[i].b;
          cout<<endl;
        }
    }
  return cp;
}
