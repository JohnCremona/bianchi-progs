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

vector<bigint> char_poly(mat_m A, long denom=1, int show_factors=0); // using NTL

int main(void)
{
 int d,max=10000;
 int np,ip,jp,nq; 
 Quad n; int mats, pols, facs, plusflag, cuspidal=1;
 cerr << "Enter field: " << flush;  cin >> d;
 if(!((d==1)||(d==2)||(d==3)||(d==7)||(d==11)))
   {
     cerr<<"field must be one of: 1, 2, 3, 7, 11!\n";
     exit(1);
   }
 Quad::field(d,max);
 Quad::displayfield(cout);
 cerr << "Plus space (0/1)? "; cin>>plusflag;
 cerr << "Cuspidal subspace (0/1)? "; cin>>cuspidal;
 cerr << "See the hecke matrices (0/1)? "; cin >> mats;
 cerr << "See the char polys (0/1)? "; cin >> pols;
 cerr << "Factor the char polys (0/1)? "; cin >> facs;
#ifdef LOOPER
 long firstn, lastn;
 cerr<<"Enter first and last norm for Quad loop: ";
 cin >> firstn >> lastn;
 if(firstn<2) firstn=2;
 for(Quadlooper alpha(d,firstn,lastn); alpha.ok(); alpha++)
#else
 Quad alpha, p; 
 while(cerr<<"Enter level: ", cin>>alpha, alpha!=0)
#endif
{
  n = makepos((Quad)alpha);  // makepos normalizes w.r.t. units
  long normn = quadnorm(n);
  cout << ">>>> Level " << ideal_label(n) <<" = ("<<n<<"), norm = "<<normn<<" <<<<" << endl;
  homspace h(n,plusflag,cuspidal,0);  //level, plusflag, cuspidal, verbose
  int d = h.h1dim();
  int den = h.h1denom();
  int cden = h.h1cdenom();
  if (cuspidal)
    {
      cout << "Cuspidal dimension = " << d << endl;
      den=cden;
    }
  else
    {
      cout << "Dimension = " << d << endl;
    }
  if(den!=1) cout << " denominator = " << den << endl;
  long hmod = h.h1hmod();
  if(hmod)
    {
      cout << "Failed to lift basis from Z/"<<hmod<<" to Z!" << endl;
      cout << "Hence characteristic polynomials are only correct modulo "<<hmod
           <<" and their factorizations are not useful."<<endl;
    }

  vector<Quad> badprimes = h.plist;
  vector<Quad>::const_iterator pr;
  nq = badprimes.size();
  vector<bigint> charpol;
  bigint MMODULUS = to_ZZ(MODULUS);
  if (d>0)
    {
      mat_m id = (den*den)*idmat(int(d));
      mat_m wq(d), wq2;
      vector<mat_m> wqlist;
      for (pr=badprimes.begin(); pr!=badprimes.end(); pr++)
	{
	  Quad q=*pr;
	  cout << "Computing W("<<q<<")...  " << flush;
	  wq =  h.heckeop(q,0,mats);
	  cout << "done. " << flush;
          // bigint lambda = to_ZZ(den);
          // int dimplus = addscalar(wq,-lambda).nullity();
          // cout << "+1 eigenspace has dimension "<<dimplus<<endl;
          // int dimminus = addscalar(wq,lambda).nullity();
          // cout << "-1 eigenspace has dimension "<<dimminus<<endl;
          charpol = char_poly(wq, den, facs&&!hmod);
          if (pols)
            cout << "char poly coeffs = " << charpol;
          cout << endl;
          wq2 = matmulmodp(wq, wq, MMODULUS);
	  if (wq2==id) cout << "Involution!" << "\n";
	  else
            {
              if(d<20) cout << "wq^2 = " << wq2 << endl;
              cout << "NOT an involution...." << "\n";
              exit(1);
            }
	  wqlist.push_back(wq);
	}
      cerr << "How many Hecke matrices T_p (max "<<nquadprimes<<")? "; 
      cin >> np;
      mat_m tp(d), tpwq(d), wqtp(d);
      vector<mat_m> tplist;
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
          charpol = char_poly(tp, den, facs&&!hmod);
          if (pols)
            cout << "char poly coeffs = " << charpol;
          cout<<endl;
	  for (int kp=0; kp<nq; kp++)
	    {
              tpwq = matmulmodp(tp, wqlist[kp], MMODULUS);
              wqtp = matmulmodp(wqlist[kp], tp, MMODULUS);
	      if (tpwq!=wqtp)
	      {
		cout << "Problem: T_p matrix for p = "<<p
		     <<" and W_q matrix "<<kp<<" do not commute!" << "\n";
	      }
	    }
	  for (jp=0; jp<ip; jp++)
	    {
              tpwq = matmulmodp(tp, tplist[jp], MMODULUS);
              wqtp = matmulmodp(tplist[jp], tp, MMODULUS);
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

vector<bigint> char_poly(mat_m A,  long denom, int show_factors) // using NTL
{
  int i, j, d = A.nrows();

  // copy into an NTL matrix:
  mat_ZZ ntl_A;
  ntl_A.SetDims(d,d);
  for(i=1; i<=d; i++)
    for(j=1; j<=d; j++)
      ntl_A(i,j)=A(i,j);

  // compute char poly in NTL:
  ZZX ntl_cp;
  CharPoly(ntl_cp, ntl_A);

  // rescale if d>1:
  if (denom>1)
    {
      // cout<<"Before rescaling, char poly = "<<ntl_cp<<endl;
      bigint dpow = to_ZZ(1);
      for(i=0; i<=d; i++)
        {
          SetCoeff(ntl_cp, d-i, coeff(ntl_cp, d-i)/dpow);
          dpow *= denom;
        }
      // cout<<"After rescaling, char poly = "<<ntl_cp<<endl;
    }

  // convert char poly back from NTL:
  vector<bigint> cp(d+1);
  for(i=0; i<=d; i++)
    cp[i] = coeff(ntl_cp,i);

  // compute and display factorization:
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
