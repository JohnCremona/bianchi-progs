// HECKETEST.CC  -- Test for Hecke operators

#include <eclib/mmatrix.h>
#include <NTL/mat_ZZ.h>
#include <NTL/mat_poly_ZZ.h>
#include "qidloop.h"
#include "homspace.h"
//#define LOOPER

// List of fields for which this has been implemented so far:
vector<long> fields = {1,2,3,7,11,19,43,67,163, 23, 31};

#define MAXPRIME 10000

vector<bigint> char_poly(mat_m A, long denom=1, int show_factors=0); // using NTL

// function to sort a factorization vector, first by degree of factor
// then exponent of factor then lexicographically

struct factor_comparison {
  bool operator()(pair_ZZX_long& fac1, pair_ZZX_long& fac2)
  {
    // first sort by degree of the factor
    int s = deg(fac1.a) - deg(fac2.a);
    if(s) return (s<0); // true if fac1 has smaller degree

    // then sort by exponent of the factor
    s = fac1.b - fac2.b;
    if(s) return (s<0); // true if fac1 is to a lower exponent

    // finally lexicographically compare the coefficient lists
    return std::lexicographical_compare(fac1.a.rep.begin(), fac1.a.rep.end(), fac2.a.rep.begin(), fac2.a.rep.end());
  }
}
    fact_cmp;

int main(void)
{
  long d, max(MAXPRIME);
 int np,ip,jp;
 Quad n; int mats, pols, facs, plusflag, cuspidal=1;
 cerr << "Enter field (one of "<<fields<<"): " << flush;  cin >> d;
 if (!check_field(d, fields))
   {
     cerr<<"field must be one of: "<<fields<<endl;
     exit(1);
   }
 Quad::field(d,max);
 Quad::displayfield(cout);
 cerr << "Plus space (0/1)? "; cin>>plusflag;
 cerr << "Cuspidal subspace (0/1)? "; cin>>cuspidal;
 cerr << "See the hecke matrices (0/1)? "; cin >> mats;
 cerr << "See the char polys (0/1)? "; cin >> pols;
 cerr << "Factor the char polys (0/1)? "; cin >> facs;
 Qideal N;
#ifdef LOOPER
 QUINT firstn, lastn;
 cerr<<"Enter first and last norm for Quad loop: ";
 cin >> firstn >> lastn;
 cerr << "How many Hecke matrices T_p? ";
 cin >> np;
 Qidealooper loop(firstn, lastn, 1, 1); // sorted within norm
 while( loop.not_finished() )
   {
     N = loop.next();
#else
 while(cerr<<"Enter level (ideal label or generator): ", cin>>N, !N.is_zero())
   {
#endif
  QUINT normn = N.norm();
  cout << ">>>> Level " << ideal_label(N) <<" = "<<gens_string(N)<<", norm = "<<normn<<" <<<<" << endl;
  homspace h(N,plusflag,cuspidal,0);  //level, plusflag, cuspidal, verbose
  int dim = h.h1dim();
  int den = h.h1denom();
  int cden = h.h1cdenom();
  if (cuspidal)
    {
      cout << "Cuspidal dimension = " << dim << endl;
      den=cden;
    }
  else
    {
      cout << "Dimension = " << dim << endl;
    }
  if(den!=1) cout << " denominator = " << den << endl;
  long hmod = h.h1hmod();
  if(hmod)
    {
      cout << "Failed to lift basis from Z/"<<hmod<<" to Z!" << endl;
      cout << "Hence characteristic polynomials are only correct modulo "<<hmod
           <<" and their factorizations are not useful."<<endl;
    }

  vector<Quadprime> badprimes = h.N.factorization().sorted_primes();
  vector<Quadprime>::const_iterator pr;
  int nq = badprimes.size();
  vector<bigint> charpol;
  bigint MMODULUS = to_ZZ(MODULUS);
  if (dim>0)
    {
      mat_m id = (den*den)*idmat(int(dim));
      mat_m wq(dim), wq2;
      vector<mat_m> wqlist;
      for (pr=badprimes.begin(); pr!=badprimes.end(); ++pr)
        {
          Quadprime Q = *pr;
          cout << "Computing W_"<<Q<<"..." << flush;
          wq =  h.wop(Q,0,mats);
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
              if(dim<20) cout << "wq^2 = " << wq2 << endl;
              cout << "NOT an involution...." << "\n";
              exit(1);
            }
	  wqlist.push_back(wq);
	}
#ifndef LOOPER
      cerr << "How many Hecke matrices T_p? ";
      cin >> np;
      cout<<endl;
#endif
      mat_m tp(dim), tpwq(dim), wqtp(dim);
      vector<mat_m> tplist;
      ip=0;
      for (pr=Quadprimes::list.begin();
	   pr!=Quadprimes::list.end() && ((pr-Quadprimes::list.begin())<np);
	   ++pr)
	{
          Quadprime P = *pr;
	  while (P.divides(N)) {++pr; P=*pr; np++;}
	  cout << "Computing T_" << P << "..."<<flush;
	  tp = h.heckeop(P,0,mats);
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
		cout << "Problem: T_"<<P
		     <<" and W_Q matrix #"<<kp<<" do not commute!" << "\n";
                exit(1);
	      }
	    }
	  for (jp=0; jp<ip; jp++)
	    {
              tpwq = matmulmodp(tp, tplist[jp], MMODULUS);
              wqtp = matmulmodp(tplist[jp], tp, MMODULUS);
	      if (tpwq!=wqtp)
		{
		  cout << "Problem: T_"<<P
		       <<" does not commute with T_P #" <<jp << "!\n";
                  exit(1);
		}
	    }
	  tplist.push_back(tp);
	}
    }      // end of if(dim>0)

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

  // NB the order the factors are found is random, so we sort them

  if(show_factors)
    {
      vec_pair_ZZX_long factors;
      ZZ cont;
      factor(cont,factors,ntl_cp);
      ::sort(factors.begin(), factors.end(), fact_cmp);
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
