// HECKETEST.CC  -- Test for Hecke operators

#include <eclib/mmatrix.h>
#include <NTL/mat_ZZ.h>
#include <NTL/mat_poly_ZZ.h>
//#include <NTL/pair_ZZX_long.h>
#include <NTL/ZZXFactoring.h>
#include "qidloop.h"
#include "homspace.h"
//#define LOOPER

#define MAXPRIME 10000

mat_ZZ mat_to_mat_ZZ(mat A);

ZZX scaled_charpoly(mat_ZZ A, const ZZ& den);

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
 Quad n; int show_mats, show_pols, show_factors, plusflag, cuspidal=1;
 cerr << "Enter field (one of "<<valid_fields<<"): " << flush;  cin >> d;
 if (!check_field(d))
   {
     cerr<<"field must be one of: "<<valid_fields<<endl;
     exit(1);
   }
 Quad::field(d,max);
 Quad::displayfield(cout);
 cerr << "Plus space (0/1)? "; cin>>plusflag;
 cerr << "Cuspidal subspace (0/1)? "; cin>>cuspidal;
 cerr << "See the hecke matrices (0/1)? "; cin >> show_mats;
 cerr << "See the char polys (0/1)? "; cin >> show_pols;
 cerr << "Factor the char polys (0/1)? "; cin >> show_factors;
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
  if (dim>0)
    {
      ZZX charpol; ZZ content; vec_pair_ZZX_long factors;
      vector<mat_ZZ> tplist, wqlist;
      for (pr=badprimes.begin(); pr!=badprimes.end(); ++pr)
        {
          Quadprime Q = *pr;
          cout << "Computing W_"<<Q<<"..." << flush;
          mat_ZZ wq =  mat_to_mat_ZZ(h.wop(Q,0,show_mats));
	  cout << "done. " << flush;
	  wqlist.push_back(wq);

          charpol = scaled_charpoly(wq, to_ZZ(den));
          if (show_pols)
            {
              cout << "char poly coeffs = " << charpol << endl;
            }
          if(show_factors)
            {
              factor(content, factors, charpol);
              ::sort(factors.begin(), factors.end(), fact_cmp);
              cout<<"Factors are:"<<endl;
              long nf = factors.length();
              for(int i=0; i<nf; i++)
                {
                  cout<<(i+1)<<":\t"<<factors[i].a
                      <<"\t(degree "<<deg(factors[i].a)<<")";
                  cout<<"\t to power "<<factors[i].b;
                  cout<<endl;
                }
            }
          cout << endl;
          if (IsDiag(power(wq,2), dim, to_ZZ(den*den)))
            {
              cout << "Involution!" << "\n";
            }
          else
            {
              cout << "NOT an involution...." << "\n";
              exit(1);
            }

	}
#ifndef LOOPER
      cerr << "How many Hecke matrices T_p? ";
      cin >> np;
      cout<<endl;
#endif
      ip=0;
      for (pr=Quadprimes::list.begin();
	   pr!=Quadprimes::list.end() && ((pr-Quadprimes::list.begin())<np);
	   ++pr)
	{
          Quadprime P = *pr;
	  while (P.divides(N)) {++pr; P=*pr; np++;}
	  cout << "Computing T_" << P << "..."<<flush;
	  mat_ZZ tp = mat_to_mat_ZZ(h.heckeop(P,0,show_mats));
	  cout << "done. " << flush;
	  tplist.push_back(tp);

          charpol = scaled_charpoly(tp, to_ZZ(den));
          if (show_pols)
            {
              cout << "char poly coeffs = " << charpol << endl;
            }
          if(show_factors)
            {
              factor(content, factors, charpol);
              ::sort(factors.begin(), factors.end(), fact_cmp);
              cout<<"Factors are:"<<endl;
              long nf = factors.length();
              for(int i=0; i<nf; i++)
                {
                  cout<<(i+1)<<":\t"<<factors[i].a
                      <<"\t(degree "<<deg(factors[i].a)<<")";
                  cout<<"\t to power "<<factors[i].b;
                  cout<<endl;
                }
            }
          cout << endl;

	  for (int kp=0; kp<nq; kp++)
	    {
              mat_ZZ tpwq = tp * wqlist[kp];
              mat_ZZ wqtp = wqlist[kp] * tp;
	      if (tpwq!=wqtp)
	      {
		cout << "Problem: T_"<<P
		     <<" and W_Q matrix #"<<kp<<" do not commute!" << "\n";
                exit(1);
	      }
	    }
	  for (jp=0; jp<ip; jp++)
	    {
              mat_ZZ tp1tp2 = tp * tplist[jp];
              mat_ZZ tp2tp1 = tplist[jp] * tp;
	      if (tp1tp2!=tp2tp1)
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

mat_ZZ mat_to_mat_ZZ(mat A)
{
  int i, j, d = A.nrows();

  // copy into an NTL matrix:
  mat_ZZ ntl_A;
  ntl_A.SetDims(d,d);
  for(i=1; i<=d; i++)
    for(j=1; j<=d; j++)
      ntl_A(i,j)=conv<ZZ>(A(i,j));
  return ntl_A;
}

ZZX scaled_charpoly(mat_ZZ A, const ZZ& den)
{
  ZZX charpol;
  CharPoly(charpol, A);
  long d = deg(charpol);
  if (den>1)
    {
      bigint dpow(1);
      for(int i=0; i<=d; i++)
        {
          SetCoeff(charpol, d-i, coeff(charpol, d-i)/dpow);
          dpow *= den;
        }
    }
  return charpol;
}
