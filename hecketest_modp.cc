// HECKETEST_MODP.CC  -- Test for Hecke operators in characteristic p>0

#include "matprocs.h"
#include "qidloop.h"
#include "homspace.h"
//#define LOOPER

#define MAXPRIME 10000

int main(void)
{
  long d, max(MAXPRIME);
 int np,ip,jp;
 Quad n; int show_mats, show_pols, show_facs, plusflag, cuspidal=1;
 cerr << "Enter field: " << flush;  cin >> d;
 if (!check_field(d))
   {
     cerr<<"field must be one of: "<<valid_fields<<endl;
     exit(1);
   }
  long ch=0;
  cerr << "Enter characteristic p (prime): " << flush;  cin >> ch;
  ZZ_p::init(ZZ(ch));
  Quad::field(d,max);
 Quad::displayfield(cout);
 cerr << "Plus space (0/1)? "; cin>>plusflag;
 cerr << "Cuspidal subspace (0/1)? "; cin>>cuspidal;
 cerr << "See the hecke matrices (0/1)? "; cin >> show_mats;
 cerr << "See the char polys (0/1)? "; cin >> show_pols;
 cerr << "Factor the char polys (0/1)? "; cin >> show_facs;
 Qideal N;
#ifdef LOOPER
 long firstn, lastn;
 cerr<<"Enter first and last norm for Quad loop: ";
 cin >> firstn >> lastn;
 cerr << "How many Hecke matrices T(P)? ";
 cin >> np;
 Qidealooper loop(firstn, lastn, 1, 1); // sorted within norm
 while( loop.not_finished() )
   {
     N = loop.next();
#else
 while(cerr<<"Enter level (ideal label or generator): ", cin>>N, !N.is_zero())
   {
#endif
  cout << ">>>> Level " << ideal_label(N) <<" = "<<gens_string(N)<<", norm = "<<N.norm()<<" <<<<" << endl;
  homspace h(N,plusflag,0, ch);  //level, plusflag, verbose, characteristic
  int dim = (cuspidal? h.h1cuspdim(): h.h1dim());
  cout << (cuspidal? "Cuspidal dimension = ": "Dimension = ") << dim << endl;

  vector<Quadprime> badprimes = h.N.factorization().sorted_primes();
  int nq = badprimes.size();
  if (dim>0)
    {
      ZZ_pX charpol;
      vector<mat_ZZ_p> tplist, wqlist;
      for ( auto& Q : badprimes)
        {
          cout << "Computing W("<<Q<<")..." << flush;
          mat_ZZ_p wq =  mat_to_mat_ZZ_p(h.calcop(AtkinLehnerQOp(Q,N),cuspidal,0,show_mats));
	  cout << "done. " << flush;
	  wqlist.push_back(wq);

          if (IsIdent(power(wq,2), dim))
            {
              cout << "Involution!" << "\n";
            }
          else
            {
              cout << "NOT an involution...." << "\n";
              exit(1);
            }
          CharPoly(charpol, wq);
          if (show_pols)
            cout << "char poly coeffs = " << charpol;
          cout << endl;

  // NB the order the factors are found is random, so we sort them

          if(show_facs)
            {
              vec_pair_ZZ_pX_long factors = berlekamp(charpol);
              // ::sort(factors.begin(), factors.end(), fact_cmp);
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
        }
#ifndef LOOPER
      cerr << "How many Hecke matrices T(P)? ";
      cin >> np;
      cout<<endl;
#endif
      ip=0;
      int ntp = 0;
      for ( auto& P : Quadprimes::list)
	{
          if (ntp>=np) break;
	  if (P.divides(N)) continue;
	  cout << "Computing T(" << P << ")..."<<flush;
	  mat_ZZ_p tp = mat_to_mat_ZZ_p(h.calcop(HeckePOp(P,N),cuspidal, 0, show_mats));
	  cout << "done. " << flush;
	  tplist.push_back(tp);

          CharPoly(charpol, tp);
          if (show_pols)
            cout << "char poly coeffs = " << charpol;
          cout << endl;

          if(show_facs)
            {
              vec_pair_ZZ_pX_long factors = berlekamp(charpol);
              // ::sort(factors.begin(), factors.end(), fact_cmp);
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

         for (int kp=0; kp<nq; kp++)
	    {
              mat_ZZ_p tpwq = tp * wqlist[kp];
              mat_ZZ_p wqtp = wqlist[kp] * tp;
	      if (tpwq!=wqtp)
                {
                  cout << "Problem: T("<<P
                       <<") and W(Q) matrix #"<<kp<<" do not commute!" << "\n";
                  exit(1);
                }
	    }
	  for (jp=0; jp<ip; jp++)
	    {
              mat_ZZ_p tp1tp2 = tp * tplist[jp];
              mat_ZZ_p tp2tp1 = tplist[jp] * tp;
	      if (tp1tp2!=tp2tp1)
		{
		  cout << "Problem: T("<<P
		       <<") does not commute with T(P) #" <<jp << "!\n";
                  exit(1);
		}
	    }
          ntp++;
	}
    }      // end of if(dim>0)

   }       // end of while()
 exit(0);
}       // end of main()
