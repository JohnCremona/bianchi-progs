// HECKETEST.CC  -- Test for Hecke operators

#include "matprocs.h"
#include "qidloop.h"
#include "newforms.h"
//#define LOOPER

#define MAXPRIME 10000

int main(void)
{
  long d, maxpnorm(MAXPRIME);
  int np, ntp;
 Quad n; int show_mats, show_pols, show_factors, plusflag, cuspidal;
 cerr << "Enter field: " << flush;  cin >> d;
 if (!check_field(d))
   {
     cerr<<"field must be one of: "<<valid_fields<<endl;
     exit(1);
   }
 Quad::field(d,maxpnorm);
 Quad::displayfield(cout);
 int n2r = Quad::class_group_2_rank;
 cerr << "Plus space (0/1)? "; cin>>plusflag;
 cerr << "Cuspidal subspace (0/1)? "; cin>>cuspidal;
 cerr << "See the hecke matrices (0/1)? "; cin >> show_mats;
 cerr << "See the char polys (0/1)? "; cin >> show_pols;
 cerr << "Factor the char polys (0/1)? "; cin >> show_factors;
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
  homspace h(N,plusflag,0);  //level, plusflag, verbose
  int dim = (cuspidal? h.h1cuspdim(): h.h1dim());
  int den = (cuspidal? h.h1cdenom(): h.h1denom());
  cout << (cuspidal?"Cuspidal dimension = ":"Dimension = ") << dim << endl;
  //if(den!=1) cout << " denominator = " << den << endl;
  long hmod = h.h1hmod();
  if(hmod)
    {
      cout << "Failed to lift basis from Z/"<<hmod<<" to Z!" << endl;
      cout << "Hence characteristic polynomials are only correct modulo "<<hmod
           <<" and their factorizations are not useful."<<endl;
    }

  if (dim>0)
    {
      vector<mat_ZZ> tplist, tpqlist, wqlist, nulist, tpwqlist;

      if (n2r)
        {
          // Compute unramified quadratic characters and check that they are involutions and commute

          vector<Qideal> t2ideals = make_nulist(N);
          mat_ZZ I = ident_mat_ZZ(dim);

          for ( auto& A : t2ideals)
            {
              cout << "Computing nu_"<< ideal_label(A) <<"..." << flush;
              mat m = h.calcop(CharOp(A, N), cuspidal, 0, 0);
              mat_ZZ nu = mat_to_mat_ZZ(m);
              cout << "done. " << endl;
              if (show_mats)
                cout << "Matrix is \n" << m << endl;

              ZZX charpol = scaled_charpoly(nu, to_ZZ(den));
              if (show_pols)
                {
                  cout << "Coefficients of characteristic polynomial are " << charpol << endl;
                }
              if(show_factors)
                {
                  display_factors(charpol);
                }
              cout << endl;
              if (!check_involution(nu,den, 1))
                {
                  exit(1);
                }
              if (!check_commute(nu, nulist))
                {
                  cout << "********* unramified character matrices do not commute with each other ***********" << endl;
                  exit(1);
                }
              nulist.push_back(nu);
            }

          for (int i=0; i<(1<<n2r); i++)
            {
              mat_ZZ A = I;
              string sgs = "";
              for (int j=0; j<n2r; j++)
                {
                  int s = bit(i,j);
                  A = A * (nulist[j] + (s? -den: den)*I);
                  sgs += (s? "-": "+");
                }
              cout << "nu-eigenspace " << sgs << ": " << flush;
              cout << rank(A) <<endl;
            }
        }

      // Compute Atkin-Lehner operators W(Q) and check that they are
      // involutions and commute with each other and with the
      // character matrices.  Restricted to Q such that the power of Q
      // dividing N has square ideal class

      vector<Quadprime> badprimes = N.factorization().sorted_primes();
      vector<Quadprime> squarebadprimes = make_squarebadprimes(N, badprimes);
      vector<Qideal> badprimepowers;
      for ( const auto& Q : badprimes)
        {
          Qideal Qe = Q;
          int e = val(Q, N);
          while(--e) Qe*=Q;
          badprimepowers.push_back(Qe);
        }

      for ( auto& Q : squarebadprimes)
        {
          int e = val(Q,N);
          Qideal Qe = Q;
          matop op;
          while (--e) Qe *= Q;
          if (Qe.is_principal()) // we compute W_Q directly
            {
              cout << "Computing W("<<Q<<")..." << flush;
              op = AtkinLehnerQOp(Q, N);
            }
          else
            if (Qe.has_square_class()) // we compute T(A,A)*W_Q
              {
                Qideal A = Qe.sqrt_coprime_to(N);
                cout << "Computing W("<<Q<<")*T(A,A) for A="<<A<<"..." << flush;
                op = AtkinLehnerQChiOp(Q, A, N);
              }
            else  // we have an odd power of an ideal with non-square ideal class and compute nothing
              continue;
          mat m = h.calcop(op,cuspidal,0,0);
          mat_ZZ wq =  mat_to_mat_ZZ(m);
	  cout << "done. " << endl;
          if (show_mats)
            cout << "Matrix is \n" << m << "\n="<< m << endl;

          ZZX charpol = scaled_charpoly(wq, to_ZZ(den));
          if (show_pols)
            {
              cout << "Coefficients of characteristic polynomial are " << charpol << endl;
            }
          if(show_factors)
            {
              display_factors(charpol);
            }
          cout << endl;
          if (!check_involution(wq,den, 1))
            {
              exit(1);
            }
          if (!check_commute(wq, nulist))
            {
              cout << "********* W(Q) does not commute with unramified character matrices ***********" << endl;
              exit(1);
            }
          if (!check_commute(wq, wqlist))
            {
              cout << "********* W(Q) matrices do not commute with each other ***********" << endl;
              exit(1);
            }
	  wqlist.push_back(wq);
	}

      // Compute Hecke operators T(P) or T(P^2) and check that they
      // commute with each other and with the character matrices and
      // with A-Ls

#ifndef LOOPER
      cerr << "How many Hecke matrices T(P)? ";
      cin >> np;
      cout<<endl;
#endif
      Quadprime P0; // when h=2 this will be the first good nonprincipal prime.
      int P0_set=0;
      ntp = 0;
      for ( auto& P : Quadprimes::list)
	{
          if (ntp>=np) break;
          mat_ZZ tp, tpq, tpwq;
	  if (P.divides(N))
            continue;
          int P_class_order=1; // will hold order of [P] mod squares
          matop op;
          if (P.is_principal())
            {
              op = HeckePOp(P, N);
              cout << "Computing T(" << P << ")..."<<flush;
            }
          else
            {
              if (P.has_square_class())
                {
                  Qideal A = P.sqrt_coprime_to(N);
                  op = HeckePChiOp(P, A, N);
                  cout << "Computing T(" << P << ")*T(A,A) for A="<<A<<"..."<<flush;
                }
              else
                {
                  P_class_order = 2;
                  if ((!P0_set) && (Quad::class_number==2))
                    {
                      P0 = P;
                      P0_set = 1;
                      //cout<<"Setting P0 to "<<P0<<endl;
                    }
                  if ((P*P).is_principal())
                    {
                      op = HeckeP2Op(P, N);
                      cout << "Computing T(" << P << "^2)..."<<flush;
                    }
                  else
                    {
                      Qideal A = P.equivalent_mod_2_coprime_to(N,1);
                      op = HeckeP2ChiOp(P, A, N);
                      cout << "Computing T(" << P << "^2)*T(A,A) for A="<<A<<"..."<<flush;
                    }
                }
            }

          mat m = h.calcop(op,cuspidal, 0, show_mats);
          tp = mat_to_mat_ZZ(m);
          cout << "done. " << endl;
          if (show_mats)
            cout << "Matrix is \n" << m <<endl;

          ntp++;
          tplist.push_back(tp);

          ZZX charpol = scaled_charpoly(tp, to_ZZ(den));
          if (show_pols)
            {
              cout << "Coefficients of characteristic polynomial are " << charpol << endl;
            }
          if(show_factors)
            {
              display_factors(charpol);
            }
          cout << endl;

          if (!check_commute(tp, nulist))
            {
              cout << "********* T(P) does not commute with unramified character matrices ***********" << endl;
              exit(1);
            }
          if (!check_commute(tp, wqlist))
            {
              cout << "********* T(P) does not commute with W(Q) matrices ***********" << endl;
              exit(1);
            }
          if (!check_commute(tp, tplist))
            {
              cout << "********* T(P) does not commute with other T(P) matrices ***********" << endl;
              exit(1);
            }

          // Computing T(P)T(P0) when P*P0 is principal (P, P0 not principal), or just of square class

          if (!(P0_set && (P0!=P) && (P_class_order==2)))
            continue;

          Qideal PP0 = P*P0;
          if (PP0.is_principal())
            {
              op = HeckePQOp(P,P0,N);
              cout << "Computing T(" << P << ") T(" << P0 << ")..."<<flush;
            }
          else
            if (PP0.has_square_class())
              {
                Qideal A = PP0.sqrt_coprime_to(N);
                op = HeckePQChiOp(P,P0,A,N);
                cout << "Computing T(" << P << ") T(" << P0 << ")*T(A,A) with A = "<<A<<"..."<<flush;
              }
            else
              continue;

          m = h.calcop(op,cuspidal, 0, 0);
          tpq = mat_to_mat_ZZ(m);
          cout << "done. " << endl;
          if (show_mats)
            cout << "Matrix is \n" << m << endl;
          tpqlist.push_back(tpq);

          charpol = scaled_charpoly(tpq, to_ZZ(den));
          if (show_pols)
            {
              cout << "Coefficients of characteristic polynomial are " << charpol << endl;
            }
          if(show_factors)
            {
              display_factors(charpol);
            }
          cout << endl;

          if (!check_commute(tpq, nulist))
            {
              cout << "********* T(PQ) does not commute with unramified character matrices ***********" << endl;
              exit(1);
            }
          if (!check_commute(tpq, wqlist))
            {
              cout << "********* T(PQ) does not commute with all W matrices ***********" << endl;
              exit(1);
            }
          if (!check_commute(tpq, tplist))
            {
              cout << "********* T(PQ) does not commute with all T matrices ***********" << endl;
              exit(1);
            }
          if (!check_commute(tpq, tpqlist))
            {
              cout << "********* T(PQ) does not commute with other T(PQ) matrices ***********" << endl;
              exit(1);
            }

          // Computing T(P)W(Q) for suitable Q^e||N

          auto qe=badprimepowers.begin();
          auto q=badprimes.begin();
          for (; q!=badprimes.end(); ++q,++qe)
            {
              Quadprime Q = *q;
              Qideal Qe = *qe;
              if ((P*Qe).is_principal())
                {
                  cout<<"Computing T("<<P<<")W("<<ideal_label(Qe)<<")..."<<flush;
                  mat m1 = h.calcop(HeckePALQOp(P,Q,N),cuspidal, 0, 0);
                  tpwq = mat_to_mat_ZZ(m1);
                  cout << "done. " << endl;
                  if (show_mats)
                    cout << "Matrix is \n" << m1 << endl;
                  tpwqlist.push_back(tpwq);
                  charpol = scaled_charpoly(tpwq, to_ZZ(den));
                  if (show_pols)
                    {
                      cout << "Coefficients of characteristic polynomial are " << charpol << endl;
                    }
                  if(show_factors)
                    {
                      display_factors(charpol);
                    }
                  cout << endl;
                  if (!check_commute(tpwq, nulist))
                    {
                      cout << "********* T(P)W(Q) does not commute with unramified character matrices ***********" << endl;
                      exit(1);
                    }
                  if (!check_commute(tpwq, wqlist))
                    {
                      cout << "********* T(P)W(Q) does not commute with all W matrices ***********" << endl;
                      exit(1);
                    }
                  if (!check_commute(tpwq, tplist))
                    {
                      cout << "********* T(P)W(Q) does not commute with all T matrices ***********" << endl;
                      exit(1);
                    }
                  if (!check_commute(tpwq, tpqlist))
                    {
                      cout << "********* T(P)W(Q) does not commute with all T(PQ) matrices ***********" << endl;
                      exit(1);
                    }
                  if (!check_commute(tpwq, tpwqlist))
                    {
                      cout << "********* T(P)W(Q) does not commute with other T(P)W(Q) matrices ***********" << endl;
                      exit(1);
                    }
                }
            }

        }
    }      // end of if(dim>0)

}       // end of while()
exit(0);
}       // end of main()
