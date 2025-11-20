// NEWHECKE.CC  -- Hecke operators: factored characteristic polynomials on the new cuspidal subspace

#include "qidloop.h"
#include "newforms.h"

//#define LOOPER

#define MAXPRIME 10000

int main(void)
{
  scalar modulus = default_modulus<scalar>();
#if (SCALAR_OPTION==3)
  //  NextPrime(modulus, pow(to_ZZ(2),256));
  NextPrime(modulus, pow(to_ZZ(2),512));
#endif
  long d, maxpnorm(MAXPRIME);
  cerr << "Enter field: " << flush;
  cin >> d;
  if (!check_field(d))
   {
     cerr<<"field must be one of: "<<valid_fields<<endl;
     exit(1);
   }
  Quad::field(d,maxpnorm);
  Quad::displayfield(cout);

  int show_only_new_pols=1;
  cerr << "See only the new char polys (0/1)? ";
  cin >> show_only_new_pols;

  int show_only_cuspidal_pols=1;
  cerr << "See only the cuspidal char polys (0/1)? ";
  cin >> show_only_cuspidal_pols;

  int np;
  cerr << "How many Hecke matrices T(P)? ";
  cin >> np;

  Qideal N;
#ifdef LOOPER
  long firstn, lastn;
  cerr<<"Enter first and last norm for Quad loop: ";
  cin >> firstn >> lastn;
  Qidealooper loop(firstn, lastn, 1, 1); // both conjugates, sorted within norm
  while( loop.not_finished() )
    {
      N = loop.next();
#else
      while(cerr<<"Enter level (ideal label or generator): ", cin>>N, !N.is_zero())
        {
#endif
          cout << ">>>> Level " << ideal_label(N) <<" = "<<gens_string(N)<<", norm = "<<N.norm()<<" <<<<" << endl;
  homspace* h = get_homspace(N, modulus);
  int dim = h->h1dim();
  int cdim = h->h1cuspdim();
  cout << "Dimension = " << dim << endl;
  cout << "Cuspidal dimension = " << cdim << endl;
  //scalar den = h->h1cdenom();
  //if(den!=1) cout << " denominator = " << den << endl;
  ZZX charpol;

  if (dim>0)
    {
      scalar hmod = h->h1hmod();
      if(hmod!=0)
        {
          cout << "Failed to lift basis from Z/"<<hmod<<" to Z -- "
               << "homology is modulo "<<hmod<<endl;
          ZZ_p::init(to_ZZ(hmod));
        }
      int ntp = 0;
      for ( auto& P : Quadprimes::list)
	{
	  if (P.divides(N)) // || !P.is_principal())
            continue;
          ntp++;
          if (ntp>np) break;
          matop T = AutoHeckeOp(P, N);
          cout << "P = " << P << " is in ideal class " << P.ideal_class()
               << ", using T = "<<T.name()<<endl;
          if (!show_only_new_pols)
            {
              if (!show_only_cuspidal_pols)
                {
                  charpol = get_poly(N, T, 0, modulus); // cuspidal=0
                  cout << "Full characteristic polynomial: "
                       << str(charpol)
                       << endl;
                  if (deg(charpol)>0)
                    {
                      cout <<"Factors:"<<endl;
                      display_factors(charpol);
                      cout << endl;
                    }
                }
              charpol = get_poly(N, T, 1, modulus); // cuspidal=1
              cout << "Full cuspidal characteristic polynomial: "
                   << str(charpol)
                   << endl;
              if (deg(charpol)>0)
                {
                  cout <<"Factors:"<<endl;
                  display_factors(charpol);
                  cout << endl;
                }
            }
          // This is commented out because the non-cuspidal oldform
          // multiplicities are not the same as for cusp forms

          // if (!show_only_cuspidal_pols)
          //   {
          //     charpol = get_new_poly(N, T, 0, modulus); // cuspidal=0
          //     int dimnew = deg(charpol);
          //     if (ntp==1)
          //       {
          //         cout << "Newspace has dimension " << dimnew << endl;
          //       }
          //     cout << "New characteristic polynomial: "
          //          << str(charpol)
          //          << endl;
          //     if (deg(charpol)>0)
          //       {
          //         cout << "Factors:" <<endl;
          //         display_factors(charpol);
          //         cout << endl;
          //       }
          //   }
          charpol = get_new_poly(N, T, 1, modulus); // cuspidal=1
          int dimnewcusp = deg(charpol);
          if (ntp==1)
            {
              cout << "Cuspidal newspace has dimension " << dimnewcusp << endl;
            }
          cout << "New cuspidal characteristic polynomial: "
               << str(charpol)
               << endl;
          if (deg(charpol)>0)
            {
              cout << "Factors:" <<endl;
              display_factors(charpol);
              cout << endl;
            }
        }
    }      // end of if(dim>0)
}       // end of while()
exit(0);
}       // end of main()
