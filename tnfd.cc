// FILE TNFD.CC:  test program for nfd (d-dimensional newform) class
//////////////////////////////////////////////////////////////////////////
//
// Adapted from the similar program (over Q) in eclib
//
//////////////////////////////////////////////////////////////////////////

#include "nfd.h"

//#define DEBUG
//#define COMPARE_OLD
//#define LOOPER
#define MAXPRIME 10000

int main()
{
  cout << "Program tnfd." << endl;
  scalar modulus = default_modulus<scalar>();
#if (SCALAR_OPTION==3)
  //  NextPrime(modulus, power_ZZ(2,256));
  NextPrime(modulus, power_ZZ(2,512));
#endif
  long d, maxpnorm(MAXPRIME);
  cerr << "Enter field (one of "<<class_number_one_fields<<"): " << flush;  cin >> d;
  if (!is_class_number_one(d))
    {
      cerr<<"field must be one of: "<<class_number_one_fields<<endl;
      exit(1);
    }
  Quad::field(d,maxpnorm);
  Quad::displayfield(cout);

  int verbose=1;
  cerr << "Verbose output? (0/1) "; cin >> verbose;

  Quad n;
  Qideal N;
#ifdef LOOPER
  long firstn, lastn;
  cerr<<"Enter first and last norm for Quad loop: ";
  cin >> firstn >> lastn;
  Qidealooper loop(firstn, lastn, 1, 1); // sorted within norm
  while( loop.not_finished() )
    {
      N = loop.next();
#else
  while(cerr<<"Enter level (ideal label or generator): ", cin>>N, !N.is_zero())
    {
#endif
     cout << ">>>> Level " << ideal_label(N) <<" = "<<gens_string(N)<<", norm = "<<N.norm()<<" <<<<" << endl;
     homspace* hplus = get_homspace(N, modulus);
     // int dimH = hplus->h1cuspdim();
     // cout << "dimension = " << dimH << endl;
     // cout << "denom     = " << hplus->h1cdenom() << endl;

     nfd forms = nfd(hplus, verbose);
     int nforms;

     // compute the splitting operator T and the multiplicity 1
     // irreducible factors of its char poly f_T(X):
     int manual = 0;
     cerr << "Choose a splitting operator yourself (1)? or automatically (0)? ";
     cin >> manual;
     if (manual)
       {
         int ok = 0;
         while (!ok)
           {
             if (verbose)
               cout<<"Computing a splitting operator T"<<endl;
             forms.find_T();
             cout<<"Computed splitting operator and its multiplicity 1 irreducible factors\n";
             nforms = forms.nfactors;
             if (nforms==0)
               cout<<"No suitable factors, use a different operator to split!"<<endl;
           }
       }
     else
       {
         INT maxpnorm(100);
         if(verbose)
           cout << "Trying all primes P of norm up to " << maxpnorm << endl;
         Quadprime P0;
         int ok = forms.find_T_auto(maxpnorm, P0, verbose);
         assert(ok);
         nforms = forms.nfactors;
         if (verbose)
           cout << "Success with P = " << ideal_label(P0) << endl;
       }
     forms.make_irreducible_subspaces();
     assert(nforms==(int)forms.factors.size());
     cout << "Found " << nforms << " irreducible components with dimensions " << forms.dimS << endl;
     for (int j=1; j<=nforms; j++)
       {
         forms.display_basis(j);
         cout<<endl;
       }
     int ip, nap=5;
     cerr<<"Number of ap? ";  cin>>nap;

     ip = 0;
     for ( auto& P : Quadprimes::list)
       {
         if (P.divides(N))
           continue;
         ip++;
         if (ip>nap)
           break;

         vector<vec> apvec = forms.ap(P);
         cout<<"a_"<<ideal_label(P)<<" : "<<apvec<<endl;
       } // end of prime loop
  cout<<endl;
    }     // end of level loop
  exit(0);
}   // end of main()

