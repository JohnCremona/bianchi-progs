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
     homspace hplus(N, modulus, 1);
     // int dimH = hplus.h1cuspdim();
     // cout << "dimension = " << dimH << endl;
     // cout << "denom     = " << hplus.h1cdenom() << endl;

     nfd forms = nfd(&hplus, verbose);

     // compute the splitting operator T and the multiplicity 1
     // irreducible factors of its char poly f_T(X):
     int ok = 0;
     while (!ok)
       {
         cout<<"Computing a splitting operator T"<<endl;
         forms.make_T();
         cout<<"Computed splitting operator and its multiplicity 1 irreducible factors\n";
         ok = forms.factors.size();
         if (!ok)
           cout<<"No suitable factors, use a different operator to split!"<<endl;
       }
     forms.make_irreducible_subspaces();
     int nforms = forms.nfactors;
     cout<<"Finished constructing "<<nforms
         <<" irreducible subspaces of dimensions "<<forms.dimS<<endl<<endl;
     for (int j=1; j<=nforms; j++)
       {
         forms.display_basis(j);
         cout<<endl;
       }
     int ip, nap=5;
     cout<<"Number of ap? ";  cin>>nap; cout<<endl;

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
    }     // end of level loop
  cout<<endl;
  exit(0);
}   // end of main()

