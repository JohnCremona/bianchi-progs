// FILE TNFD.CC:  test program for Newforms (d-dimensional newform) class
//////////////////////////////////////////////////////////////////////////
//
// Adapted from the similar program (over Q) in eclib
//
//////////////////////////////////////////////////////////////////////////

//#define LOOPER
#ifdef LOOPER
#include "qidloop.h"
#endif
#include "nfd.h"

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
  cerr << "Enter field: " << flush;  cin >> d;
  Quad::field(d,maxpnorm);
  Quad::displayfield(cout);
  //int n2r = Quad::class_group_2_rank;

  int verbose=1;
  cerr << "Verbose output? (0/1) "; cin >> verbose;
  int nap=5;
  cerr<<"Number of ap? ";
  cin>>nap;

  Quad n;
  Qideal N;
#ifdef LOOPER
  long firstn, lastn;
  cerr<<"Enter first and last norm for Quad loop: ";
  cin >> firstn >> lastn;
  Qidealooper loop(firstn, lastn, 0, 1); // 0 = not both conjugates; 1 = sorted within norm
  while( loop.not_finished() )
    {
      N = loop.next();
#else
  while(cerr<<"Enter level (ideal label or generator): ", cin>>N, !N.is_zero())
    {
#endif
      cout << endl;
      cout << ">>>> Level " << ideal_label(N) <<" = "<<gens_string(N)<<", norm = "<<N.norm()<<" <<<<" << endl;
     homspace* hplus = get_homspace(N, modulus);

     INT maxpnorm(100);
     Newforms forms(hplus, maxpnorm, verbose);
     if (!forms.ok())
       {
         cout << "Failed to find a splitting operator with maxpnorm = " << maxpnorm << endl;
         continue; // to next level
       }
     cout << "Success with " << forms.splitopname() << endl;
     int nnf = forms.nforms();
     cout << "Found " << nnf << " newforms";
     if (nnf)
       cout << " with dimensions " << forms.dimensions();
     cout << endl;
     if (!nnf)
       continue;
     forms.display_bases();
     int ip = 0;
     for ( auto& P : Quadprimes::list)
       {
         if (P.divides(N))
           continue;
         ip++;
         if (ip>nap)
           break;
         matop T = AutoHeckeOp(P,N);
         vector<vec> apvec = forms.eig(T);
         cout<<T.name() << ":\t" <<apvec<<endl;
       } // end of prime loop
     //  cout<<endl;
    }     // end of level loop
  exit(0);
}   // end of main()

