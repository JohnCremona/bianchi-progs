// FILE RNFD.CC:  Read and display precomputed Newspaces (d-dimensional newforms)
/////////////////////////////////////////////////////////////////////////////////

//#define LOOPER
#ifdef LOOPER
#include "qidloop.h"
#endif
#include "nfd.h"
#include "field.h"
#include <eclib/pari_init.h>

#define MAXPRIME 10000

int main()
{
  cout << "Program rnfd: read and display precomputed Bianchi newforms of arbitrary dimension." << endl;
  eclib_pari_init();

  scalar modulus = default_modulus<scalar>();
#if (SCALAR_OPTION==3)
  //  NextPrime(modulus, power_ZZ(to_ZZ(2),256));
  NextPrime(modulus, pow(to_ZZ(2),512));
#endif
  long d, maxpnorm(MAXPRIME);
  cerr << "Enter field: " << flush;
  cin >> d;
  if (d<1)
    {
      cout << "Field parameter d must be positive for Q(sqrt(-d))" << endl;
      exit(0);
    }
  Quad::field(d,maxpnorm);
  Quad::displayfield(cout);

  Quad n;
  Qideal N;
#ifdef LOOPER
  long firstn, lastn;
  cerr<<"Enter first and last norm for levels: ";
  cin >> firstn >> lastn;

  Qidealooper loop(firstn, lastn, 1, 1); // 1 = both conjugates; 1 = sorted within norm
  while( loop.not_finished() )
    {
      N = loop.next();
#else
  while(cerr<<"Enter level (ideal label or generator): ", cin>>N, !N.is_zero())
    {
#endif
      cout << endl;

      string Nlabel = ideal_label(N);
      cout << ">>>> Level " << Nlabel <<" = "<<gens_string(N)<<", norm = "<<N.norm()<<" <<<<" << endl;

      Newspace NS;
      NS.input_from_file(N);
      int nnf = NS.nforms();
      cout << nnf << (nnf==1? " newform" : " newforms") << endl;
      NS.display_newforms(1, 1); // aP: yes;  AL: yes (principal eigs, triv_char_only: no)
    }     // end of level loop
  cout << endl;
  exit(0);
}   // end of main()

