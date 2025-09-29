// NEWHECKE_MODP.CC -- Hecke operators: factored characteristic
// polynomials on the new cuspidal subspace in characteristic p

#include "eclib.h"

#include "matprocs.h"
#include "qidloop.h"
#include "newforms.h"
#include "nfd.h"

//#define LOOPER

#define MAXPRIME 10000

int main(void)
{
  long d, maxpnorm(MAXPRIME);
  cerr << "Enter field: " << flush;
  cin >> d;
  if (!check_field(d))
   {
     cerr<<"field must be one of: "<<valid_fields<<endl;
     exit(1);
   }
  scalar ch(0);
  cerr << "Enter characteristic p (prime): " << flush;  cin >> ch;
  ZZ_p::init(ZZ(ch));
  Quad::field(d,maxpnorm);
  Quad::displayfield(cout);

  int show_pols=1;
  cerr << "See the full char polys (0/1)? ";
  cin >> show_pols;
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
          homspace* h = get_homspace_modp(N, ch);
  int dim = h->h1cuspdim();
  cout << "Cuspidal dimension = " << dim << endl;

  if (dim>0)
    {
      int ntp = 0;
      for ( auto& P : Quadprimes::list)
	{
	  if (P.divides(N) || !P.is_principal())
            continue;
          ntp++;
          if (ntp>np) break;
          if (show_pols)
            {
              cout << "Characteristic polynomial of T(" << P << ")"<<endl;
              ZZ_pX charpol = get_full_poly_modp(N,P,ch);
              cout << "Coefficients: " << charpol << endl;
              cout<<"Factors:"<<endl;
              display_factors(charpol);
              cout << endl;
            }
          ZZ_pX newpol = get_new_poly_modp(N,P, ch);
          int dimnew = deg(newpol);
          if (dimnew==0)
            {
              cout << "Newspace is trivial" << endl;
              break;
            }
          if (ntp==1)
            {
              cout << "Newspace has dimension " << dimnew << endl;
            }
          cout << "New characteristic polynomial of T(" << P << "):" << endl;
          cout << "Coefficients: " << newpol << endl;
          cout<<"Factors:"<<endl;
          display_factors(newpol);
          cout << endl;
        }
    }      // end of if(dim>0)
}       // end of while()
exit(0);
}       // end of main()
