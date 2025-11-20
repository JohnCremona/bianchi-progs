// NEWHECKE_MODP.CC -- Hecke operators: factored characteristic
// polynomials on the new cuspidal subspace in characteristic p

#include "qidloop.h"
#include "newforms.h"

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

  // Read new poly dict cache
  int show_pols=1;
  cerr << "See the full char polys (0/1)? ";
  cin >> show_pols;
  int np;
  cerr << "How many Hecke matrices T(P)? ";
  cin >> np;
  Qideal N;

  ostringstream s;
  s << "new_poly_mod_p_dict_" << field_label();
  string poly_dict_filename = s.str();

  cout << "About to read cache of new polynomials from " << poly_dict_filename << endl;
  ifstream poly_dict_in;
  poly_dict_in.open(poly_dict_filename.c_str());
  if (poly_dict_in.is_open())
    {
      new_poly_modp_dict = input_poly_dict(poly_dict_in, ZZ(ch));
      poly_dict_in.close();
      cout << "Read cache of " << new_poly_modp_dict.size() << " new polynomials from " << poly_dict_filename << endl;
      //output_poly_dict(cout, new_poly_modp_dict);
    }
  else
    {
      cout << "No cache file " << poly_dict_filename << " exists for field " << field_label() << " yet" << endl;
    }
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
          cout << "Characteristic p = " << ch << endl;
          homspace* h = get_homspace_modp(N, ch);
  int dim = h->h1cuspdim();
  int dimnew = -1;
  cout << "Cuspidal dimension = " << dim << endl;
  ZZ_pX pol1; set(pol1); // sets to constant polynomial 1

  if (1) // dim>0)
    {
      int ntp = 0;
      for ( auto& P : Quadprimes::list)
	{
	  if (P.divides(N) || !P.is_principal())
            continue;
          ntp++;
          if (ntp>np)
            break;
          if (dimnew==0) // just to fill the cache, so higher levels
                         // do not have to recreate this hmomspace
            {
              Qideal NN=N; // copy as N is const, for alldivs()
              Quadprime PP=P; // copy as P is const, for ideal_label
              string NP = NPmodpkey(NN,PP,ch);
              new_poly_modp_dict[NP] = pol1;
              continue;
            }
          if (show_pols)
            {
              cout << "Characteristic polynomial of T(" << P << "): "<<flush;
              ZZ_pX charpol = get_full_poly_modp(N,P,ch);
              cout << str(charpol) << endl;
              cout<<"Factors:"<<endl;
              display_factors(charpol);
              cout << endl;
            }
          ZZ_pX newpol = get_new_poly_modp(N,P, ch);
          dimnew = deg(newpol);
          if (ntp==1)
            cout << "Newspace has dimension " << dimnew << endl;
          if (dimnew>0)
            {
              cout << "New characteristic polynomial of T(" << P << "): "
                   << str(newpol) << endl;
              cout<<"Factors:"<<endl;
              display_factors(newpol);
              cout << endl;
            }
        }
    }      // end of if(dim>0)

  // We output the cache after every level so we don't lose data by killing mid-level:

  // Write new poly dict cache
  cout << "About to write cache of new polynomials to " << poly_dict_filename << endl;
  ofstream poly_dict_out;
  poly_dict_out.open(poly_dict_filename.c_str());
  output_poly_dict(poly_dict_out, new_poly_modp_dict);
  poly_dict_out.close();
  cout << "Written cache of " << new_poly_modp_dict.size() << " new polynomials to " << poly_dict_filename << endl;
  }       // end of while()

  exit(0);
}       // end of main()
