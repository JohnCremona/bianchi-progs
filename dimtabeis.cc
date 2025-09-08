// DIMTABEIS.CC  -- Table of dimensions of homology at level N with trivial character

// Columns:  Field Weight(2) Level dimall dimcusp dimeis

// where dimall  = dimension of homology subspace with trivial character
//       dimcusp = dimension of cuspidal homology subspace with trivial character
//       dimeis  = dimension of non-cuspidal homology subspace with trivial character

// Inputs (prompted for): d (field)
//                        both_conj (flag to include both conjugates in loop over levels)
//                        min_norm (lower bound on level norm in loop over levels)
//                        max_norm (upper bound on level norm in loop over levels)

#include "qidloop.h"
#include "homspace.h"
//#define MODP

scalar modulus = default_modulus<scalar>();

int main ()
{
  long d, maxpnorm(1000);
  cerr << "Enter field: " << flush;  cin >> d;
  if (!check_field(d))
    {
      cerr<<"field must be one of: "<<valid_fields<<endl;
      exit(1);
    }
  scalar ch(0);
#ifdef MODP
  cerr << "Enter characteristic (0 or prime): " << flush;  cin >> ch;
#endif
  long min_norm, max_norm; Quad n;
  int both_conj;
  cerr<<"Both conjugates? (0/1) "; cin >> both_conj;
  int verbose=0;

  cerr<<"Enter first and last norm for Quad loop: ";
  cin >> min_norm >> max_norm;
  cerr<<endl;
  Quad::field(d,maxpnorm);


 cout << "# Table of dimensions of ";
 if (ch!=0) cout<<"mod "<<ch<<" ";
 cout<<"weight 2 Bianchi cuspidal and Eisenstein forms for GL2 over Q(sqrt(-"<<d<<"))" << endl;
 if (Quad::class_group_2_rank>0)
   cout<<"# (with trivial character)"<<endl;
 cout << "# Field\tWeight\tLevel\t";
 cout << "dim(all)\tdim(cuspidal)\tdim(eisenstein)" << endl;

 Qidealooper loop(min_norm, max_norm, both_conj, 1); // sorted within norm
 while( loop.not_finished() )
   {
     Qideal N = loop.next();
     cout << "\t"<< d << "\t2\t";                  // field and weight
     cout << ideal_label(N)<<"\t"; // level
     homspace hplus(N, modulus, 1, verbose, ch);  // plusflag=1
     pair<int,int> dims = hplus.trivial_character_subspace_dimensions();
     int dimcusp = dims.second;
     int dimall = dims.first;
     int dimeis = dimall-dimcusp;
     cout << dimall << "\t\t" << dimcusp << "\t\t" << dimeis << endl;
   }

}
