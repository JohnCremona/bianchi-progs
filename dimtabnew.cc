// DIMTABNEW.CC  -- Table of dimensions of Bianchi cusp forms at level N with trivial character

// Columns:  Field Weight(2) Level dimcusp dimcuspnew

// where dimcusp = dimension of cusp forms (with trivial character)
//       dimcuspnew = dimension of new cusp forms (with trivial character)

// Inputs (prompted for): d (field)
//                        max_norm (upper bound on level norm)

// This program requires precomputed data for all levels of norm up to
// max_norm.  The new (cuspidal, trivial character) homology
// dimnesions are read from data files, converted to new Bianchi
// dimensions and stored in a dict (with level labels as keys); the
// old dimensions are computed from the new dimensions at lower
// levels.

#include "qidloop.h"
#include "newforms.h"
//#define MODP
scalar modulus = default_modulus<scalar>();

int main ()
{
  long d, maxpnorm(2000);
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

  long min_norm=1, max_norm; Quad n;
  int both_conj=1;
  //cerr<<"Both conjugates? (0/1) "; cin >> both_conj;

  cerr<<"Enter max norm for Quad loop: ";
  cin >> max_norm;
  cerr<<endl;

  Quad::field(d,maxpnorm);

  // Each new eigensystem stored which is stored gives rise to either
  // 2^r Bianchi newforms (if not self-twist) or 2^{r-1} Bianchi
  // newforms (if self-twist by a quadratic genus character), which
  // are all twists of eachother by the group of unramified quadratic
  // twists.
  int nchi = 1<<(Quad::class_group_2_rank);
  int nchi2 = nchi/2; // only used when nchi is even

  cout << "# Table of dimensions of ";
  if (ch!=0) cout<<"mod "<<ch<<" ";
  cout<<"weight 2 Bianchi cuspidal forms and newforms for GL2 over Q(sqrt(-"<<d<<"))" << endl;
  if (Quad::class_group_2_rank>0)
    cout<<"# (with trivial character)"<<endl;
  cout << "# Field\tWeight\tLevel\t";
  cout << "dim(cuspidal)\tdim(cuspidal, new)" << endl;

  // dictionary with key = level label, value = new dimension
  map<string, int> dims;

  Qidealooper loop(min_norm, max_norm, both_conj, 1); // sorted within norm
  while( loop.not_finished() )
   {
     Qideal N = loop.next();
     string level_label = ideal_label(N);
     vector<Qideal> Ndivisors = alldivs(N);
     // cout<<"Genus characters of divisors of "<<level_label<<":"<<endl;
     // cout << "\t\t"<<Quad::all_disc_factors <<endl;
     // for ( const auto& I : Ndivisors)
     //   {
     //     cout<<ideal_label(I)<<" --> \t";
     //     for (const auto& d : Quad::all_disc_factors)
     //       cout<<I.genus_character(d)<<" ";
     //     cout<<endl;
     //   }
     cout << "\t"<< d << "\t2\t";                  // field and weight
     cout << level_label<<"\t\t"; // level

     newforms nfdata(N, modulus, 0, ch); // verbose=0
     if (nfdata.read_from_file())
       {
         int dimcuspnew, dimcuspold, dimcusp;
         if (nchi == 1)
           {
             dimcuspnew = nfdata.n1ds + nfdata.n2ds;
           }
         else
           {
             // multiplicity nchi for non-self-twist eigensystems
             dimcuspnew = nchi*(nfdata.new1dims[0] + nfdata.new2dims[0]);
             // multiplicity nchi/2 for self-twist eigensystems
             for (int i=1; i<nchi; i++)
               dimcuspnew += nchi2*(nfdata.new1dims[i] + nfdata.new2dims[i]);
           }
         dims[level_label] = dimcuspnew;
         dimcuspold = 0;
         for ( auto D: Ndivisors)
           {
             if (N==D)
               break;
             Qideal M = N/D;
             dimcuspold += dims[ideal_label(D)] * ndivs(M);
           }
         dimcusp = dimcuspold + dimcuspnew;
         cout << dimcusp << "\t\t";
         cout << dimcuspnew << endl;
       }
     else
       {
         cout<<"no newforms file for level " << level_label<<"! "<<endl;
         break;
       }
   }
}
