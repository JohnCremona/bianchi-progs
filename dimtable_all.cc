// dimtable_all
//
// For all implemented fields output a file with the following data for level 1
//
// D SL2-all SL2-cusp GL2-all GL2-cusp NGL2-all NGL2-cusp
//
// where D = field discriminant, G-all and G-cusp are dimensions of
// the homology and cuspidal homology for G = SL(2,O_K), G =
// GL(2,O_K), and G = NGL(2,O_K), the latter being the normaliser of
// GL(2,O_K) in GL(2,K), hence only relevant for even class number
// (otherwise GL=NGL).

#include "qidloop.h"
#include "homspace.h"

int main ()
{
  vector<long> fields = valid_field_discs(0); // all
  long nfields = fields.size();
  long min_disc, max_disc = *std::max_element(fields.begin(), fields.end());
  cerr << nfields << " fields available, discriminants <= " << max_disc << endl;

  cerr << "Enter min and max abs disc (0 0 for all): " << flush;  cin >> min_disc >> max_disc;
  fields = valid_field_discs(max_disc);
  long max(2000);
  cout << "Table of dimensions of ";
  cout<<"level 1 homology over Q(sqrt(-d))" << endl;
  cout<<"\n";
  cout<<"   (For fields with even class number, The GL2 dimensions include all unramified quadratic character subspaces, while"<<endl;
  cout<<"   the NGL2 dimensions are of just the trivial character subspace.)"<<endl;

  cout << "\nDiscriminant\tClass Number\t\t";
  cout << "SL2" << "\t\t" << "GL2" << "\t\tNGL2" << "\n\t\t\t\t";
  cout<<"all cuspidal\tall cuspidal\tall cuspidal" << endl;

  for (auto D : fields)
    {
      if (D<min_disc)
        continue;
      long d = (D%4==0? D/4: D);
      Quad::field(d,max);
      int n2r = Quad::class_group_2_rank>0;
      Qideal N(ONE);
      cout << Quad::disc << "\t\t";        // field discriminant
      cout << Quad::class_number << "\t\t"; // class number

      // SL2

      homspace h1(N, 0);  //level, plusflag=0 for SL2
      long dimall = h1.h1dim();
      long dimcusp = h1.h1cuspdim();
      cout << dimall << "   " << dimcusp << "\t\t";

      // GL2

     homspace h1plus(N, 1);  //level, plusflag=1 for GL2
     dimall = h1plus.h1dim();
     dimcusp = h1plus.h1cuspdim();
     cout << dimall << "   " << dimcusp;

     if (n2r)
       {
         dimall = h1plus.trivial_character_subspace_dimension(0);
         dimcusp = h1plus.trivial_character_subspace_dimension(1);
         cout << "\t\t" << dimall << "   " << dimcusp;
       }
     cout << endl;
    }
}
