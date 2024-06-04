#include "geometry.h"
#include "swan_tess.h" // for remake_triangle
#include "swan_hom.h"  // for integral_homology
#include "swan.h"  // for integral_homology

#define MAX_DISC 100
#define MIN_DISC 0

#define VERBOSE 0 // verbose setting to use if not overridden locally
#define DEBUG 0   // verbose setting to use if not overridden locally

int main ()
{
  int verbose = VERBOSE;
  int debug = DEBUG;

  long d, f, maxpnorm=100;
  vector<long> fields = valid_field_discs();
  cerr << "Enter field (0 for all up to "<<MAX_DISC<<", -1 for all): " << flush;  cin >> f;

  cerr << endl;

  if (f>0)
    fields = {f};
  for (auto D: fields)
    {
      if ((f==0) && (D>MAX_DISC))
        break;
      if (D<MIN_DISC) continue;
      d = (D%4==0? D/4: D);
      Quad::field(d,maxpnorm);

      if (D!=fields.front())
        cout << "-------------------------------------" <<endl;

      debug = DEBUG;
      verbose = VERBOSE;

      if (verbose)
        Quad::displayfield(cout);
      else
        cout << "Field Q(sqrt("<<-d<<"))\tdiscriminant = "<<D<<endl;

      // Read precomputed geometry data

      if (verbose)
        cout << "Reading previously computed sigmas and alphas from geodata/geodata_"<<d<<".dat..." <<flush;
      Quad::setup_geometry(verbose);
      if (verbose)
        cout << "done..."<<endl;

      // Compute integral homology

      debug = DEBUG;
      verbose = VERBOSE;
      if (verbose)
        cout << "Computing integral homology..." <<endl;

      vector<vector<int>> invariants = Quad::SD.integral_homology(3, debug);

      cout << "GL2 integral homology: "; show_invariants(invariants[0]); cout << endl;
      cout << "SL2 integral homology: "; show_invariants(invariants[1]); cout << endl;

      cout<<"----------------------------------------------------------------------------------\n";
    }
}
