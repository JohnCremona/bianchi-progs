// FILE INT_HOM.CC: compute GL2 and SL2 integral homology (assuming that a geodata file exists)

// Define GP_CODE to not compute it but only output a gp script to compute it

#include "geometry.h"
#include "swan_hom.h"  // for show_invariants
#include "swan.h"  // for integral_homology

#define MAX_DISC 300
#define MIN_DISC 0

#define VERBOSE 1 // verbose setting to use if not overridden locally
#define DEBUG 1   // verbose setting to use if not overridden locally
#define PRETTY_INVARIANTS 0

#define GP_SCRIPT 1

int main ()
{
  int verbose = VERBOSE;
  int debug = DEBUG;

  long d, f, maxpnorm=100;
  vector<long> fields = valid_field_discs();
  cerr << "Enter field (0 for all up to "<<MAX_DISC<<", -1 for all): " << flush;  cin >> f;

  cerr << endl;

  if (f>0)
    fields = {(f%4==3?f:4*f)};
  for (auto D: fields)
    {
      if ((f==0) && (D>MAX_DISC))
        break;
      if (D<MIN_DISC) continue;
      d = (D%4==0? D/4: D);
      Quad::field(d,maxpnorm);

      if (D!=fields.front())
        cout << "-------------------------------------------" <<endl;

      debug = DEBUG;
      verbose = VERBOSE;

      if (verbose)
        Quad::displayfield(cout);
      else
        cout << "Field Q(sqrt("<<-d<<"))\tdiscriminant = -"<<D<<endl;

      // Read precomputed geometry data

      debug = DEBUG;
      verbose = VERBOSE;

      if (Quad::SD.read("geodata", 0)<2)
        {
          cout << "No geodata file exists for D = " << D << endl;
          cout << "Run make_geodata first!" << endl;
          continue;
        }
      Quad::SD.set_timer_status(verbose);
      Quad::setup_geometry("geodata", verbose);

      debug = DEBUG;
      verbose = VERBOSE;

#if(GP_SCRIPT)
      // Output a gp script to compute H1
      if (verbose)
        cout << "Outputting GP script to compute integral homology..." <<endl;
      Quad::SD.integral_homology_gp_code(3, debug);
      if (verbose)
        cout << "done." << endl;
#else
      // Compute integral homology

      if (verbose)
        cout << "Computing integral homology..." <<endl;

      vector<vector<INT>> invariants = Quad::SD.integral_homology(3, debug);
      cout << -D << " ";
      if (PRETTY_INVARIANTS)
        cout << "GL2 integral homology: ";
      else
        cout << "GL2 ";
      show_invariants(invariants[0], PRETTY_INVARIANTS);
      cout << endl;
      cout << -D << " ";
      if (PRETTY_INVARIANTS)
        cout << "SL2 integral homology: ";
      else
        cout << "SL2 ";
      show_invariants(invariants[1], PRETTY_INVARIANTS);
      cout << endl;

      debug = 0; //DEBUG;
      if (debug)
        {
          cout<<"Recomputing using old code..."<<flush;
          vector<vector<INT>> old_invariants = Quad::SD.old_integral_homology(3, debug);
          int ok=1;
          if (invariants[0] != old_invariants[0])
            {
              ok=0;
              cout<<"GL2: old invariants " << old_invariants[0] << " not " << invariants[0] <<endl;
            }
          if (invariants[1] != old_invariants[1])
            {
              ok=0;
              cout<<"SL2: new invariants " << old_invariants[1] << " not " << invariants[1] <<endl;
            }
          if(ok)
            cout<<"...new code agrees!"<<endl;
          else
            cout<<"...new code does not agree **********************!"<<endl;
        }
#endif
    }
}
