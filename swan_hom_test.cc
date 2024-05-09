#include "geometry.h"
#include "swan_tess.h" // for remake_triangle
#include "swan_hom.h"  // for integral_homology

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

      if (verbose)
        {
          Quad::displayfield(cout);
        }
      else
        cout << "Field Q(sqrt("<<-d<<"))\tdiscriminant = "<<D<<endl;

      // Read precomputed alphas and sigmas
      if (verbose)
        cout << "Reading previously computed sigmas and alphas from geodata/geodata_"<<d<<".dat..." <<flush;
      Quad::setup_geometry(); // this sets lots of globals including alphas and sigmas, and M_alphas
      if (verbose)
        cout << "done..."<<endl;

      verbose = VERBOSE;
      if (verbose)
        cout << "Reading encoded faces from geodata/geodata_"<<d<<".dat..." <<flush;
      vector<vector<POLYGON>> all_polys = read_polygons(verbose);
      if (verbose)
        cout << "done..."<<endl;

      verbose = VERBOSE;
      if (verbose)
        cout << "Converting encoded faces to vertex lists..." <<flush;

      // convert POLYGON encodings to actual faces, by type:
      int nT = all_polys[0].size(), nU = all_polys[1].size(), nQ = all_polys[2].size(), nH = all_polys[3].size();
      vector<CuspList> all_faces(nT+nU+nQ+nH);
      std::transform(all_polys[0].begin(), all_polys[0].end(), all_faces.begin(),
                     [] (const POLYGON& P) {return remake_triangle(P, alphas, sigmas, 0);});
      std::transform(all_polys[1].begin(), all_polys[1].end(), all_faces.begin()+nT,
                     [] (const POLYGON& P) {return remake_triangle(P, alphas, sigmas, 1);});
      std::transform(all_polys[2].begin(), all_polys[2].end(), all_faces.begin()+nT+nU,
                     [] (const POLYGON& P) {return remake_quadrilateral(P, alphas);});
      std::transform(all_polys[3].begin(), all_polys[3].end(), all_faces.begin()+nT+nU+nQ,
                     [] (const POLYGON& P) {return remake_hexagon(P, alphas);});
      if (verbose)
        cout << "done..."<<endl;

      // Compute integral homology

      debug = DEBUG;
      verbose = VERBOSE;
      if (verbose)
        cout << "Computing integral homology..." <<endl;

      vector<vector<int>> invariants = integral_homology(all_faces, 3, debug);

      cout << "GL2 integral homology: "; show_invariants(invariants[0]); cout << endl;
      cout << "SL2 integral homology: "; show_invariants(invariants[1]); cout << endl;


      cout<<"----------------------------------------------------------------------------------\n";
    }
}
