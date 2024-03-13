#include "swan_sigmas.h"
#include "swan_alphas.h"
#include "swan_tess.h"
#include "swan_hom.h"
#include "swan.h"
#include "geometry.h"

#define MAX_DISC 100
#define MIN_DISC 910

#define VERBOSE 0 // verbose setting to use if not overridden locally
#define DEBUG 0   // verbose setting to use if not overridden locally

int main ()
{
  int verbose = VERBOSE;
  int debug = DEBUG;

  long d, f, max=100;
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
      Quad::field(d,max);

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
      vector<CuspList> all_faces;

      // convert POLYGON encodings to actual faces, by type:
      for (const auto& P : all_polys[0]) // T triangles
        {
          all_faces.push_back(remake_triangle(P, alphas, sigmas, 0));
        }
      for (const auto& P : all_polys[1]) // U triangles
        {
          all_faces.push_back(remake_triangle(P, alphas, sigmas, 1));
        }
      for (const auto& P : all_polys[2]) // Q squares
        {
          all_faces.push_back(remake_quadrilateral(P, alphas));
        }
      for (const auto& P : all_polys[3]) // H hexagons
        {
          all_faces.push_back(remake_hexagon(P, alphas));
        }
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
