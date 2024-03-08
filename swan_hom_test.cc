#include "swan_sigmas.h"
#include "swan_alphas.h"
#include "swan_tess.h"
#include "swan_hom.h"
#include "swan.h"
#include "geometry.h"

#define MAX_DISC 100

#define VERBOSE 0 // verbose setting to use if not overridden locally
#define DEBUG 0   // verbose setting to use if not overridden locally

int main ()
{
  int verbose = VERBOSE;
  int debug = DEBUG;
  int to_file=1;
  int to_screen=0;

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
      cout << "Reading previously computed sigmas and alphas from geodata/geodata_"<<d<<".at..." <<flush;
      Quad::setup_geometry(); // this sets lots of globals including alphas and sigmas, and M_alphas
      cout << "done..."<<endl;

      // Sort alphas:
      vector<vector<Quad>> pluspairs, minuspairs, fours;
      alphas = sort_alphas(alphas, pluspairs, minuspairs, fours, verbose, debug);
      sigmas = sort_singular_points(sigmas);

      // compute the M_alphas and alpha_inv from our lists:
      M_alphas = all_M_alphas(alphas, alpha_inv);

      // Find all principal polyhedra:
      verbose = VERBOSE;
      int n;
      cout << "Constructing principal polyhedra from alphas";
      if (debug) cout << ": "<<alphas;
      cout << "..." << flush;
      if (verbose) cout<<endl;
      vector<POLYHEDRON> princ_polys = principal_polyhedra(alphas, verbose);
      n = princ_polys.size();
      if (n==1)
        cout << "done: 1 principal polyhedron constructed:"<<endl;
      else
        cout << "done: " << n << " principal polyhedra constructed:"<<endl;
      map<string,int> poly_counts;
      for (const auto& P: princ_polys)
        poly_counts[poly_name(P)]++;
      for (const auto& pc : poly_counts)
        cout<<pc.second<<" "<<pc.first << (pc.second>1?"s":"") << endl;
      // Find all singular polyhedra:
      verbose = VERBOSE;
      cout << "Constructing singular polyhedra..."<<flush;
      vector<POLYHEDRON> sing_polys = singular_polyhedra(sigmas, alphas, verbose);
      n = sing_polys.size();
      if (n==1)
        cout << "done: 1 singular polyhedron constructed" <<endl;
      else
        cout << "done: " << n << " singular polyhedra constructed" <<endl;
      map<string,int> spoly_counts;
      for (const auto& P: sing_polys)
        spoly_counts[poly_name(P)]++;
      for (const auto& pc : spoly_counts)
        cout<<pc.second<<" "<<pc.first << (pc.second>1?"s":"") << endl;

      vector<POLYHEDRON> all_polys = sing_polys;
      all_polys.insert(all_polys.end(), princ_polys.begin(), princ_polys.end());

      cout << "\nFinding all faces up to GL2-equivalence" << endl;
      vector<vector<int>> M32;
      vector<int> redundant_faces;
      verbose = VERBOSE;
      auto all_faces = get_faces(all_polys, alphas, sigmas, M32, redundant_faces, verbose);
      //int nfaces = all_faces.size();

      verbose = VERBOSE;

      // split up faces into 4 types for reporting and output:
      vector<CuspList> aaa_triangles, aas_triangles, squares, hexagons;
      cout << "Faces up to GL2-action and reflection:\n";
      int sing, red, i=0;
      for (const auto& face: all_faces)
        {
          sing = is_face_singular(face, sigmas);
          red = std::find(redundant_faces.begin(), redundant_faces.end(), i)!=redundant_faces.end();
          int n = face.size();
          string s;
          if (n==4)
            {
              if (!red) squares.push_back(face);
              s = "square";
              if (red) s+= " (redundant)";
            }
          else
            {
              if (n==6)
                {
                  if (!red) hexagons.push_back(face);
                  s = "hexagon";
                  if (red) s+= " (redundant)";
                }
              else
                {
                  if (sing)
                    {
                      if (!red) aas_triangles.push_back(face);
                      s = "aas triangle";
                      if (red) s+= " (redundant)";
                    }
                  else
                    {
                      if (!red) aaa_triangles.push_back(face);
                      s = "aaa triangle";
                      if (red) s+= " (redundant)";
                    }
                }
            }
          if (verbose)
            cout<<i<<" ("<<s<<"): "<<face<<endl;
          i++;
        }

      verbose = VERBOSE;
      int all_ok = 1;
      cout<<aaa_triangles.size()<<" aaa-triangles\n";
      for ( const auto& face : aaa_triangles)
        {
          if (verbose) cout <<face << " --> ";
          POLYGON P = make_polygon(face, alphas, sigmas, sing);
          if (verbose) cout <<face << " -->  ["<<P.first<<","<<P.second<<"]"<<endl;
          int ok = check_aaa_triangle(P, verbose);
          if (!ok)
            cout<<"aaa-triangle "<<face<<" --> ["<<P.first<<","<<P.second<<"] fails"<<endl;
          all_ok = ok &&all_ok;
        }
      cout<<aas_triangles.size()<<" aas-triangles\n";
      for ( const auto& face : aas_triangles)
        {
          POLYGON P = make_polygon(face, alphas, sigmas, sing);
          int ok = check_aas_triangle(P, verbose);
          if (!ok)
            cout<<"aas-triangle "<<face<<" --> ["<<P.first<<","<<P.second<<"] fails"<<endl;
          all_ok = ok &&all_ok;
        }
      cout<<squares.size()<<" squares\n";
      for ( const auto& face : squares)
        {
          POLYGON P = make_polygon(face, alphas, sigmas, sing);
          int ok = check_square(P);
          if (!ok)
            cout<<"square "<<face<<" --> ["<<P.first<<","<<P.second<<"] fails"<<endl;
          all_ok = ok &&all_ok;
        }
      cout<<hexagons.size()<<" hexagons\n";
      for ( const auto& face : hexagons)
        {
          POLYGON P = make_polygon(face, alphas, sigmas, sing);
          int ok = check_hexagon(P);
          if (!ok)
            cout<<"hexagon "<<face<<" --> ["<<P.first<<","<<P.second<<"] fails"<<endl;
          all_ok = ok &&all_ok;
        }
      if (all_ok)
        cout<<"all encodings check out OK" << endl;
      else
        {
          cout<<"*****************not all encodings check out OK" << endl;
          exit(1);
        }
      if (to_file||to_screen) cout << "geodata encodings of faces";
      if (to_file) cout << " output to geodata file";
      if (to_file||to_screen) cout << "\n";
      output_faces({aaa_triangles, squares, hexagons, aas_triangles},
                   alphas, sigmas, to_file, to_screen);

      // Compute integral homology

      debug = DEBUG;
      if (debug)
        {
          cout<<"alphas: "<<alphas<<endl;
          cout<<"sigmas: "<<sigmas<<endl;
          cout<<"pluspairs: "<<pluspairs<<endl;
          cout<<"minuspairs: "<<minuspairs<<endl;
          cout<<"fours: "<<fours<<endl;
          cout<<"faces: "<<all_faces<<endl;
        }
      vector<vector<int>> invariants = integral_homology(all_faces,
                                                         alphas, sigmas,
                                                         pluspairs, minuspairs, fours,
                                                         3, debug);
      cout << "GL2 integral homology invariants: " << invariants[0] << endl;
      cout << "SL2 integral homology invariants: " << invariants[1] << endl;


      cout<<"----------------------------------------------------------------------------------\n";
    }
}
