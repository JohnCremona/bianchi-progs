#include "swan_sigmas.h"
#include "swan_alphas.h"
#include "swan_tess.h"
#include "swan_hom.h"
#include "swan.h"
#include "geometry.h"

#define MAX_DISC 100

#define VERBOSE 0 // verbose setting to use if not overridden locally
#define DEBUG 0   // verbose setting to use if not overridden locally

vector<RatQuad> test_singular_points(int output_level=0)
{
  if (output_level>=3)
    for (auto I : Quad::class_group)
      {
        cout<<"Ideal class ["<<I<<"]: ";
        cout<<"singular points "<<singular_points_in_class(I,(output_level>3))<<endl;
      }
  auto sigmas = singular_points();
  if (output_level>=1)
    cout << "Number of singular points, including oo: "<<sigmas.size()<<endl;
  if (output_level>=2)
    cout << "Unsorted singular points: "<<sigmas<<endl;
  sigmas = sort_singular_points(sigmas);
  if (output_level>=1)
    cout << "Sorted singular points: "<<sigmas<<endl;
  int to_file=0; //(output_level>=1);
  int to_screen=0; //(output_level>=2);
  output_singular_points(sigmas, to_file, to_screen);
  return sigmas;
}

void test_principal_cusps(int n1=20, int n2=100)
{
  for (int n=1; n<=n1; n++)
    {
      auto alist = principal_cusps_of_dnorm(n);
      cout << "Principal cusps of denominator norm "<<n<<": "<<alist<<endl;
    }
  auto alist = principal_cusps_of_dnorm_up_to(n2);
  cout << "Principal cusps of denominator norm up to "<<n2<<": "<<alist<<endl;
}

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
    {
      // convert to an absolute discriminant if necessary
      if (f%4 !=3) f*=4;
      fields = {f};
    }
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

      //test_principal_cusps(20,20);

      cout << "Finding sigmas and alphas..."<<flush;

      verbose = VERBOSE;
      debug = DEBUG;
      auto new_sigmas = singular_points();
      auto sorted_sigmas = sort_singular_points(new_sigmas);
      if (debug)
        cout<<"Before sorting, "<<new_sigmas.size()<<" sigmas:\n"<<new_sigmas<<endl;
      if (debug||verbose)
        cout<<"After  sorting, "<<sorted_sigmas.size()<<" sigmas:\n"<<sorted_sigmas<<endl;
      new_sigmas = sorted_sigmas;

      //test_principal_cusps(20, 30);
      verbose = VERBOSE;
      debug = DEBUG;
      auto new_alphas = covering_alphas(new_sigmas, verbose);
      INT maxn = max_dnorm(new_alphas);
      if (verbose)
        {
          cout << new_alphas.size() << " covering alphas";
          if (debug) cout << ": " << new_alphas || "\n";
          cout << " with max dnorm = " << maxn <<endl;
        }

      verbose = VERBOSE;
      debug = DEBUG;
      new_alphas = saturate_covering_alphas(new_alphas, new_sigmas, maxn, debug, verbose);
      maxn = max_dnorm(new_alphas);
      if (verbose)
        cout << new_alphas.size() << " saturated alphas, max denom norm = " << maxn <<endl;
      if (debug || (verbose && new_alphas.size()<20))
        cout << new_alphas << endl;
      auto points = triple_intersections(new_alphas);
      RAT minht = points[0].t2;
      std::for_each(points.begin(), points.end(),
                    [&minht](H3point P) {minht = min(minht, P.t2);});
      if (verbose)
        cout << points.size() << " vertices, min square height = " << minht <<endl;

      // Check that all alphas are in the rectangle:
      assert (std::all_of(new_alphas.begin(), new_alphas.end(), [](const RatQuad& a) {return a.in_rectangle();}));
      // NB sigmas not in rectangle in the following cases (not changed, for consistency with old code):
      // d=7(8), we use (1-w)/2 though (w-1)/2 is in rectangle
      // d=3(12), we use (-1-w)/3 though (2-w)/3 is in rectangle
      // assert (std::all_of(new_sigmas.begin(), new_sigmas.end(), [](const RatQuad& a) {return a.in_rectangle();}));

      // Sort and output alphas:
      vector<vector<Quad>> pluspairs, minuspairs, fours;
      verbose = VERBOSE;
      debug = DEBUG;
      auto sorted_alphas = sort_alphas(new_alphas, pluspairs, minuspairs, fours, verbose, debug);
      if (debug)
        cout<<"Before sorting, alphas:\n"<<new_alphas<<endl;
      if (debug||verbose)
        cout<<"After  sorting, alphas:\n"<<sorted_alphas<<endl;
      new_alphas = sorted_alphas;
      cout<<"...done, now outputting"<<endl;
      output_alphas(pluspairs, minuspairs, fours, to_file, to_screen);
      output_singular_points(new_sigmas, to_file, to_screen);


      if (D<1000)
        {
      // Compare with precomputed alphas and sigmas
      cout << "Testing newly computed sigmas and alphas with old..." <<flush;
      Quad::setup_geometry(); // this sets lots of globals including alphas and sigmas, and M_alphas
      cout << "read in old data..."<<flush;
      Quad t;
      if (compare_CuspLists_as_sets(sigmas, new_sigmas))
        cout << "sigmas agree..." <<flush;
      else
        {
          if (compare_CuspLists_as_sets_mod_translation(sigmas, new_sigmas))
            cout << "sigmas agree (up to translation)..." <<flush;
          else
            {
              cout << "sigmas DO NOT agree:\n";
              cout << sigmas.size() << " old sigmas";
              // if (verbose)
                cout << ": " <<sigmas;
              cout << endl;
              cout << new_sigmas.size() << " new sigmas";
              // if (verbose)
                cout << ": " <<new_sigmas;
              cout << endl;
              cout << "old not in new:\n";
              for (const auto& a : sigmas)
                if (cusp_index_with_translation(a,new_sigmas,t)==-1)
                  cout << a << " " << a.coords(1) << endl;
              cout << "new not in old:\n";
              for (const auto& a : new_sigmas)
                if (cusp_index_with_translation(a,sigmas,t)==-1)
                  cout << a << " " << a.coords(1) << endl;
              exit(1);
            }
        }
      if (compare_CuspLists_as_sets(alphas, new_alphas))
        cout << "alphas agree!" <<endl;
      else
        {
          if (compare_CuspLists_as_sets_mod_translation(alphas, new_alphas))
            cout << "alphas agree (up to translation)!" <<endl;
          else
            {
              cout << "alphas DO NOT agree:\n";
              cout << alphas.size() << " old alphas";
              if (verbose) cout << ": " <<alphas;
              cout << endl;
              cout << new_alphas.size() << " new alphas";
              if (verbose) cout << ": " <<new_alphas;
              cout << endl;
              cout << "old not in new:\n";
              for (const auto& a : alphas)
                if (cusp_index_with_translation(a,new_alphas,t)==-1)
                  cout << a << " " << a.coords(1) << endl;
              cout << "new not in old:\n";
              for (const auto& a : new_alphas)
                if (cusp_index_with_translation(a,alphas,t)==-1)
                  cout << a << " " << a.coords(1) << endl;
              exit(1);
            }
        }
        } // end of code to compare new and old

      alphas = new_alphas; // overwrite the globals with our lists
      sigmas = new_sigmas; //

      // compute the M_alphas and alpha_inv from our lists:
      M_alphas = all_M_alphas(alphas, alpha_inv);

# if(0)
      // Look at neighbours of each finite singular point
      cout << "Neighbours of each finite singular point:" << endl;
      for (const auto& sigma: sigmas)
        {
          if (sigma.is_infinity()) continue;
          auto nbrs = sorted_neighbours(sigma, alphas);
          auto n = nbrs.size();
          cout << sigma << " has " << n << " neighbours (sorted): " << nbrs << endl;

          // cout << sigma << "\t has coords\t " << sigma.coords(1) <<endl;
          // for (auto a : nbrs)
          //   {
          //     RatQuad b = a-sigma;
          //     cout <<a<<" - "<<sigma<<" = "<<b<< "\t has coords\t " << b.coords(1) <<endl;
          //   }
        }
#endif
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
          if (verbose) cout <<face << " -->  ["<<P.indices<<","<<P.shifts<<"]"<<endl;
          int ok = check_aaa_triangle(P, verbose);
          if (!ok)
            cout<<"aaa-triangle "<<face<<" --> ["<<P.indices<<","<<P.shifts<<"] fails"<<endl;
          all_ok = ok &&all_ok;
        }
      cout<<aas_triangles.size()<<" aas-triangles\n";
      for ( const auto& face : aas_triangles)
        {
          POLYGON P = make_polygon(face, alphas, sigmas, sing);
          int ok = check_aas_triangle(P, verbose);
          if (!ok)
            cout<<"aas-triangle "<<face<<" --> ["<<P.indices<<","<<P.shifts<<"] fails"<<endl;
          all_ok = ok &&all_ok;
        }
      cout<<squares.size()<<" squares\n";
      for ( const auto& face : squares)
        {
          POLYGON P = make_polygon(face, alphas, sigmas, sing);
          int ok = check_square(P);
          if (!ok)
            cout<<"square "<<face<<" --> ["<<P.indices<<","<<P.shifts<<"] fails"<<endl;
          all_ok = ok &&all_ok;
        }
      cout<<hexagons.size()<<" hexagons\n";
      for ( const auto& face : hexagons)
        {
          POLYGON P = make_polygon(face, alphas, sigmas, sing);
          int ok = check_hexagon(P);
          if (!ok)
            cout<<"hexagon "<<face<<" --> ["<<P.indices<<","<<P.shifts<<"] fails"<<endl;
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
      cout << "GL2 integral homology: "; show_invariants(invariants[0]); cout << endl;
      cout << "SL2 integral homology: "; show_invariants(invariants[1]); cout << endl;


      cout<<"----------------------------------------------------------------------------------\n";
    }
}
