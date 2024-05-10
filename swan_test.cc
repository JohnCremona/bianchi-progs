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

  long f, maxpnorm=100;
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
      int verbose = VERBOSE;
      int debug = DEBUG;
      int to_file=1;
      int to_screen=0;

      if ((f==0) && (D>MAX_DISC))
        break;
      long d = (D%4==0? D/4: D);
      Quad::field(d,maxpnorm);

      if (D!=fields.front())
        cout << "-------------------------------------" <<endl;

      if (verbose)
        Quad::displayfield(cout);
      else
        cout << "Field Q(sqrt("<<-d<<"))\tdiscriminant = "<<D<<endl;

      //test_principal_cusps(20,20);

////////////////////////////////////////////////////////////////////////////////////////////////
//
// Old code for finding alphas and sigmas using swan_alphas and swan_sigmas, not SwanData class
//
////////////////////////////////////////////////////////////////////////////////////////////////

      cout << "Finding sigmas and alphas (old method)"<<endl;

      timer tim;
      string method = "old-alpha-sigma-finder";
      tim.start(method);

      verbose = VERBOSE;
      debug = DEBUG;
      string step = "finding sigmas";
      tim.start(step);
      auto new_sigmas = singular_points();
      auto sorted_sigmas = sort_singular_points(new_sigmas);
      if (debug)
        cout<<"Before sorting, "<<new_sigmas.size()<<" sigmas:\n"<<new_sigmas<<endl;
      if (debug||verbose)
        cout<<"After  sorting, "<<sorted_sigmas.size()<<" sigmas:\n"<<sorted_sigmas<<endl;
      new_sigmas = sorted_sigmas;
      tim.stop(step);
      if (debug||verbose) tim.show(1, step);

      //test_principal_cusps(20, 30);
      verbose = VERBOSE;
      debug = DEBUG;
      step = "covering alphas";
      tim.start(step);
      auto new_alphas = covering_alphas(new_sigmas, verbose);
      INT maxn = max_dnorm(new_alphas);
      if (verbose)
        {
          cout << new_alphas.size() << " covering alphas";
          if (debug) cout << ": " << new_alphas || "\n";
          cout << " with max dnorm = " << maxn <<endl;
        }
      tim.stop(step);
      if (debug||verbose) {cout<<step<<"...done: "; tim.show(1, step);}

      verbose = VERBOSE;
      debug = DEBUG;
      step = "saturating alphas";
      tim.start(step);
      new_alphas = saturate_covering_alphas(new_alphas, new_sigmas, maxn, debug, verbose);
      maxn = max_dnorm(new_alphas);
      if (verbose)
        cout << new_alphas.size() << " saturated alphas, max denom norm = " << maxn <<endl;
      if (debug || (verbose && new_alphas.size()<20))
        cout << new_alphas << endl;
      tim.stop(step);
      if (debug||verbose) {cout<<step<<"...done: "; tim.show(1, step);}

      step = "intersection points";
      tim.start(step);
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
      tim.stop(step);
      if (debug||verbose) {cout<<step<<"...done: "; tim.show(1, step);}

      // Sort and output alphas:

      step = "alpha sorting";
      tim.start(step);
      vector<vector<Quad>> pluspairs, minuspairs, fours;
      verbose = VERBOSE;
      debug = DEBUG;
      auto sorted_alphas = sort_alphas(new_alphas, pluspairs, minuspairs, fours, verbose, debug);
      if (debug)
        cout<<"Before sorting, alphas:\n"<<new_alphas<<endl;
      if (debug||verbose)
        cout<<"After  sorting, alphas:\n"<<sorted_alphas<<endl;
      new_alphas = sorted_alphas;
      tim.stop(step);
      if (debug||verbose) {cout<<step<<"...done: "; tim.show(1, step);}
      tim.stop(method);
      cout<<method<<"...done: "; tim.show(1, method);
      // tim.showAll();
      output_alphas(pluspairs, minuspairs, fours, to_file, to_screen);
      output_singular_points(new_sigmas, to_file, to_screen);

////////////////////////////////////////////////////////////////////////////////////////////////
//
// New code for finding alphas and sigmas using SwanData class (not swan_alphas and swan_sigmas)
//
////////////////////////////////////////////////////////////////////////////////////////////////

      // Create SwanData object
      verbose = 1; // VERBOSE;
      debug = DEBUG;
      cout<<"----------------------------------------------------------------------------------\n";
      cout << "Creating SwanData object"<<endl;
      method = "SwanData";
      tim.start(method);
      SwanData SD(verbose); // 1 means show times for each step
      cout << "Using SwanData object to create sigmas:"<<endl;
      auto SDsigmas = SD.get_sigmas();
      cout << SDsigmas.size() << " sigmas found by SwanData: "<<SDsigmas<<endl;
      cout << "Using SwanData object to create (covering, saturated) alphas:"<<endl;
      auto SDalphas = SD.get_alphas(debug);
      cout << SDalphas.size() << " alphas found by SwanData"<<endl;
      auto SDcorners = SD.get_corners();
      cout << SDcorners.size() << " corners found by SwanData"<<endl;
      tim.stop(method);
      cout<<"...done: "; tim.show(1, method);
      string subdir = "xgeodata";
      cout<<"SwanData A- and S- data output to " << subdir << "/geodata_" << Quad::d << ".dat" << endl;
      SD.output_alphas_and_sigmas(0,subdir);
      // tim.showAll();

      new_alphas = SDalphas; // overwrite the SD lists for comparison with stored data
      new_sigmas = SDsigmas; //

////////////////////////////////////////////////////////////////////////////////////////////////
//
// Compare alphas and sigmas found using old and new code
//
////////////////////////////////////////////////////////////////////////////////////////////////

      if (D<1000)
        {
      // Compare with precomputed alphas and sigmas
      cout << "Testing newly computed sigmas and alphas with old..." <<flush;
      Quad::setup_geometry(); // this sets lots of globals including alphas and sigmas, and M_alphas
      cout << "read in old data..."<<flush;
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
              Quad temp;
              for (const auto& a : sigmas)
                if (cusp_index_with_translation(a,new_sigmas,temp)==-1)
                  cout << a << " " << a.coords(1) << endl;
              cout << "new not in old:\n";
              for (const auto& a : new_sigmas)
                if (cusp_index_with_translation(a,sigmas,temp)==-1)
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
              // if (verbose)
                cout << ": " <<alphas;
              cout << endl;
              cout << new_alphas.size() << " new alphas";
              // if (verbose)
                cout << ": " <<new_alphas;
              cout << endl;
              cout << "old not in new:\n";
              Quad temp;
              for (const auto& a : alphas)
                if (cusp_index_with_translation(a,new_alphas,temp)==-1)
                  cout << a << " " << a.coords(1) << endl;
              cout << "new not in old:\n";
              for (const auto& a : new_alphas)
                if (cusp_index_with_translation(a,alphas,temp)==-1)
                  cout << a << " " << a.coords(1) << endl;
              exit(1);
            }
        }
        } // end of code to compare new and old

      // continue; // while testing SwanData only

////////////////////////////////////////////////////////////////////////////////////////////////
//
// Use alphas and sigmas to find tesseallation
//
////////////////////////////////////////////////////////////////////////////////////////////////

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
          auto sigma_nbrs = sorted_neighbours(sigma, alphas);
          auto n = sigma_nbrs.size();
          cout << sigma << " has " << n << " neighbours (sorted): " << sigma_nbrs << endl;
        }
#endif

      // Find all principal polyhedra:

      verbose = VERBOSE;
      cout << "Constructing principal polyhedra from alphas";
      if (debug) cout << ": "<<alphas;
      cout << "..." << flush;
      if (verbose) cout<<endl;
      vector<POLYHEDRON> princ_polys = principal_polyhedra(alphas, verbose);
      int npp = princ_polys.size();
      if (npp==1)
        cout << "done: 1 principal polyhedron constructed:"<<endl;
      else
        cout << "done: " << npp << " principal polyhedra constructed:"<<endl;
      map<string,int> poly_counts;
      for (const auto& P: princ_polys)
        poly_counts[poly_name(P)]++;
      for (const auto& pc : poly_counts)
        cout<<pc.second<<" "<<pc.first << (pc.second>1?"s":"") << endl;

      // Find all singular polyhedra:

      verbose = VERBOSE;
      cout << "Constructing singular polyhedra..."<<flush;
      vector<POLYHEDRON> sing_polys = singular_polyhedra(sigmas, alphas, verbose);
      int nsp = sing_polys.size();
      if (nsp==1)
        cout << "done: 1 singular polyhedron constructed" <<endl;
      else
        cout << "done: " << nsp << " singular polyhedra constructed" <<endl;
      map<string,int> spoly_counts;
      for (const auto& P: sing_polys)
        spoly_counts[poly_name(P)]++;
      for (const auto& pc : spoly_counts)
        cout<<pc.second<<" "<<pc.first << (pc.second>1?"s":"") << endl;

      // Merge principal and singular polyhedra:

      vector<POLYHEDRON> all_polys = sing_polys;
      all_polys.insert(all_polys.end(), princ_polys.begin(), princ_polys.end());


      // Get faces and M32 (= matrix of delta: 3-cells --> 2-cells) from polyhedra:

      cout << "\nFinding all faces up to GL2-equivalence" << endl;
      vector<vector<int>> M32;
      vector<int> redundant_faces;
      verbose = VERBOSE;
      auto all_faces = get_faces(all_polys, alphas, sigmas, M32, redundant_faces, verbose);
      //int nfaces = all_faces.size();

      verbose = VERBOSE;

      // Split up faces into 4 types for reporting and output:

      vector<CuspList> aaa, aas, sqs, hexs;
      cout << "Faces up to GL2-action and reflection:\n";
      int sing, i=0;
      for (const auto& face: all_faces)
        {
          sing = is_face_singular(face, sigmas);
          int red = std::find(redundant_faces.begin(), redundant_faces.end(), i)!=redundant_faces.end();
          int n = face.size();
          string s;
          if (n==4)
            {
              if (!red) sqs.push_back(face);
              s = "square";
              if (red) s+= " (redundant)";
            }
          else
            {
              if (n==6)
                {
                  if (!red) hexs.push_back(face);
                  s = "hexagon";
                  if (red) s+= " (redundant)";
                }
              else
                {
                  if (sing)
                    {
                      if (!red) aas.push_back(face);
                      s = "aas triangle";
                      if (red) s+= " (redundant)";
                    }
                  else
                    {
                      if (!red) aaa.push_back(face);
                      s = "aaa triangle";
                      if (red) s+= " (redundant)";
                    }
                }
            }
          if (verbose)
            cout<<i<<" ("<<s<<"): "<<face<<endl;
          i++;
        }

      // Report on faces found, and check their encodings for consistency:

      verbose = VERBOSE;
      int all_ok = 1;
      cout<<aaa.size()<<" aaa-triangles\n";
      for ( const auto& face : aaa)
        {
          if (verbose) cout <<face << " --> ";
          POLYGON P = make_polygon(face, alphas, sigmas, sing);
          if (verbose) cout <<face << " -->  ["<<P.indices<<","<<P.shifts<<"]"<<endl;
          int ok = check_aaa_triangle(P, verbose);
          if (!ok)
            cout<<"aaa-triangle "<<face<<" --> ["<<P.indices<<","<<P.shifts<<"] fails"<<endl;
          all_ok = ok &&all_ok;
        }
      cout<<aas.size()<<" aas-triangles\n";
      for ( const auto& face : aas)
        {
          POLYGON P = make_polygon(face, alphas, sigmas, sing);
          int ok = check_aas_triangle(P, verbose);
          if (!ok)
            cout<<"aas-triangle "<<face<<" --> ["<<P.indices<<","<<P.shifts<<"] fails"<<endl;
          all_ok = ok &&all_ok;
        }
      cout<<sqs.size()<<" squares\n";
      for ( const auto& face : sqs)
        {
          POLYGON P = make_polygon(face, alphas, sigmas, sing);
          int ok = check_square(P);
          if (!ok)
            cout<<"square "<<face<<" --> ["<<P.indices<<","<<P.shifts<<"] fails"<<endl;
          all_ok = ok &&all_ok;
        }
      cout<<hexs.size()<<" hexagons\n";
      for ( const auto& face : hexs)
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
      output_faces({aaa, sqs, hexs, aas},
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
