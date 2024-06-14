#include "swan_sigmas.h"
#include "swan_alphas.h"
#include "swan_tess.h"
#include "swan_hom.h"
#include "swan.h"
#include "geometry.h"

#define MAX_DISC 100

#define VERBOSE 1 // verbose setting to use if not overridden locally
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
// New code for finding alphas and sigmas using SwanData class
//
////////////////////////////////////////////////////////////////////////////////////////////////

      // Create SwanData object
      verbose = VERBOSE;
      debug = DEBUG;
      cout<<"----------------------------------------------------------------------------------\n";
      cout << "Creating SwanData object"<<endl;
      method = "SwanData";
      tim.start(method);
      SwanData SD(1); // 1 means show times for each step
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
      std::system(("mkdir -p "+subdir).c_str());
      string geodata_filename = subdir+"/geodata_"+std::to_string(Quad::d)+".dat";
      cout<<"SwanData A- and S- data output to " << geodata_filename << endl;
      SD.output_alphas_and_sigmas(0,subdir);
      // tim.showAll();

      auto alphas = SDalphas;
      auto sigmas = SDsigmas;

////////////////////////////////////////////////////////////////////////////////////////////////
//
// Use alphas and sigmas to find tessellation
//
////////////////////////////////////////////////////////////////////////////////////////////////

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

////////////////////////////////////////////////////////////////////////////////////////////////
//
// New code for finding tessellating polyhedra using SwanData class
//
////////////////////////////////////////////////////////////////////////////////////////////////

      verbose = VERBOSE;
      cout << "\nFinding all polyhedra using SwanData object" << endl;
      SD.make_all_polyhedra(verbose);


////////////////////////////////////////////////////////////////////////////////////////////////
//
// Compare polyhedra found using old and new code
//
////////////////////////////////////////////////////////////////////////////////////////////////

      if (sing_polys==SD.singular_polyhedra)
        cout<<"Singular polyhedra agree!\n";
      else
        {
          cout<<"********************************Discrepancy in singular polyhedra! (d="<<Quad::d<<")\n";
          exit(1);
        }
      if (princ_polys==SD.principal_polyhedra)
        cout<<"Principal polyhedra agree!\n";
      else
        {
          cout<<"********************************Discrepancy in principal polyhedra! (d="<<Quad::d<<")\n";
          exit(1);
        }
      if (all_polys==SD.all_polyhedra)
        cout<<"All polyhedra agree!\n";
      else
        cout<<"********************************Discrepancy in all_polyhedra!\n";

////////////////////////////////////////////////////////////////////////////////////////////////
//
// Get faces and M32 (= matrix of delta: 3-cells --> 2-cells) from polyhedra
//
////////////////////////////////////////////////////////////////////////////////////////////////

      cout << "\nFinding all faces up to GL2-equivalence" << endl;
      verbose = VERBOSE;
      SD.make_all_faces(verbose);

      // Encode and (if verbose) report on faces found, and check
      // their encodings for consistency:
      verbose = VERBOSE;
      int ok = SD.encode_all_faces(1, verbose);
      if (ok)
        cout<<"all encodings check out OK" << endl;
      else
        {
          cout<<"*****************not all encodings check out OK" << endl;
          exit(1);
        }

      verbose = VERBOSE;
      if (to_file||to_screen) cout << "geodata encodings of faces";
      if (to_file)
        cout<<" output to " << geodata_filename;
      if (to_file||to_screen) cout << "\n";
      if (to_file)
        SD.output_face_data(subdir, to_screen);

      // Test that reading in the data is consistent:

      cout<<"\n***************** reading data back in from "<<geodata_filename<<" ************************\n"<<endl;

      Quad::setup_geometry(subdir, 0);
      int all_ok = 1;
      int sigmas_ok = (SD.slist == Quad::SD.slist);
      cout << "sigmas " << (sigmas_ok? "agree" : "DO NOT agree:") << endl;
      if (!sigmas_ok)
        {
          cout << "sigmas output: " << SD.slist << endl;
          cout << "sigmas input:  " << Quad::SD.slist << endl;
          all_ok = 0;
        }
      int alphas_ok = (SD.alist == Quad::SD.alist);
      cout << "alphas " << (alphas_ok? "agree" : "DO NOT agree:") << endl;
      if (!alphas_ok)
        {
          cout << "alphas output: " << SD.alist << endl;
          cout << "alphas input:  " << Quad::SD.alist << endl;
          all_ok = 0;
        }
      int mats_ok = (SD.Mlist == Quad::SD.Mlist);
      cout << "M_alphas " << (mats_ok? "agree" : "DO NOT agree:") << endl;
      if (!mats_ok) all_ok = 0;
      int a_inv_ok = (SD.a_inv == Quad::SD.a_inv);
      cout << "alpha_inv " << (a_inv_ok? "agree" : "DO NOT agree:") << endl;
      if (!a_inv_ok) all_ok = 0;
      int a_flip_ok = (SD.a_flip == Quad::SD.a_flip);
      cout << "alpha_flip " << (a_flip_ok? "agree" : "DO NOT agree:") << endl;
      if (!a_flip_ok) all_ok = 0;
      int s_flip_ok = (SD.s_flip == Quad::SD.s_flip);
      cout << "sigma_flip " << (s_flip_ok? "agree" : "DO NOT agree:") << endl;
      if (!s_flip_ok)
        {
          cout << "sigma_flip output: " << SD.s_flip << endl;
          cout << "sigma_flip input:  " << Quad::SD.s_flip << endl;
          all_ok = 0;
        }
      int epp_ok = (SD.edge_pairs_plus == Quad::SD.edge_pairs_plus);
      cout << "edge_pairs_plus " << (epp_ok? "agree" : "DO NOT agree:") << endl;
      if (!epp_ok) all_ok = 0;
      int epm_ok = (SD.edge_pairs_minus == Quad::SD.edge_pairs_minus);
      cout << "edge_pairs_minus " << (epm_ok? "agree" : "DO NOT agree:") << endl;
      if (!epm_ok) all_ok = 0;
      int ef_ok = (SD.edge_fours == Quad::SD.edge_fours);
      cout << "edge_fours " << (ef_ok? "agree" : "DO NOT agree:") << endl;
      if (!ef_ok) all_ok = 0;

      all_ok = all_ok && alphas_ok;

      if (!all_ok)
        {
          cout << "!!!!!!!!!!!!!!!!! Inconsistency in data output and input!" << endl;
          exit(1);
        }

      // Compute integral homology

      cout<<"\n***************** integral homology for d = "<<Quad::d<<" ************************\n"<<endl;

      debug = DEBUG;
      if (0)
        {
          cout<<"alphas: "<<alphas<<endl;
          cout<<"sigmas: "<<sigmas<<endl;
          cout<<"pluspairs: "<<pluspairs<<endl;
          cout<<"minuspairs: "<<minuspairs<<endl;
          cout<<"fours: "<<fours<<endl;
          cout<<"faces: "<<SD.all_faces<<endl;
        }
      vector<vector<int>> invariants = SD.integral_homology(3, debug);
      cout << Quad::disc << " GL2 integral homology: "; show_invariants(invariants[0]); cout << endl;
      cout << Quad::disc << " SL2 integral homology: "; show_invariants(invariants[1]); cout << endl;

      vector<vector<int>> invariants_old = integral_homology(SD.all_faces,
                                                         alphas, sigmas,
                                                         pluspairs, minuspairs, fours,
                                                         3, debug);
      cout << "GL2 integral homology (old): "; show_invariants(invariants_old[0]); cout << endl;
      cout << "SL2 integral homology (old): "; show_invariants(invariants_old[1]); cout << endl;
      if (invariants != invariants_old)
        cout<<"******* disagreement in integral homology (d="<<Quad::d<<") **********"<<endl;
      cout<<"----------------------------------------------------------------------------------\n";
    }
}
