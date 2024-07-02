#include "geometry.h"
#include "swan_tess.h" // for remake_triangle
#include "swan_hom.h"  // for integral_homology
#include "swan_alphas.h"  // for integral_homology
#include "swan.h"  // for integral_homology

#define MAX_DISC 100
#define MIN_DISC 0

#define VERBOSE 1 // verbose setting to use if not overridden locally
#define DEBUG 0   // verbose setting to use if not overridden locally

//#define COMPARE_OLD

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

      debug = DEBUG;
      verbose = VERBOSE;

      if (verbose)
        cout << "Reading data from geodata/geodata_"<<d<<".dat if possible, else creating from scratch..." <<flush;
      Quad::setup_geometry("geodata", verbose);
      auto alphas = Quad::SD.get_alphas();
      auto sigmas = Quad::SD.get_sigmas();
      if (verbose)
        cout << "done: "<<alphas.size()<<" alphas, "<<sigmas.size()<<" sigmas"<<endl;

      // Compute integral homology

      debug = DEBUG;
      verbose = VERBOSE;
      if (verbose)
        cout << "Computing integral homology..." <<endl;

      vector<vector<int>> invariants = Quad::SD.integral_homology(3, debug);
#ifdef COMPARE_OLD
      cout << "NEW CODE" << endl;
#endif
      cout << -D << " GL2 integral homology: "; show_invariants(invariants[0]); cout << endl;
      cout << -D << " SL2 integral homology: "; show_invariants(invariants[1]); cout << endl;

#ifdef COMPARE_OLD
      vector<vector<POLYGON>> all_polys = read_polygons("geodata", verbose);
      int nT = all_polys[0].size(), nU = all_polys[1].size(), nQ = all_polys[2].size(), nH = all_polys[3].size();
      vector<CuspList> all_faces(nT+nU+nQ+nH);
      std::transform(all_polys[0].begin(), all_polys[0].end(), all_faces.begin(),
                     [alphas, sigmas] (const POLYGON& P) {return remake_triangle(P, alphas, sigmas, 0);});
      std::transform(all_polys[1].begin(), all_polys[1].end(), all_faces.begin()+nT,
                     [alphas, sigmas] (const POLYGON& P) {return remake_triangle(P, alphas, sigmas, 1);});
      std::transform(all_polys[2].begin(), all_polys[2].end(), all_faces.begin()+nT+nU,
                     [alphas] (const POLYGON& P) {return remake_quadrilateral(P, alphas);});
      std::transform(all_polys[3].begin(), all_polys[3].end(), all_faces.begin()+nT+nU+nQ,
                     [alphas] (const POLYGON& P) {return remake_hexagon(P, alphas);});
      static const CuspList tri0 = {{0,0,1}, {1,0,0}, {1,0,1}}; // {0,oo,1}
      all_faces.push_back(tri0);
      // if ((d==19)||(d==43)||(d==67)||(d==163))
      //   {
      //     // {w/2,oo,(w-1)/2}
      //     CuspList tri1 = {{0,1,2}, {1,0,0}, {-1,1,2}};
      //     all_faces.push_back(tri1);
      //   }
      // if (d==2)
      //   {
      //     // {0,oo,w,w/2}
      //     CuspList sq = {{0,0,1}, {1,0,0}, {0,1,1}, {Quad(-1),Quad::w}};
      //     all_faces.push_back(sq);
      //   }
      // if (d==7)
      //   {
      //     // {0,oo,w-1,(w-1)/2}
      //     CuspList sq = {{0,0,1}, {1,0,0}, {-1,1,1}, {Quad(-1),Quad::w}};
      //     all_faces.push_back(sq);
      //   }
      // if (d==11)
      //   {
      //     // {0,oo,w-1,2(w-1)/3,(w-1)/2,(w-1)/3}
      //     CuspList hex = {{0,0,1}, {1,0,0}, {-1,1,1},  {Quad(-2),Quad::w}, {-1,1,2}, {Quad(-1),Quad::w}};
      //     all_faces.push_back(hex);
      //   }

      vector<vector<Quad>> pluspairs, minuspairs, fours;
      auto sorted_alphas = sort_alphas(alphas, pluspairs, minuspairs, fours, verbose, debug);
      vector<vector<int>> old_invariants = integral_homology(all_faces, alphas, sigmas,
                                                             pluspairs, minuspairs, fours,
                                                             3, debug);

      cout << "OLD CODE" << endl;
      cout <<  -D << " GL2 integral homology: "; show_invariants(old_invariants[0]); cout << endl;
      cout <<  -D << " SL2 integral homology: "; show_invariants(old_invariants[1]); cout << endl;

      if (invariants==old_invariants)
        {
          cout << "Old and new agree!" << endl;
        }
      else
        {
          cout << "!!!!!!!!!!!!! Old and new disagree!" << endl;
        }
#endif
      cout<<"----------------------------------------------------------------------------------\n";
    }
}
