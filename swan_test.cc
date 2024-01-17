#include "swan.h"
//#include "ratquads.h"
//#include "qideal.h"
#include "geometry.h"

#define MAX_DISC 100

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
      auto alphas = principal_cusps_of_dnorm(n);
      cout << "Principal cusps of denominator norm "<<n<<": "<<alphas<<endl;
    }
  auto alphas = principal_cusps_of_dnorm_up_to(n2);
  cout << "Principal cusps of denominator norm up to "<<n2<<": "<<alphas<<endl;
}

int main ()
{
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
      cout << "The field is ";
      Quad::displayfield(cout);

      //test_principal_cusps(20,20);

      cout << "Finding sigmas and alphas..."<<flush;

      int verbose = 0;
      int debug = 0;
      auto new_sigmas = test_singular_points(0);
      output_singular_points(new_sigmas, 1, verbose);
      //test_principal_cusps(20, 30);
      verbose=0;
      debug=0;
      auto new_alphas = covering_alphas(new_sigmas, verbose);
      INT maxn = max_dnorm(new_alphas);
      if (verbose)
        {
          cout << new_alphas.size() << " covering alphas";
          if (debug) cout << ": " << new_alphas || "\n";
          cout << " with max dnorm = " << maxn <<endl;
        }

      // auto corners = triple_intersections(new_alphas, 1);
      // cout << corners << endl;
      verbose=0;
      debug=0;
      new_alphas = saturate_covering_alphas(new_alphas, new_sigmas, maxn, debug, verbose);
      maxn = max_dnorm(new_alphas);
      if (verbose)
        cout << new_alphas.size() << " saturated alphas, max denom norm = " << maxn <<endl;
      if (debug || (verbose && new_alphas.size()<20))
        cout << new_alphas << endl;
      auto points = triple_intersections(new_alphas);
      RAT minht = points[0].second;
      std::for_each(points.begin(), points.end(),
                    [&minht](H3point P) {minht = min(minht, P.second);});
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
      verbose = 0;
      debug = 0;
      auto sorted_alphas = sort_alphas(new_alphas, pluspairs, minuspairs, fours, verbose, debug);
      output_alphas(pluspairs, minuspairs, fours, 1, verbose);

      cout<<"...done"<<endl;

      // Compare with precomputed alphas and sigmas
      cout << "Testing newly computed sigmas and alphas with old..." <<flush;
      Quad::setup_geometry();
      cout << "read in old data..."<<flush;
      Quad t;
      if (compare_CuspLists_as_sets(sigmas, new_sigmas))
        cout << "sigmas agree..." <<flush;
      else
        {
          if (compare_CuspLists_as_sets_mod_translation(sigmas, new_sigmas))
            cout << "sigmas agree (up to translation)..." <<endl;
          else
            {
              cout << "sigmas DO NOT agree:\n";
              cout << sigmas.size() << " old sigmas";
              if (verbose) cout << ": " <<sigmas;
              cout << endl;
              cout << new_sigmas.size() << " new sigmas";
              if (verbose) cout << ": " <<new_sigmas;
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
    }
}
