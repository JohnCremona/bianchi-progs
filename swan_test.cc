#include "swan.h"
#include "ratquads.h"
#include "qideal.h"

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
  cout << "Number of singular points, including oo: "<<sigmas.size()<<endl;
  if (output_level>=2)
    cout << "Unsorted singular points: "<<sigmas<<endl;
  sigmas = sort_singular_points(sigmas);
  if (output_level>=1)
    cout << "Sorted singular points: "<<sigmas<<endl;
  int to_file=(output_level>=1);
  int to_screen=(output_level>=2);
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
      cout << endl;

      auto sigmas = test_singular_points(0);
      //test_principal_cusps(20, 30);
      int verbose = 1;
      auto alphas = covering_alphas(sigmas, verbose);
      INT maxn = max_dnorm(alphas);
      cout << alphas.size() << " covering alphas: " << alphas << endl;
      cout << "max dnorm = " << maxn <<endl;

      // auto corners = triple_intersections(alphas, 1);
      // cout << corners << endl;

      int debug = 0;
      alphas = saturate_covering_alphas(alphas, sigmas, maxn, debug, verbose);
      maxn = max_dnorm(alphas);
      cout << alphas.size() << " saturated alphas, max denom norm = " << maxn <<endl;
      if (debug)
        cout << alphas << endl;
      auto points = triple_intersections(alphas);
      RAT minht = points[0].second;
      std::for_each(points.begin(), points.end(),
                    [&minht](H3point P) {minht = min(minht, P.second);});
      cout << points.size() << " vertices, min square height = " << minht <<endl;
    }
}
