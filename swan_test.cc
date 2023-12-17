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
  output_singular_points(sigmas, 1, 1);
  return sigmas;
}

void test_principal_cusps(int n1=20, int n2=100)
{
  for (int n=1; n<=n1; n++)
    {
      auto alphas = principal_cusps_of_norm(n);
      cout << "Principal cusps of denominator norm "<<n<<": "<<alphas<<endl;
    }
  auto alphas = principal_cusps_up_to(n2);
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

      if (D!=fields[0])
        cout << "-------------------------------------" <<endl;
      cout << "The field is ";
      Quad::displayfield(cout);
      cout << endl;

      auto sigmas = test_singular_points(1);
      test_principal_cusps(20, 100);
      auto alphas = covering_alphas(sigmas, 1);
    }
}
