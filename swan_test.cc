#include "swan.h"
#include "ratquads.h"
#include "qideal.h"

void test_singular_points(int output_level=0)
{
  if (output_level>=3)
    for (auto I=Quad::class_group.begin(); I!=Quad::class_group.end(); ++I)
      {
        cout<<"Ideal class ["<<*I<<"]: ";
        cout<<"singular points "<<singular_points_in_class(*I,(output_level>3))<<endl;
      }
  vector<RatQuad> sigmas = singular_points();
  cout << "Number of singular points, including oo: "<<sigmas.size()<<endl;
  if (output_level>=2)
    cout << "Unsorted singular points: "<<sigmas<<endl;
  sigmas = sort_singular_points(sigmas);
  if (output_level>=1)
    cout << "Sorted singular points: "<<sigmas<<endl;
  output_singular_points(sigmas, 0, 1);
}

int main ()
{
  long d, max=100;
  cout << "Enter field: " << flush;  cin >> d;
  //cout << "Enter max. norm for primes: " << flush;  cin >> max;
  Quad::field(d,max);
  cout << "The field is "; Quad::displayfield(cout); cout << endl;

  test_singular_points(3);
}
