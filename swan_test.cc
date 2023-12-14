#include "swan.h"
#include "ratquads.h"
#include "qideal.h"

#define MAX_DISC 100

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
  output_singular_points(sigmas, 1, 1);
}

int main ()
{
  long d, D, f, max=100;
  vector<long> fields = valid_field_discs();
  cerr << "Enter field (0 for all up to "<<MAX_DISC<<", -1 for all): " << flush;  cin >> f;

  cerr << endl;

  if (f>0)
    fields = {f};
  for (auto di = fields.begin(); di!=fields.end(); ++di)
    {
      D = *di;
      if ((f==0) && (D>MAX_DISC))
        break;
      d = (D%4==0? D/4: D);
      Quad::field(d,max);

      if (di!=fields.begin())
        cout << "-------------------------------------" <<endl;
      cout << "The field is ";
      Quad::displayfield(cout);
      cout << endl;

      test_singular_points(1);
    }
}
