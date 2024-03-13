#include "swan.h"
//#include "ratquads.h"
//#include "qideal.h"
#include "geometry.h"
#include "looper.h"

int main ()
{
  long d=391, max=100, maxdnorm=1580;
  Quad::field(d,max);
  Quad::displayfield(cout);

  int i;

  CuspList sigmas = singular_points();
  cout << "Singular points: "<<sigmas<<endl;
  CuspList alist;
  cout << "Finding principal cusps of dnorm up to "<<maxdnorm<<"..."<<flush;
  auto alist1 = principal_cusps_of_dnorm_up_to(maxdnorm);
  cout << "done: " <<alist1.size()<<" principal cusps"<<endl;
  cout << "Eliminating any whose circles lie inside previous circles..."<<flush;
  for (const auto& a: alist1)
    {
      if (!circle_inside_any_circle(a, alist, 0)) // strict=0
        {
          alist.push_back(a);
        }
    }
  cout << "done: using "<<alist.size()<<" alphas of dnorm up to "<<maxdnorm<<endl;

  cout << "Enter: " << endl;
  cin >> i;
  vector<CuspPair> pairs_ok;
  //  RatQuad a0 = RatQuad(Quad(-32,7), Quad(31,2));

  int allok = 1;
  for (const auto& a0 : alist)
    {
      cout << "testing whether "<<a0<<" is surrounded..."<<flush;
      int ok = is_alpha_surrounded(a0, alist, sigmas, pairs_ok, 1);
      if (ok)
        cout << "YES" << endl;
      else
        {
          cout << "No" << endl;
          exit(1);
        }
      allok = allok && ok;
    }
  if (allok)
    cout << "all ARE surrounded";
  else
    cout << "all are NOT surrounded";
}
