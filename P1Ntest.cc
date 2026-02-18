// FILE P1Ntest.cc  -- for testing the P1N class

#include <iostream>
#include "looper.h"
#include "qidloop.h"
#include "mat22.h"
#include "primes.h"
#include "P1N.h"

// Run P1N_test() for all ideals up to this norm bound:
#define MAXN 30
// for fields of absolute discriminant up to this:
#define MAX_DISC 50

void multi_index_test(const vector<long>& nlist, int verbose=0)
{
  long i, n=1;
  for (i=0; i<(long)nlist.size(); i++) n*=nlist[i];
  if (verbose)
    cout << "nlist = "<<nlist<<" with product "<<n<<endl;
  for (i=0; i<n; i++)
    {
      vector<long> isplit = split_indices(nlist, i);
      long j = merge_indices(nlist, isplit);
      if (verbose)
        {
          cout << i << " --> "<<isplit << " --> " << j << endl;
        }
      else
        assert(i==j);
    }
}

void P1N_test(const Qideal& N, int verbose=0)
{
  P1N P1(N);
  P1.check(verbose);
  P1.check_lifts(verbose);
}

int main(void)
{
  cout << endl << "P1N TEST PROGRAM" << endl;

  long f;
  vector<long> fields = valid_field_discs();
  cout << "Enter field (0 for all): " << flush;  cin >> f;
  if (f)
    fields = {f};
  for (auto di = fields.begin(); di!=fields.end(); ++di)
    {
      long D = *di;
      if (D>MAX_DISC)
        break;
      cout << "===============================================================\n";
      long d = (D%4==0? D/4: D);
      Quad::field(d);
      Quad::displayfield(cout);

      cout << "testing split/merge of multi-indices..." << endl;
      multi_index_test({2,4,2});
      multi_index_test({1,5,6});
      int both=1, sorted=1;
      long maxn(MAXN);
      cout << "testing P1(N) symbol-index bijections for ideals of norm up to "<<maxn<<"..." << endl;
      Qidealooper loop(1, maxn, both, sorted);
      while( loop.not_finished() )
        {
          P1N_test(loop.next(), 0);
        }
      cout << "testing P1(N) symbol-index bijections for Quads of norm up to "<<maxn<<"..." << endl;
      Quadlooper alpha(1, maxn, both);
      while( alpha.ok() )
        {
          P1N_test(Qideal(alpha), 0);
          ++alpha;
        }
      cout<<"field " << d << " done"<<endl;
    }
}

// END OF FILE P1Ntest.cc
