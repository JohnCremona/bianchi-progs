// FILE P1Ntest.cc  -- for testing the P1N class

#include <iostream>
#include "looper.h"
#include "qidloop.h"
#include "mat22.h"
#include "primes.h"
#include "P1N.h"

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

void symb_index_test(Qideal N, int verbose=0)
{
  P1N P1(N);
  P1.check(verbose);
}

void init()
{
  long d;
  cout << "Enter field: " << flush;  cin >> d;
  Quad::field(d);
  Quad::displayfield(cout);
}

int main(void)
{
  cout << endl << "P1N TEST PROGRAM" << endl;
  init();
  cout << "testing split/merge of multi-indices..." << endl;
  multi_index_test({2,4,2});
  multi_index_test({1,5,6});
  int both=1, sorted=1;
  long maxn=100;
  cout << "testing P1(N) symbol-index bijections for ideals of norm up to "<<maxn<<"..." << endl;
  Qidealooper loop(1, maxn, both, sorted);
  while( loop.not_finished() )
    {
      symb_index_test(loop.next(), 0);
    }
  cout << "testing P1(N) symbol-index bijections for Quads of norm up to "<<maxn<<"..." << endl;
  Quadlooper alpha(1, maxn, both);
  while( alpha.ok() )
    {
      symb_index_test((Quad)alpha, 0);
      ++alpha;
    }
  cout<<"done"<<endl;
}

// END OF FILE P1Ntest.cc
