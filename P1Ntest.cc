// FILE P1Ntest.cc  -- for testing the P1N class

#include <iostream>
#include "qidloop.h"
#include "mat22.h"
#include "primes.h"
#include "P1N.h"

void multi_index_test(const vector<long>& nlist, int verbose=0)
{
  long i, j, n=1;
  for (i=0; i<(long)nlist.size(); i++) n*=nlist[i];
  if (verbose)
    cout << "nlist = "<<nlist<<" with product "<<n<<endl;
  for (i=0; i<n; i++)
    {
      vector<long> isplit = split_indices(nlist, i);
      j = merge_indices(nlist, isplit);
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
  cout << "Testing P1(N) for N = " << N << ":" <<flush;
  P1N P1(N);
  long psi = P1.size();
  cout << "psi(N) = "<<psi<<endl;
  long i, j;
  for (i=0; i<psi; i++)
    {
      pair<Quad, Quad> cd = P1.symb(i);
      j = P1.index(cd);
      if (verbose)
        {
          cout << i << " --> ("<<cd.first<<":"<<cd.second << ") --> " << j << endl;
          assert(i==j);
        }
      else
        {
          if (i!=j)
            {
              cout << i << " --> ("<<cd.first<<":"<<cd.second << ") --> " << j << endl;
            }
        }
    }
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
  Qidealooper loop(1, 100, both, sorted);
  while( loop.not_finished() )
    {
      Qideal N = loop.next();
      symb_index_test(N);
    }
}

// END OF FILE P1Ntest.cc
