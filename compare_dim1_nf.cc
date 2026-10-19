// FILE COMPARE_DIM1_NF.CC:  Read and compare precomputed (1-dimensional) newforms in newforms/ and newspaces/
//
/////////////////////////////////////////////////////////////////////////////////

//#define LOOPER
#ifdef LOOPER
#include "qidloop.h"
#endif
#include "newspace.h"

#define MAXPRIME 10000

scalar modulus = default_modulus<scalar>();

int main()
{
#ifndef LOOPER
  cout << "Program list_nsp: read and list precomputed Bianchi newforms of arbitrary dimension." << endl;
#endif
  eclib_pari_init();

  long d, maxpnorm(MAXPRIME);
  cerr << "Enter field: " << flush;
  cin >> d;
  if (d<1)
    {
      cout << "Field parameter d must be positive for Q(sqrt(-d))" << endl;
      exit(0);
    }
  Quad::field(d,maxpnorm);

  Quad n;
  Qideal N;

#ifdef LOOPER
  long firstn, lastn;
  cerr<<"Enter first and last norm for levels: ";
  cin >> firstn >> lastn;

  Qidealooper loop(firstn, lastn, 1, 1); // 1 = both conjugates; 1 = sorted within norm
  while( loop.not_finished() )
    {
      N = loop.next();
#else
  while(cerr<<"Enter level (ideal label or generator): ", cin>>N, !N.is_zero())
    {
#endif

      newforms nf(N, modulus, 0);
      nf.read_from_file();
      const vector<newform>& nflist = nf.the_newforms();
      if (nflist.empty())
        continue;

      Newspace NS;
      NS.input_from_file(N, 0);
      const vector<Newform>& NFlist = NS.the_newforms();

      auto  f = nflist.begin();
      auto  F = NFlist.begin();
      while (f!=nflist.end())
        {
          int ok = compare_dim1_newforms(*f, *F, 1);
          if (ok) cout << f->long_label() << ": checks ok" << endl;
          ++f;
          ++F;
        }
    }
      flint_cleanup_master();
      exit(0);
 }   // end of main()
