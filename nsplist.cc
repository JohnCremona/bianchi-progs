// FILE NSPLIST.CC:  Read and compactly display precomputed Newspaces (d-dimensional newforms)
//
// Output format same as nflist.cc (which is only for 1-d newforms)
//
/////////////////////////////////////////////////////////////////////////////////

//#define LOOPER
#ifdef LOOPER
#include "qidloop.h"
#endif
#include "nfd.h"

#define MAXPRIME 10000

// switch on//off different formats for eigenvalue display

const int SHOW_AP_RELATIVE = 0;
const int SHOW_AP_ABSOLUTE = 0;
const int SHOW_AP_INT_COORDS = 1;

int main()
{
  cout << "Program nsplist";
  #ifdef LOOPER
  cout << "_loop";
#endif
  cout << ": read and list precomputed Bianchi newforms of arbitrary dimension." << endl;
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
  Quad::displayfield(cout);
  //  int C4 = is_C4();
  int triv_char_only = 1;
  // cerr << "Newspaces with trivial character only? \n(Otherwise -- for C4 class group only -- also include newspaces with unramified character.  Hence ignored except for C4 fields)" <<endl;
  // cin >> triv_char_only;
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
      Newspace NS;
      NS.input_from_file(N, 0);
      NS.list_newforms(triv_char_only);
    }
      flint_cleanup_master();
      exit(0);
 }   // end of main()
