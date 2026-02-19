// FILE RNFD.CC:  Read and display precomputed Newspaces (d-dimensional newforms)
/////////////////////////////////////////////////////////////////////////////////

//#define LOOPER
#ifdef LOOPER
#include "qidloop.h"
#endif
#include "nfd.h"

#define MAXPRIME 10000

int main()
{
  cout << "Program rnfd: read and display precomputed Bianchi newforms of arbitrary dimension." << endl;
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
  int n2r = Quad::class_group_2_rank;
  int C4 = is_C4();
  int triv_char_only = !C4;

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
      cout << endl;

      string Nlabel = label(N);
      cout << ">>>> Level " << Nlabel <<" = "<<gens_string(N)<<", norm = "<<N.norm()<<" <<<<" << endl;

      Newspace NS;
      int verbose = 0;
      NS.input_from_file(N, verbose);
      int nnf = NS.nforms();
      cout << nnf << (nnf==1? " newform" : " newforms");
      if (n2r) cout << " up to unramified quadratic twist";
      int nnf_triv_char = NS.nforms_triv_char();
      if (n2r)
        cout << ", of which " << nnf_triv_char
             << (nnf_triv_char==1? " has": " have")
             << " trivial character";
      cout << endl <<endl ;

      int show_aP = 1; // do display aP
      int show_AL = 1; // do display AL
      int show_princ = 0; // do not display principal eigs
      int show_traces = 0; // do not display traces
      NS.display_newforms(show_aP, show_AL, show_princ, show_traces, triv_char_only);

      if (n2r && nnf_triv_char)
        {
          NS.add_unram_quadratic_twists();
          cout << "Full eigensystems for forms with trivial character"
               << ", including unramified quadratic twists"
               << endl<<endl;
          NS.display_newforms(show_aP, show_AL, show_princ, show_traces, 1);
        }
    }     // end of level loop
  cout << endl;
  exit(0);
}   // end of main()
