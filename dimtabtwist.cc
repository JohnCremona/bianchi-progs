// DIMTABTWIST.CC  -- Table of dimensions: cuspidal
//                                         cuspidal with trivial character
//                                         cuspidal with trivial character and self-twist

// This version only handles one unramified quadratic character, hence
// is useful only for fields where the 2-rank of the class group is 1.

#include "matprocs.h"
#include "homspace.h"
#define LOOPER

#define MAXPRIME 10000

int main(void)
{
  long d, max(MAXPRIME);
  Quad n;
  int plusflag=1, cuspidal=1;
  cerr << "Enter field (one of "<<valid_fields<<"): " << flush;  cin >> d;
  if (!check_field(d))
    {
      cerr<<"field must be one of: "<<valid_fields<<endl;
      exit(1);
    }
  cout << "# Table of dimensions of weight 2 Bianchi cuspidal spaces for GL2 over Q(sqrt(-"<<d<<"))" << endl;
  cout << "# Field\tLevel\t";
  cout << "dim(cuspidal)\tdim(cuspidal, trivial unramified character)\tdim(cuspidal self-twist)" << endl;
  Quad::field(d,max);
  Quad::displayfield(cout);
  int n2r = Quad::class_group_2_rank;
  Qideal N;
#ifdef LOOPER
  QUINT firstn, lastn;
  cerr<<"Enter first and last norm for Quad loop: ";
  cin >> firstn >> lastn;
  Qidealooper loop(firstn, lastn, 1, 1); // sorted within norm
  while( loop.not_finished() )
    {
      N = loop.next();
#else
  while(cerr<<"Enter level (ideal label or generator): ", cin>>N, !N.is_zero())
    {
#endif
      cout << "\t"<< field_label() << "\t";     // field
      cout << ideal_label(N)<<"\t"; // level
      homspace h(N,plusflag,cuspidal,0);  //level, plusflag, cuspidal, verbose
      int cdim = h.h1dim();
      int den = h.h1cdenom();

      cout << cdim << "\t";

      if (n2r == 0)
        {
          cout << "0\t0" << endl;
          continue;
        }

      vector<Qideal> nulist = make_nulist(N);
      vector<Qideal>::iterator nui = nulist.begin();
      smat m = h.s_calcop(CharOp(*nui++, N), 0, 0);
      ssubspace s = eigenspace(m, den);
      int dimtriv = dim(s);

      for (; nui!=nulist.end() && dimtriv>0; ++nui)
        {
          m = h.s_calcop(CharOp(*nui, N), 0, 0);
          s = subeigenspace(m, den, s);
          dimtriv = dim(s);
        }

      cout << dimtriv << "\t";
      if (n2r == 0)
        {
          cout << "0" << endl;
          continue;
        }

      cout << endl;

    }       // end of while()
  exit(0);
}       // end of main()

