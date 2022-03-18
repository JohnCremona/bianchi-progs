// DIMTABTWIST.CC  -- Table of dimensions: cuspidal
//                                         cuspidal with trivial character
//                                         cuspidal with trivial character and self-twist

// This version only handles one unramified quadratic character, hence
// is useful only for fields where the 2-rank of the class group is 1.

#include "matprocs.h"
#include "homspace.h"

#define MAXPRIME 10000

int main(void)
{
  long d, max(MAXPRIME);
  //  Quad n;
  int plusflag=1, cuspidal=1;
  Qideal N;
  long firstn, lastn;

  cerr << "Enter field (one of "<<valid_fields<<"): " << flush;  cin >> d;
  if (!check_field(d))
    {
      cerr<<"field must be one of: "<<valid_fields<<endl;
      exit(1);
    }
  cerr<<"Enter first and last norm for Quad loop: ";
  cin >> firstn >> lastn;

  Quad::field(d,max);
  int n2r = Quad::class_group_2_rank;
  int nchi = (1<<n2r);

  cout << "# Table of dimensions of weight 2 Bianchi cuspidal spaces for GL2 over Q(sqrt(-"<<d<<"))\n";
  cout << "# Field\t\tLevel\tdim\tdim\t\tdim(s)\n";
  cout << "# \t\t\tcusp\tcusp\t\tcusp\n";
  cout << "# \t\t\t    \ttriv char\ttriv char,\n";
  cout << "# \t\t\t    \t         \tself-twist by\n";
  cout << "# \t\t\t    \t         \t";
  cout<<"no "<<setw(4);
  for (auto Di=Quad::all_disc_factors.begin()+1; Di!=Quad::all_disc_factors.end(); Di++)
    cout<<(*Di)<<" ";
  cout << "any" << endl;

  Qidealooper loop(firstn, lastn, 1, 1); // sorted within norm
  while( loop.not_finished() )
    {
      N = loop.next();
      cout << field_label()<<"\t"<<ideal_label(N)<<"\t"; // level
      homspace h(N,plusflag,cuspidal,0);  //level, plusflag, cuspidal, verbose

      int cdim = h.h1dim();
      cout << cdim << "\t";

      if (n2r == 0 || cdim==0)
        {
          cout << "0\t\t0   ";
          for(int i=0; i<nchi; i++)
            cout<<"0   ";
          cout << endl;
          continue;
        }

      int dimtriv = h.trivial_character_subspace_dimension(cuspidal);
      cout << dimtriv << "\t";
      if (dimtriv == 0)
        {
          cout << "\t0   ";
          for (int i=0; i<nchi; i++)
            cout << "0   ";
          cout << endl;
          continue;
        }

      // Now we find the intersection of ker(T(P^2)-N(P)) for P
      // running over primes with genus character chi(P)=-1 for each
      // nontrivial unramified quadratic character chi. The number of
      // these is 2^n2r-1.

      // Here dimlist[0] is the dimension of the subspace with no
      // self-twist (by an unramified quadratic character)

      vector<int> dimlist = h.trivial_character_subspace_dimension_by_twist(cuspidal);

      for(auto di = dimlist.begin(); di!=dimlist.end(); ++di)
        {
          if (di==dimlist.begin())
            cout << "\t";
          cout << *di << "   ";
        }
      int stdim = std::accumulate(dimlist.begin()+1, dimlist.end(), 0, std::plus<int>());
      cout << stdim << endl;
    }       // end of while()
  exit(0);
}       // end of main()

