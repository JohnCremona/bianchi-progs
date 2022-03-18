// DIMTABTWIST.CC  -- Table of dimensions: cuspidal
//                                         cuspidal with trivial character
//                                         cuspidal with trivial character and self-twist

#include "matprocs.h"
#include "homspace.h"
#include "oldforms.h"

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
  cout << "#   showing full (respectively old, new) dimensions of the cupidal space and its splitting\n";
  cout << "#   into subspaces with no self-twist and self-twist by each unramified quadratic character.\n\n";

  cout << "# Field\t\tLevel\tdim cusp\tdim cusp triv char\t" <<Quad::all_disc_factors << endl;
  cout << "\t\t\t\t\tall\t\told\t\tnew" <<endl;
  // keys: ideal labels for levels M
  // values: new dimension at level M, by self-twist character
  map<string, vector<int>> newdimlists;

  Qidealooper loop(firstn, lastn, 1, 1); // sorted within norm
  while( loop.not_finished() )
    {
      N = loop.next();
      string Nlabel = ideal_label(N);
      cout << field_label() << "\t" << Nlabel << "\t"; // field, level
      homspace h(N,plusflag,cuspidal,0);  //level, plusflag, cuspidal, verbose

      int cdim = h.h1dim();
      cout << cdim << "\t\t";

      int dimtriv = h.trivial_character_subspace_dimension(cuspidal);
      cout << dimtriv << " ";

      // Now we find the intersection of ker(T(P^2)-N(P)) for P
      // running over primes with genus character chi(P)=-1 for each
      // nontrivial unramified quadratic character chi. The number of
      // these is 2^n2r-1.

      // Here dimlist[0] is the dimension of the subspace with no
      // self-twist (by an unramified quadratic character)

      vector<int> dimlist = h.trivial_character_subspace_dimension_by_twist(cuspidal);

      // full dimensions:
      cout << dimlist<< "\t";

      // Compute old dimensions from proper divisors and subtract:
      vector<int> olddims(nchi,0);
      vector<int> newdims = dimlist;
      vector<Qideal> DD = alldivs(N);
      for(vector<Qideal>::iterator Di = DD.begin(); Di!=DD.end(); ++Di)
        {
          Qideal D = *Di;
          if (D==N)
            continue;
          vector<int> olddimsD = old_multiplicities(D, newdimlists[ideal_label(D)], N);
          for (int i=0; i<nchi; i++)
            {
              olddims[i] += olddimsD[i];
              newdims[i] -= olddimsD[i];
            }
        }

      // old dimensions:
      int olddim = std::accumulate(olddims.begin(), olddims.end(), 0, std::plus<int>());
      cout << olddim << " " << olddims << "\t";

      // new dimensions:
      int newdim = std::accumulate(newdims.begin(), newdims.end(), 0, std::plus<int>());
      cout << newdim << " " << newdims << endl;

      // store new dimensions for use in multiple levels:
      newdimlists[Nlabel] = newdims;

    }       // end of while() level loop
  exit(0);
}       // end of main()

