// DIMTABTWIST.CC  -- Table of dimensions of Bianchi cusp forms at level N with trivial character
//                    showing split by self-twist characters

// Columns:  Field Level homdimcusp homdimcusptriv dimcusp dimcusptrivall [*] dimcusptrivold [*] dimcusptrivnew [*]

// where homdimcusp = dimension of cuspidal homology with quadratic unramified or trivial character
//       homdimcusptriv = dimension of cuspidal homology with trivial character
//       dimcusptrivall = dimension of all Bianchi cuspforms with trivial character
//       dimcusptrivold = dimension of old Bianchi cuspforms with trivial character
//       dimcusptrivnew = dimension of new Bianchi cuspforms with trivial character
//
// and for the last 3 columns, also the subdimension split by self-twist character.

// This is similar to dimtabnew.cc.  Here we also list the full
// cuspidal homology dimension and the trivial character subspace
// homology dimension (only different in even class number).  For the
// Bianchi cuspforms (with trivial character only), we give the full,
// old and new dimensions, and also show the dimension split between
// the non-self-twist and all possible self-twist characters (only
// relevant in even class number)

// Inputs (prompted for): d (field)
//                        min_norm (lower bound on level norm)
//                        max_norm (upper bound on level norm)

// This program requires precomputed data for all levels of norm up to
// max_norm.  The new (cuspidal, trivial character) homology
// dimnesions are read from data files, converted to new Bianchi
// dimensions and stored in a dict (with level labels as keys); the
// old dimensions are computed from the new dimensions at lower
// levels.


#include "matprocs.h"
#include "newforms.h"

#define MAXPRIME 10000

int main(void)
{
  long d, max(MAXPRIME);
  //  Quad n;
  int plusflag=1;
  Qideal N;
  long min_norm=1, max_norm;

  cerr << "Enter field: " << flush;  cin >> d;
  if (!check_field(d))
    {
      cerr<<"field must be one of: "<<valid_fields<<endl;
      exit(1);
    }
  cerr<<"Enter max norm for Quad loop: ";
  cin >> max_norm;

  Quad::field(d,max);

  // Each new eigensystem stored which is stored gives rise to either
  // 2^r Bianchi newforms (if not self-twist) or 2^{r-1} Bianchi
  // newforms (if self-twist by a quadratic genus character), which
  // are all twists of eachother by the group of unramified quadratic
  // twists.
  int n2r = Quad::class_group_2_rank;
  int nchi = 1<<n2r;
  int nchi2 = nchi/2; // only used when nchi is even

  cout << "# Table of dimensions of weight 2 Bianchi cuspidal spaces for GL2 over Q(sqrt(-"<<d<<"))\n";
  cout << "#   showing full (respectively old, new) dimensions of the cuspidal space and its splitting\n";
  cout << "#   into subspaces with no self-twist and self-twist by each unramified quadratic character.\n\n";

  cout << "# Field\t\tLevel\thom-cusp-dim\tBMF-cusp-dim (triv char)" <<Quad::all_disc_factors << endl;
  cout << "\t\t\tall\ttriv\tall\t\told\t\tnew" <<endl;
  // keys: ideal labels for levels M
  // values: new dimension at level M, by self-twist character
  map<string, vector<int>> newdimlists;

  Qidealooper loop(min_norm, max_norm, 1, 1); // sorted within norm
  while( loop.not_finished() )
    {
      N = loop.next();
      string Nlabel = ideal_label(N);
      cout << field_label() << "\t" << Nlabel << "\t"; // field, level

      // compute homology dimensions directly: cuspidal and trivial character cuspidal
      homspace h(N,plusflag,0);  //level, plusflag, verbose
      int cdim = h.h1cuspdim();
      cout << cdim << "\t";
      int dimtriv = h.trivial_character_subspace_dimension(/*cuspidal*/ 1);
      cout << dimtriv << "\t";

      // Read new dimensions from data file:

      vector<int> newdims(nchi,0);
      newforms nfdata(N);
      if (nfdata.read_from_file())
        {
          if (nchi == 1)
            {
              newdims[0] = nfdata.n1ds + nfdata.n2ds;
            }
          else
            {
              // multiplicity nchi for non-self-twist eigensystems
              newdims[0] = nchi*(nfdata.new1dims[0] + nfdata.new2dims[0]);
              // multiplicity nchi/2 for self-twist eigensystems
              for (int i=1; i<nchi; i++)
                newdims[i] = nchi2*(nfdata.new1dims[i] + nfdata.new2dims[i]);
            }
          newdimlists[Nlabel] = newdims;
        }
      else
        {
          cout<<"no newforms file for level " << Nlabel<<"! "<<endl;
          break;
        }

      // Compute old dimensions from proper divisors:

      vector<int> olddims(nchi,0);
      vector<Qideal> DD = alldivs(N);
      for( auto& D : DD)
        {
          if (D==N)
            continue;
          Qideal M = N/D;
          int mult = ndivs(M);
          vector<int> dimsD = newdimlists[ideal_label(D)];
          for (int i=0; i<nchi; i++)
            {
              olddims[i] += mult*dimsD[i];
            }
        }

      // add old+new:

      vector<int> fulldims(nchi,0);
      for (int i=0; i<nchi; i++)
        {
          fulldims[i] = olddims[i] + newdims[i];
        }

      // total dims, not split by self-twist character:

      int newdim = std::accumulate(newdims.begin(), newdims.end(), 0, std::plus<int>());
      int olddim = std::accumulate(olddims.begin(), olddims.end(), 0, std::plus<int>());
      int fulldim = newdim+olddim;

      // full dimensions:
      cout << fulldim<< " " << fulldims << "\t";

      // old dimensions:
      cout << olddim << " " << olddims << "\t";

      // new dimensions:
      cout << newdim << " " << newdims << endl;
    }       // end of while() level loop
  exit(0);
}       // end of main()

