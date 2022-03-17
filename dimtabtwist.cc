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
  QUINT firstn, lastn;

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
  vector<QUINT> Dlist;
  for(int chi_index=1; chi_index<nchi; chi_index++)
    {
      QUINT D(1);
      for (int i=0; i<n2r; i++)
        if (bit(chi_index,i)==1)
          D *= Quad::discfactors[i];
      if (D>0)
        D = Quad::disc/D;
      Dlist.push_back(D);
    }

  cout << "# Table of dimensions of weight 2 Bianchi cuspidal spaces for GL2 over Q(sqrt(-"<<d<<"))\n";
  cout << "# Field\t\tLevel\tdim\tdim\t\tdim(s)\n";
  cout << "# \t\t\tcusp\tcusp\t\tcusp\n";
  cout << "# \t\t\t    \ttriv char\ttriv char,\n";
  cout << "# \t\t\t    \t         \tself-twist by\n";
  cout << "# \t\t\t    \t         \t";
  for (auto Di=Dlist.begin(); Di!=Dlist.end(); Di++)
    cout<<setw(3)<<(*Di)<<" ";
  cout << "any" << endl;

  Qidealooper loop(firstn, lastn, 1, 1); // sorted within norm
  while( loop.not_finished() )
    {
      N = loop.next();
      cout << field_label()<<"\t"<<ideal_label(N)<<"\t"; // level
      homspace h(N,plusflag,cuspidal,0);  //level, plusflag, cuspidal, verbose
      int cdim = h.h1dim();
      int den = h.h1cdenom();

      cout << cdim << "\t";

      if (n2r == 0 || cdim==0)
        {
          cout << "0\t0" << endl;
          continue;
        }

      ssubspace s = h.trivial_character_subspace(cuspidal,0); // not dual
      int dimtriv = dim(s);
      cout << dimtriv << "\t\t";
      if (dimtriv == 0)
        {
          for (int i=0; i<nchi; i++)
            cout << "0   ";
          cout << endl;
          continue;
        }

      // Now we find the intersection of ker(T(P^2)-N(P)) for P
      // running over primes with genus character chi(P)=-1 for each
      // nontrivial unramified quadratic character chi. The number of
      // these is 2^n2r-1.

      int stdim = 0;
      for(auto Di = Dlist.begin(); Di!=Dlist.end(); ++Di)
        {
          QUINT D = *Di;
          ssubspace sD = s;
          QuadprimeLooper Pi(N); // loop over primes not dividing N
          int ip = 0, np = 10;   // only use 10 non-square-class primes
          int subdim = dimtriv;
          while (ip<np && subdim>0 && Pi.ok())
            {
              Quadprime P = Pi;
              if (P.genus_character(D) == -1)
                {
                  ip++;
                  long Pnorm = I2long(P.norm());
                  long eig = -den*Pnorm;
                  Qideal P2 = P*P;
                  matop op;
                  if (P2.is_principal())
                    op = HeckeP2Op(P, N);
                  else   // compute T(P^2)*T(A,A)
                    {
                      Qideal A = P.equivalent_coprime_to(N, 1);
                      op = HeckeP2ChiOp(P,A,N);
                    }
                  smat m = h.s_calcop(op, 0, 0); // not dual, no display
                  sD = subeigenspace(m, eig, sD);
                  subdim = dim(sD);
                }
              ++Pi;
            }
          stdim += subdim;
          // cout << subdim << "("<<D<<")   ";
          cout << subdim << "   ";
        }
      cout << stdim << endl;
    }       // end of while()
  exit(0);
}       // end of main()

