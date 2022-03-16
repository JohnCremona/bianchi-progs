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

      // Now we find the intersection of ker(T(P^2)-N(P)) for P
      // running over primes with nontrivial genus class

      // Or (when n2r>1) we could consider over all unramified
      // quadratic characters chi and do this for P such that
      // chi(P)=-1

      // Loop over all Quadprimes P, while dim(s)>0, with a bound on the number of P
      // skip if P.divides(N)
      // skip if P.genus_class()==0
      // m = h.s_calcop(HeckePOp(P, N), 0, 0);
      // s = subeigenspace(m, den*I2long(P.norm()), s);

      QuadprimeLooper Pi(N);
      int ip = 0, np = 10;
      int subdim = dim(s);
      while (ip<np && subdim>0 && Pi.ok())
        {
          Quadprime P = Pi;
          if (P.genus_class()!=0)
            {
              ip++;
              Qideal P2 = P*P;
              matop op;
              if (P2.is_principal())
                op = HeckeP2Op(P, N);
              else   // compute T(P^2)*T(A,A)
                {
                  Qideal A = P.equivalent_coprime_to(N, 1);
                  op = HeckeP2ChiOp(P,A,N);
                }
              m = h.s_calcop(op, 0, 0);
              s = subeigenspace(m, -den*I2long(P.norm()), s);
              subdim = dim(s);
            }
          ++Pi;
        }
      cout << subdim << endl;

    }       // end of while()
  exit(0);
}       // end of main()

