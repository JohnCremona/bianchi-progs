// FILE EDGE_RELATIONS.CC: Implemention of the edge relations for class homspace

#include <eclib/msubspace.h>
#include <eclib/xmod.h>
#include "homspace.h"
#include "pseuclid.h"
#include <assert.h>

// 2-term (edge) relations

void homspace::edge_relations()    // computes coordindex
{
  int field = Quad::d;

  if(verbose)
    cout<<"Edge relations for type 0 symbols (denominator 1)\n";
  edge_relations_1();
  if (field<19) return;

  if(verbose)
    cout<<"Edge relations for type 1,2 symbols (denominator 2)\n";
  edge_relations_2();
  if (field<43) return;

  if(verbose)
    cout<<"Edge relations for type 3,4,5,6,7,8 symbols (denominator 3)\n";
  edge_relations_3();
  if (field<67) return;

  cout<<"edge relations not yet fully implemented for fields 67, 163" << endl;
}

void homspace::edge_relations_1()    // basic edge relations for alpha = 0
{
  Quad unit = fundunit;
  long lenrel = Quad::nunits;
  if(!plusflag) {unit=fundunit*fundunit; lenrel/=2;}
  symbop eps(this,unit,0,0,1);  assert (eps.det()==unit);
  symbop sof(this, mat22::S);
  vector<int> a(lenrel), b(lenrel);
  vector<int> done(nsymb, 0);
  int j, k, triv;
  if(verbose) cout << "About to start on 2-term (edge) relations.\n";
  if(verbose && n_alphas>1)
    cout<<"Generic edge relations for type 0 symbols\n";
  for (j=nsymb-1; j>=0; j--)
    {
      if (!done[j])
        {
	  if(verbose>1) cout << "j = " << j << ":\t";
          a[0]=j; b[0]=sof(j); triv=(j==b[0]);
          for(k=1; k<lenrel; k++)
            {
              a[k]= eps(a[k-1]);
              b[k]= eps(b[k-1]);
              triv= triv | (j==b[k]);
            }
          for (k=0; k<lenrel; k++) done[a[k]]=done[b[k]]=1;
	  if(verbose>1)
	    {
	      cout<<"+:\t";
	      for (k=0; k<lenrel; k++) cout<<a[k]<<" ";
              cout<<endl;
	      cout<<"\t-:\t";
	      for (k=0; k<lenrel; k++) cout<<b[k]<<" ";
              cout<<endl;
	    }
          if (triv)
            for (k=0; k<lenrel; k++) coordindex[a[k]]=coordindex[b[k]]=0;
          else
            {
              gens[++ngens] = j;
              for(k=0; k<lenrel; k++)
                {
                  coordindex[a[k]] =  ngens;
                  coordindex[b[k]] = -ngens;
                }
            }
        }
    }
}

void homspace::edge_relations_2()    // extra edge relations for alphas with denominator 2
{
  Quad w = Quad::w;
  int j, k, l, m;

  // relevant alphas are  {1:w/2, 2:(w-1)/2}

  symbop K(this, M_alphas[1]);
  symbop L(this, -1,w-1,0,1); assert (L.det()==-1);

  // K maps w/2 --> oo --> (w-1)/2, where K = [w-1,u;2,-w], det=1,  order 3
  // so (g)_(w-1/2) = {g((w-1)/2),g(oo)} = {gK(oo),gK(w/2)} = -(gK)_w/2.
  //
  // Also (g)_(w-1)/2 = (gL)_(w-1)/2      with L = [-1,w-1;0,1],  det=-1, order 2, if plus
  //                  = -(gLK)_w/2

  vector<int> done(nsymb, 0);
  int offset1 = nsymb, offset2 = 2*nsymb;
  for (j=0; j<nsymb; j++) // index of a type 2 symbol
    {
      if (!done[j])
        {
          k = K(j);      // index of type 1 symbol
          l = L(j);      // index of type 2 symbol
          m = K(l);      // index of type 1 symbol
          done[j] = done[l] = 1;
          gens[++ngens] = offset1+k;
          coordindex[offset1 + k] = ngens;
          coordindex[offset2 + j] = -ngens;

          // if plusflag=0 we have a new gen, unless k=m
          if (!plusflag && k!=m)
            {
              gens[++ngens] = offset1+m;
            }
          coordindex[offset1 + m] = ngens;
          coordindex[offset2 + l] = -ngens;
        }
    }
}

void homspace::edge_relations_3()    // extra edge relations for alphas with denominator 3
{
  int j, k, l, m;

  // relevant alphas are  {3:w/3, 4:-w/3, 5:(1-w)/3, 6:(w-1)/3, 7:(1+w)/3, 8:(-1-w)/3}
  int offset7 = 7*nsymb, offset8 = 8*nsymb, offset3 = 3*nsymb, offset4 = 4*nsymb, offset5 = 5*nsymb, offset6 = 6*nsymb;

  // J has det -1, order 2, swaps {a,oo} and {-a,oo}
  symbop J(this, mat22::J);

  // (1) type 7, alpha=(1+w)/3, antisymmetric via M77:

  // M77 has order 2, swaps {(1+w)/3,oo} and its reverse
  symbop M77(this, M_alphas[7]);

  vector<int> done(nsymb, 0);
  for (j=0; j<nsymb; j++)
    {
      if (!done[j])
        {
          k = M77(j);
          done[j] = done[k] = 1;
          if (j==k) // symbol trivial
            {
              coordindex[offset7+j] = 0;
            }
          else
            {
              gens[++ngens] = offset7+j;
              coordindex[offset7+j] = ngens;
              coordindex[offset7+k] = -ngens;
            }
        }
    }

  // (2) type 8, alpha = -(1+w)/3:
  // either (+) pair with type 7's via J
  // or     (0) impose antisymmetry by M88

  if(plusflag)
    {
      for (j=0; j<nsymb; j++)
        {
          k = J(j);
          coordindex[offset8+k] = coordindex[offset7+j];
        }
    }
  else
    {
      // M88 has order 2, det +1, swaps {-(1+w)/3,oo} and its reverse, i.e. alpha_8 -> oo -> alpha_8
      symbop M88(this, M_alphas[8]);
      std::fill(done.begin(), done.end(), 0);
      for (j=0; j<nsymb; j++)
        {
          if (!done[j])
            {
              k = M88(j);
              done[j] = done[k] = 1;
              if (j==k) // symbol trivial
                {
                  coordindex[offset8+j] = 0;
                }
              else
                {
                  gens[++ngens] = offset8+j;
                  coordindex[offset8+j] = ngens;
                  coordindex[offset8+k] = -ngens;
                }
            }
        }
    }

  // (3) types 3,4,5,6: identify in 4-tuples up to sign if (+), else in pairs if (0)

  // M53 has det 1, maps {(1-w)/3,oo} to {oo,  w/3}, i.e. alpha_5 -> oo -> alpha_3
  symbop M53(this, M_alphas[5]);
  // J swaps {alpha_3, alpha_4} and {alpha_5, alpha_6}

  for (j=0; j<nsymb; j++) // index of type 3 symbol
    {
      k = M53(j); // index of type 5 symbol: (M53)_5 = - (I)_3
      m = J(j);   // index of type 4 symbol: (I)_3 = (J)_4 if plsuflag
      l = J(k);   // index of type 6 symbol
      gens[++ngens] = offset3+j;
      coordindex[offset3+j] = ngens;
      coordindex[offset5+k] = -ngens;
      if (!plusflag)
        {
          gens[++ngens] = offset4+m;
        }
      coordindex[offset4+m] = ngens;
      coordindex[offset6+l] = -ngens;
    }
}

