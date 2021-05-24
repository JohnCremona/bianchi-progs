// FILE EDGE_RELATIONS.CC: Implemention of the edge relations for class homspace

#include <eclib/msubspace.h>
#include <eclib/xmod.h>
#include "homspace.h"
#include <assert.h>

// 2-term (edge) relations

void homspace::edge_relations()    // computes coordindex
{
  Quad unit = fundunit;
  long lenrel = Quad::nunits;
  if(!plusflag) {unit=fundunit*fundunit; lenrel/=2;}
  symbop eps(this,unit,0,0,1);  assert (eps.det()==unit);
  symbop sof(this,0,-1,1,0);  assert (sof.det()==1);
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
  int field = Quad::d;
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
  int field = Quad::d;
  Quad w = Quad::w;
  int j, jj, k, l, u=-(field+5)/12; // u=-4, -6, -14 for 43,67,163

  // relevant alphas are  {3:w/3, 4:-w/3, 5:(1-w)/3, 6:(w-1)/3, 7:(1+w)/3, 8:(-1-w)/3}

  // J has det -1, order 2, swaps {a,oo} and {-a,oo}
  symbop J(this, -1,0,0,1);
  assert (J.det()==-1);

  // (1) type 7, alpha=(1+w)/3, antisymmetric via K3:

  // K3 has order 2, swaps {(1+w)/3,oo} and its reverse
  symbop K3(this, 1+w,-(w+u+1),3,-(1+w));
  assert (K3.det()==1);
  int offset7 = 7*nsymb;
  vector<int> done(nsymb, 0);
  for (j=0; j<nsymb; j++)
    {
      if (!done[j])
        {
          done[j] = 1;
          k = K3(j);
          if (j==k) // symbol trivial
            {
              coordindex[offset7+j] = 0;
            }
          else
            {
              done[k] = 1;
              gens[++ngens] = offset7+j;
              coordindex[offset7+j] = ngens;
              coordindex[offset7+k] = -ngens;
            }
        }
    }

  // (2) type 8, alpha = -(1+w)/3:
  // either (+) pair with type 7's via J
  // or     (0) impose antisymmetry by K4

  int offset8 = 8*nsymb;
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
      // K4 = J*K3*K has order 2, det +1, swaps {-(1+w)/3,oo} and its reverse
      symbop K4(this, -(1+w),-(w+u+1),3,1+w);
      assert (K4.det()==1);
      std::fill(done.begin(), done.end(), 0);
      for (j=0; j<nsymb; j++)
        {
          if (!done[j])
            {
              done[j] = 1;
              k = K4(j);
              if (j==k) // symbol trivial
                {
                  coordindex[offset8+j] = 0;
                }
              else
                {
                  done[k] = 1;
                  gens[++ngens] = offset8+j;
                  coordindex[offset8+j] = ngens;
                  coordindex[offset8+k] = -ngens;
                }
            }
        }
    }

  // (3) types 3,4,5,6: identify in 4-tuples up to sign if (+), else in pairs if (0)

  int offset3 = 3*nsymb, offset4 = 4*nsymb, offset5 = 5*nsymb, offset6 = 6*nsymb;
  // M1 has det 1, maps {(1-w)/3,oo} to {oo, w/3}
  symbop M1(this, w,u,3,w-1);
  assert (M1.det()==1);
  // M3 = J*M1*J has det 1, maps {(w-1)/3,oo} to {oo, -w/3}
  symbop M3(this, w,-u,-3,w-1);
  assert (M3.det()==1);

  std::fill(done.begin(), done.end(), 0);
  for (j=0; j<nsymb; j++) // index of type 3 symbol
    {
      if (!done[j])
        {
          jj = J(j); // index of type 4 symbol
          k = M1(j); // index of type 5 symbol
          l = M3(jj); // index of type 6 symbol
          done[j] = done[jj] = 1;
          gens[++ngens] = offset3+j;
          coordindex[offset3+j] = ngens;
          coordindex[offset5+k] = -ngens;
          if (!plusflag)
            {
              gens[++ngens] = offset4+jj;
            }
          coordindex[offset4+jj] = ngens;
          coordindex[offset6+l] = -ngens;
        }
    }
}

