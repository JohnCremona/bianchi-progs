// FILE EDGE_RELATIONS.CC: Implemention of the edge relations for class homspace

#include <eclib/msubspace.h>
#include <eclib/xmod.h>
#include "homspace.h"
#include "geometry.h"
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
    cout<<"General edge pair relations\n";
  for (vector<int>::const_iterator i=edge_pairs.begin(); i!=edge_pairs.end(); i++)
    edge_pairing(*i);
  if(verbose)
    cout<<"General edge quadruple relations\n";
  for (vector<int>::const_iterator i=edge_fours.begin(); i!=edge_fours.end(); i++)
    edge_pairing_double(*i);

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

// For use when alpha[i]=r/s with r^2=-1 (mod s)

void homspace::edge_pairing(int i)
{
  int j, k, j2, k2;
  int offset_a = i*nsymb, offset_b = (i+1)*nsymb;
  symbop J(this, mat22::J);
  symbop M(this, M_alphas[i]);
  vector<int> done(nsymb, 0);

  for (j=0; j<nsymb; j++)
    {
      if (!done[j])
        {
          k = M(j);
          j2 = J(j);
          k2 = J(k);
          done[j] = done[k] = 1;
          if (j==k) // symbol trivial
            {
              coordindex[offset_a+j] = 0;
              coordindex[offset_b+j2] = 0;
            }
          else
            {
              gens[++ngens] = offset_a+j;
              coordindex[offset_a+j] = ngens;
              coordindex[offset_a+k] = -ngens;
              if (!plusflag)
                {
                  gens[++ngens] = offset_b+j2;
                }
              coordindex[offset_b+j2] = ngens;
              coordindex[offset_b+k2] = -ngens;
            }
        }
    }
}

// For use wieth alpha[i]=r1/s, alpha[i+1]=-r1/s, alpha[i+2]=r2/s, alpha[i+3]=-r2/s, where r1*r2=-1 (mod s)

void homspace::edge_pairing_double(int i)
{
  int j, k, j2, k2;
  int offset_a = i*nsymb, offset_b = (i+1)*nsymb, offset_c = (i+2)*nsymb, offset_d = (i+3)*nsymb;

  // M has det 1, maps {alpha[i+2],oo} to {oo,  alpha[i]}
  symbop M(this, M_alphas[i+2]);
  symbop J(this, mat22::J);

  for (j=0; j<nsymb; j++) // index of type i symbol
    {
      k = M(j); // index of type i+2 symbol: (M)_{i+2} = - (I)_i
      j2 = J(j); // index of type i+1 symbol: (I)_I = (J)_{i+1} if plusflag
      k2 = J(k); // index of type i+3 symbol
      gens[++ngens] = offset_a+j;
      coordindex[offset_a+j] = ngens;
      coordindex[offset_c+k] = -ngens;
      if (!plusflag)
        {
          gens[++ngens] = offset_b+j2;
        }
      coordindex[offset_b+j2] = ngens;
      coordindex[offset_d+k2] = -ngens;
    }
}
