// FILE EDGE_RELATIONS.CC: Implemention of the edge relations for class homspace

#include "mat22.h"
#include "ratquads.h"
#include "edge_relations.h"
#include <assert.h>

// 2-term (edge) relations

edge_relations::edge_relations(P1N* s, int plus, int verb)
  :P1(s), plusflag(plus), verbose(verb)
{
  int field = Quad::d;
  nsymb = P1->size();
  nsymbx = nsymb*n_alphas;
  ngens=0;
  coordindex.resize(nsymbx);
  gens.reserve(1+nsymbx);  //NB start of gens array is at 1 not 0
  gens.push_back(0);

  if(verbose)
    {
      cout << "About to start on 2-term (edge) relations.\n";
      cout<<"Edge relations for type 0 symbols (denominator 1)\n";
    }
  edge_relations_1();

  if (Quad::is_Euclidean)
    {
      if (verbose)
        report();
      return;
    }

  if(verbose)
    cout<<"Edge relations for type 1,2 symbols (denominator 2)\n";
  edge_relations_2();
  if (field<43)
    {
      if (verbose)
        report();
      return;
    }

  if(verbose)
    cout<<"General edge pair relations\n";
  for (vector<int>::const_iterator i=edge_pairs.begin(); i!=edge_pairs.end(); i++)
    edge_pairing(*i);
  if(verbose)
    cout<<"General edge quadruple relations\n";
  for (vector<int>::const_iterator i=edge_fours.begin(); i!=edge_fours.end(); i++)
    edge_pairing_double(*i);

  if (verbose)
    report();
}

void edge_relations::report()
{
  cout << "After 2-term relations, ngens = "<<ngens<<endl;
  cout << "gens = ";
  int i, j;
  for (i=1; i<=ngens; i++) cout << gens[i] << " ";
  cout << endl;
  cout << "coordindex = \n";
  for (j=0; j<n_alphas; j++)
    {
      int offset = j*nsymb;
      if(n_alphas>1)
        {
          mat22 M = M_alphas[j];
          RatQuad alpha(-M.entry(1,1), M.entry(1,0));
          cout << "Type " << j << ", alpha = "<< alpha<<", M_alpha = "<<M<<":\n";
        }
      Quad c, d;
      for (i=0; i<nsymb; i++)
        {
          P1->make_symb(i, c, d);
          cout << i<<":\t("<<c<<":"<<d<<")\t"<<coordindex[i+offset] << "\n";
        }
    }
  cout << endl;
}

void edge_relations::edge_relations_1()    // basic edge relations for alpha = 0
{
  Quad unit = fundunit;
  long lenrel = Quad::nunits;
  if(!plusflag) {unit=fundunit*fundunit; lenrel/=2;}
  action eps(P1,unit,0,0,1);  assert (eps.det()==unit);
  action sof(P1, mat22::S);
  vector<int> a(lenrel), b(lenrel);
  vector<int> done(nsymb, 0);
  int j, k, triv;
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
              ++ngens;
              gens.push_back(j);
              for(k=0; k<lenrel; k++)
                {
                  coordindex[a[k]] =  ngens;
                  coordindex[b[k]] = -ngens;
                }
            }
        }
    }
}

// edge relations for alphas & sigmas with denominator 2
void edge_relations::edge_relations_2()
{
  int d = Quad::d;
  switch (d%4) {
  case 1:
  case 2:
    edge_relations_2_d12mod4();
    return;
  case 3:
  default:
    switch (d%8) {
    case 3:
      edge_relations_2_d3mod8();
      return;
    case 7:
    default:
      edge_relations_2_d7mod8();
      return;
    } // switch on d%8
  } // switch on d%4
}

// edge relations for one alpha, one sigma with denominator 2 when 2 is ramified d%4=1,2 (d>2)
void edge_relations::edge_relations_2_d12mod4()
{
  Quad w = Quad::w;
  Quad a=w, b=w;
  ((Quad::d)%4==1? b: a) += 1;

  action M(P1, M_alphas[1]);
  action L(P1, -1,a,0,1); assert (L.det()==-1);
  action K(P1, -1,b,0,1); assert (K.det()==-1);

  // alpha#1 = w/2 (d%4=1), (w+1)/2 (d%4=2)

  vector<int> done(nsymb, 0);
  int offset = nsymb;
  int i, m, l, k;
  for (i=0; i<nsymb; i++)
    {
      if (done[i])
        continue;
      m = M(i);
      l = L(i);
      k = M(l); // = L(m)
      done[i] = done[k] = done[l] = done[m] = 1;

      if ((i==m) || (plusflag&&(i==k))) // i==m iff l==k
        {
          coordindex[offset + i] = 0;
          coordindex[offset + l] = 0;
        }
      else
        {
          ++ngens;
          gens.push_back(offset+i);
          coordindex[offset + i] = ngens;
          coordindex[offset + m] = -ngens;
          if (!plusflag)
            {
              ++ngens;
              gens.push_back(offset+l);
            }
          coordindex[offset + l] = ngens;
          coordindex[offset + k] = -ngens;
        }
    }

  // sigma#1 = (1+w)/2  (d%4=1), w/2 (d%4=2)

  if (!plusflag)
    return;

  std::fill(done.begin(), done.end(), 0);
  offset = n_alphas*nsymb;
  for (i=0; i<nsymb; i++)
    {
      if (done[i])
        continue;
      k = K(i);
      done[i] = done[k] = 1;
      ++ngens;
      gens.push_back(offset+i);
      coordindex[offset + i] = ngens;
      coordindex[offset + k] = ngens;
    }
}

// edge relations for two sigmas with denominator 2 whenever 2 is split, i.e. d%8=7 (d>7)
void edge_relations::edge_relations_2_d7mod8()
{
  // sigma#1 = w/2,  sigma#2 = (w+1)/2
  for (int t=0; t<2; t++)
    {
      action L(P1, -1, Quad::w + t, 0,1);
      vector<int> done(nsymb, 0);
      int i, l, offset = nsymb*(n_alphas+t);
      for (i=0; i<nsymb; i++)
        {
          if (done[i])
            continue;
          l = L(i);
          done[i] = done[l] = 1;
          ++ngens;
          gens.push_back(offset+i);
          coordindex[offset + i] = ngens;
          if (!plusflag)
            {
              ++ngens;
              gens.push_back(offset+l);
            }
          coordindex[offset + l] = ngens;
        }
    }
}

// edge relations for two alphas with denominator 2 whenever 2 is inert, i.e. d%8=3 (d>3)
void edge_relations::edge_relations_2_d3mod8()
{
  Quad w = Quad::w;
  int j, k, l, m;

  // relevant alphas are  {1:w/2, 2:(w-1)/2}

  action K(P1, M_alphas[1]);
  action L(P1, -1,w-1,0,1); assert (L.det()==-1);

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
          ++ngens;
          gens.push_back(offset1+k);
          coordindex[offset1 + k] = ngens;
          coordindex[offset2 + j] = -ngens;

          // if plusflag=0 we have a new gen, unless k=m
          if (!plusflag && k!=m)
            {
              ++ngens;
              gens.push_back(offset1+m);
            }
          coordindex[offset1 + m] = ngens;
          coordindex[offset2 + l] = -ngens;
        }
    }
}

// For use when alpha[i]=r/s with r^2=-1 (mod s)

void edge_relations::edge_pairing(int i)
{
  int j, k, j2, k2;
  int offset_a = i*nsymb, offset_b = (i+1)*nsymb;
  action J(P1, mat22::J);
  action M(P1, M_alphas[i]);
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
              ++ngens;
              gens.push_back(offset_a+j);
              coordindex[offset_a+j] = ngens;
              coordindex[offset_a+k] = -ngens;
              if (!plusflag)
                {
                  ++ngens;
                  gens.push_back(offset_b+j2);
                }
              coordindex[offset_b+j2] = ngens;
              coordindex[offset_b+k2] = -ngens;
            }
        }
    }
}

// For use with alpha[i]=r1/s, alpha[i+1]=-r1/s, alpha[i+2]=r2/s, alpha[i+3]=-r2/s, where r1*r2=-1 (mod s)

void edge_relations::edge_pairing_double(int i)
{
  int j, k, j2, k2;
  int offset_a = i*nsymb, offset_b = (i+1)*nsymb, offset_c = (i+2)*nsymb, offset_d = (i+3)*nsymb;

  // M has det 1, maps {alpha[i+2],oo} to {oo,  alpha[i]}
  action M(P1, M_alphas[i+2]);
  action J(P1, mat22::J);

  for (j=0; j<nsymb; j++) // index of type i symbol
    {
      k = M(j); // index of type i+2 symbol: (M)_{i+2} = - (I)_i
      j2 = J(j); // index of type i+1 symbol: (I)_I = (J)_{i+1} if plusflag
      k2 = J(k); // index of type i+3 symbol
      ++ngens;
      gens.push_back(offset_a+j);
      coordindex[offset_a+j] = ngens;
      coordindex[offset_c+k] = -ngens;
      if (!plusflag)
        {
          ++ngens;
          gens.push_back(offset_b+j2);
        }
      coordindex[offset_b+j2] = ngens;
      coordindex[offset_d+k2] = -ngens;
    }
}
