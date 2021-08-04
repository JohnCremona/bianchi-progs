// FILE EDGE_RELATIONS.CC: Implemention of the edge relations for class homspace

#include "mat22.h"
#include "ratquads.h"
#include "edge_relations.h"
#include <assert.h>

// The base points are the initial points alpha or sigma of the base
// edges {alpha,oo} or {sigma,0}.  Here the alpha are principal cusps
// with alpha[0]=0 while the sigmas (if any) are singular points.

// NB We always set sigma[0]=oo as a filler, so the base path for type
// t is either {alpha[t],oo} if t>=0 or {sigma[-t],oo} if t<0.  For
// this to work, we cannot have a sigma with index 0.

// NB The indices for coordindex are
//
// i+nsymb*t (0<=i<nsymb) for 0 <= t < n_alphas, with "offset" nsymb*t,
//
// and
//
// i+nsymb*(n_alphas-t-1) (0<=i<nsymb) for 0 < -t < n_sigmas, with "offset" nsymb*(n_alphas-t-1).
//
// The method offset(t) gvies the correct offfset for type t.
//
// For example for d=5 with n_alphas=6 and only one singular point, we
// have sigmas[0]=oo, sigmas[1]=(w+1)/2, so edges of type t=-1 map to
// indices i+nsymb*n_alphas.

RatQuad base_point(int t)
{
  if (t>=0) // alpha
    {
      mat22 M = M_alphas[t];
      return RatQuad(-M.entry(1,1), M.entry(1,0));
    }
  else // sigma
    return sigmas[-t];
}

// 2-term (edge) relations

edge_relations::edge_relations(P1N* s, int plus, int verb)
  :P1(s), plusflag(plus), verbose(verb)
{
  nsymb = P1->size();
  long nsymbx = nsymb*(n_alphas+n_sigmas-1);
  ngens=0;
  coordindex.resize(nsymbx);
  gens.reserve(1+nsymbx);  //NB start of gens array is at 1 not 0
  gens.push_back(0);

  if(verbose)
    {
      cout << "About to start on 2-term (edge) relations.\n";
      if (n_alphas>1)
        {
          cout<<"alphas: ";
          for (int i=0; i<n_alphas; i++)
            {
              mat22 M = M_alphas[i];
              RatQuad alpha(-M.entry(1,1), M.entry(1,0));
              cout<<alpha<<" ";
            }
          cout<<endl;
        }
      if (n_sigmas>1)
        {
          cout<<"sigmas: ";
          for (vector<RatQuad>::const_iterator si = sigmas.begin()+1; si!=sigmas.end(); si++)
            {
              cout<<(*si)<<" ";
            }
          cout<<endl;
        }
      cout<<"Edge relations (denominator 1)\n";
    }
  edge_relations_1();

  if (Quad::is_Euclidean)
    {
      if (verbose)
        report();
      return;
    }

  if(verbose)
    {
      cout << "After denominator 1 relations, ngens = "<<ngens<<endl;
      cout<<"Edge relations (denominator 2)\n";
    }
  edge_relations_2();
  if(verbose)
    {
      cout << "After denominator 2 relations, ngens = "<<ngens<<endl;
    }

  if(!edge_pairs_minus.empty())
    {
      if(verbose)
        cout<<"General edge pair relations (-)\n";
      for (vector<int>::const_iterator i=edge_pairs_minus.begin(); i!=edge_pairs_minus.end(); i++)
        {
          if(verbose) cout<<" pair "<< (*i)<<endl;
          edge_pairing_minus(*i);
        }
      if(verbose)
        cout << "After edge pair (-) relations, ngens = "<<ngens<<endl;
    }

  if(!edge_pairs_plus.empty())
    {
      if(verbose)
        cout<<"General edge pair relations (+)\n";
      for (vector<int>::const_iterator i=edge_pairs_plus.begin(); i!=edge_pairs_plus.end(); i++)
        {
          if(verbose) cout<<" pair "<< (*i)<<endl;
          edge_pairing_plus(*i);
        }
      if(verbose)
        cout << "After edge pair (+) relations, ngens = "<<ngens<<endl;
    }

  if(!edge_fours.empty())
    {
      if(verbose)
        cout<<"General edge quadruple relations\n";
      for (vector<int>::const_iterator i=edge_fours.begin(); i!=edge_fours.end(); i++)
        {
          if(verbose) cout<<" quadruple "<< (*i)<<endl;
          edge_pairing_double(*i);
        }
    }
  if (verbose)
    report();
}

void edge_relations::report()
{
  long i;
  int j, t, off;
  RatQuad alpha;
  string name;
  Quad c, d;

  cout << "After 2-term relations, ngens = "<<ngens<<endl;
  cout << "gens = [";
  for (i=1; i<=ngens; i++)
    cout << gens[i] << " ";
  cout << "]" << endl;
  cout << "coordindex = \n";
  for (j=0; j<n_alphas+n_sigmas-1; j++)
    {
      if (j<n_alphas)
        {
          t = j;
          name = "alpha";
        }
      else
        {
          t = n_alphas-j-1;
          name = "sigma";
        }
      off = offset(t);
      alpha = base_point(t);
      if(n_alphas>1)
        {
          cout << "Type " << t << ", "<<name<<" = "<< alpha;
          if (t>=0)
            cout<<", M_alpha = "<<M_alphas[j];
          cout << "\n";
        }
      for (i=0; i<nsymb; i++)
        {
          P1->make_symb(i, c, d);
          cout << i<<":\t("<<c<<":"<<d<<")\t"<<coordindex[i+off] << "\n";
        }
    }
  cout << endl;
}

void edge_relations::edge_relations_1()    // basic edge relations for alpha = 0
{
  Quad unit = fundunit;
  int lenrel = Quad::nunits;
  if(!plusflag) {unit=fundunit*fundunit; lenrel/=2;}
  action eps(P1,unit,0,0,1);  assert (eps.det()==unit);
  action sof(P1, mat22::S);
  vector<int> a(lenrel), b(lenrel);
  vector<int> done(nsymb, 0);
  long j, k;
  int triv;
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
  long off = offset(1);
  long i, m, l, k;
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
          coordindex[off + i] = 0;
          coordindex[off + l] = 0;
        }
      else
        {
          ++ngens;
          gens.push_back(off+i);
          coordindex[off + i] = ngens;
          coordindex[off + m] = -ngens;
          if (!plusflag)
            {
              ++ngens;
              gens.push_back(off+l);
            }
          coordindex[off + l] = ngens;
          coordindex[off + k] = -ngens;
        }
    }

  // sigma#1 = (1+w)/2  (d%4=1), w/2 (d%4=2)

  if (!plusflag)
    return;

  std::fill(done.begin(), done.end(), 0);
  off = offset(-1);
  for (i=0; i<nsymb; i++)
    {
      if (done[i])
        continue;
      k = K(i);
      done[i] = done[k] = 1;
      ++ngens;
      gens.push_back(off+i);
      coordindex[off + i] = ngens;
      coordindex[off + k] = ngens;
    }
}

// edge relations for two sigmas with denominator 2 whenever 2 is split, i.e. d%8=7 (d>7)
void edge_relations::edge_relations_2_d7mod8()
{
  // sigma#1 = w/2,  sigma#2 = (w+1)/2
  for (int t=0; t<2; t++) // types -1, -2, i.e. -1-t
    {
      action L(P1, -1, Quad::w + t, 0,1);
      vector<int> done(nsymb, 0);
      long i, l, off = offset(-1-t);
      for (i=0; i<nsymb; i++)
        {
          if (done[i])
            continue;
          l = L(i);
          done[i] = done[l] = 1;
          ++ngens;
          gens.push_back(off+i);
          coordindex[off + i] = ngens;
          if (!plusflag)
            {
              ++ngens;
              gens.push_back(off+l);
            }
          coordindex[off + l] = ngens;
        }
    }
}

// edge relations for two alphas with denominator 2 whenever 2 is inert, i.e. d%8=3 (d>3)
void edge_relations::edge_relations_2_d3mod8()
{
  Quad w = Quad::w;
  long j, k, l, m;

  // relevant alphas are  {1:w/2, 2:(w-1)/2}

  action K(P1, M_alphas[1]);
  action L(P1, -1,w-1,0,1); assert (L.det()==-1);

  // K maps w/2 --> oo --> (w-1)/2, where K = [w-1,u;2,-w], det=1,  order 3
  // so (g)_(w-1/2) = {g((w-1)/2),g(oo)} = {gK(oo),gK(w/2)} = -(gK)_w/2.
  //
  // Also (g)_(w-1)/2 = (gL)_(w-1)/2      with L = [-1,w-1;0,1],  det=-1, order 2, if plus
  //                  = -(gLK)_w/2

  vector<int> done(nsymb, 0);
  long off1 = offset(1), off2 = offset(2);
  for (j=0; j<nsymb; j++) // index of a type 2 symbol
    {
      if (!done[j])
        {
          k = K(j);      // index of type 1 symbol
          l = L(j);      // index of type 2 symbol
          m = K(l);      // index of type 1 symbol
          done[j] = done[l] = 1;
          ++ngens;
          gens.push_back(off1+k);
          coordindex[off1 + k] = ngens;
          coordindex[off2 + j] = -ngens;

          // if plusflag=0 we have a new gen, unless k=m
          if (!plusflag && k!=m)
            {
              ++ngens;
              gens.push_back(off1+m);
            }
          coordindex[off1 + m] = ngens;
          coordindex[off2 + l] = -ngens;
        }
    }
}

// For use when alpha[i]=r/s with r^2=-1 (mod s)

void edge_relations::edge_pairing_minus(int i)
{
  long j, k, j2, k2;
  long off1 = offset(i), off2 = offset(i+1);
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
              coordindex[off1+j] = 0;
              coordindex[off2+j2] = 0;
            }
          else
            {
              ++ngens;
              gens.push_back(off1+j);
              coordindex[off1+j] = ngens;
              coordindex[off1+k] = -ngens;
              if (!plusflag)
                {
                  ++ngens;
                  gens.push_back(off2+j2);
                }
              coordindex[off2+j2] = ngens;
              coordindex[off2+k2] = -ngens;
            }
        }
    }
}

// For use when alpha[i]=r/s with r^2=+1 (mod s)

void edge_relations::edge_pairing_plus(int i)
{
  long i1, j1, j2, i2;
  long off1 = offset(i), off2 = offset(i+1);
  action J(P1, mat22::J);
  action M(P1, M_alphas[i]);
  vector<int> done(nsymb, 0);

  for (j1=0; j1<nsymb; j1++)
    {
      if (done[j1])
        continue;
      i1 = M(j1);
      j2 = J(j1);
      i2 = J(i1);
      done[j1] = done[i1] = 1;
      ++ngens;
      gens.push_back(off1+j1);
      coordindex[off1+j1] = ngens;
      coordindex[off1+i1] = -ngens;
      if (!plusflag)
        {
          ++ngens;
          gens.push_back(off2+j2);
        }
      coordindex[off2+j2] = ngens;
      coordindex[off2+i2] = -ngens;
    }
}

// For use with alpha[i]=r1/s, alpha[i+1]=-r1/s, alpha[i+2]=r2/s, alpha[i+3]=-r2/s, where r1*r2=-1 (mod s)

void edge_relations::edge_pairing_double(int i)
{
  long j, k, j2, k2;
  long off1 = offset(i), off2 = offset(i+1), off3 = offset(i+2), off4 = offset(i+3);

  // M has det 1, maps {alpha[i+2],oo} to {oo,  alpha[i]}
  action M(P1, M_alphas[i+2]);
  action J(P1, mat22::J);

  for (j=0; j<nsymb; j++) // index of type i symbol
    {
      k = M(j); // index of type i+2 symbol: (M)_{i+2} = - (I)_i
      j2 = J(j); // index of type i+1 symbol: (I)_I = (J)_{i+1} if plusflag
      k2 = J(k); // index of type i+3 symbol
      ++ngens;
      gens.push_back(off1+j);
      coordindex[off1+j] = ngens;
      coordindex[off3+k] = -ngens;
      if (!plusflag)
        {
          ++ngens;
          gens.push_back(off2+j2);
        }
      coordindex[off2+j2] = ngens;
      coordindex[off4+k2] = -ngens;
    }
}
