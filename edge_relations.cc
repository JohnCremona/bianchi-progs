// FILE EDGE_RELATIONS.CC: Implemention of the edge relations for class homspace

#include "mat22.h"
#include "ratquads.h"
#include "homspace.h"
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

edge_relations::edge_relations(P1N* s, int plus, int verb, long ch)
  :P1(s), plusflag(plus), verbose(verb), characteristic(ch)
{ //cout<<"In edge_relations constructor with P^1("<<P1->level()<<"), plus="<<plus<<endl;
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
          for ( const auto& s : sigmas)
            if (!s.is_infinity())
              cout<<s<<" ";
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
  edge_relations_2(); // alphas and sigmas with denom 2
  if(verbose)
    {
      cout << "After denominator 2 relations, ngens = "<<ngens<<endl;
    }

  if(!edge_pairs_minus.empty())
    {
      if(verbose)
        cout<<"General edge pair relations (-)\n";
      for ( const auto& e : edge_pairs_minus)
        {
          if(verbose) cout<<" pair "<< e <<flush;
          edge_pairing_minus(e);
          if(verbose) cout<<": ngens now "<< ngens<<endl;
        }
      if(verbose)
        cout << "After edge pair (-) relations, ngens = "<<ngens<<endl;
    }

  if(!edge_pairs_plus.empty())
    {
      if(verbose)
        cout<<"General edge pair relations (+)\n";
      for ( const auto& e : edge_pairs_plus)
        {
          if(verbose) cout<<" pair "<< e <<flush;
          edge_pairing_plus(e);
          if(verbose) cout<<": ngens now "<< ngens<<endl;
        }
      if(verbose)
        cout << "After edge pair (+) relations, ngens = "<<ngens<<endl;
    }

  if(!edge_fours.empty())
    {
      if(verbose)
        cout<<"General edge quadruple relations\n";
      for ( const auto& e : edge_fours)
        {
          if(verbose) cout<<" quadruple "<< e <<flush;
          edge_pairing_double(e);
          if(verbose) cout<<": ngens now "<< ngens<<endl;
        }
    }
  if(verbose)
    cout<<"Singular edge relations...";
  sigma_relations();
  if(verbose)
    {
      cout<<" ngens now "<< ngens<<endl;
      report();
    }
  //  assert(check());
}

void edge_relations::report()
{
  long i;
  int j, t;
  RatQuad alpha;
  string name;
  Quad c, d;

  cout << "After 2-term relations, ngens = "<<ngens<<endl;
  cout << "gens = [";
  for (i=1; i<=ngens; i++)
    cout << gens[i] << " ";
  cout << "]" << endl;
  cout << "coordindex = "<<coordindex<<"\n";
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
      int off = offset(t);
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

int edge_relations::check() // not completely implemented, don't use
{
  if (verbose)
    cout<<"Checking edge relations..."<<endl;
  int i, i2, n;
  pair<long, int> st, st2;
  action J(P1, mat22::J);
  Quad c,d, c2, d2;

  for (i=0; i<(int)coordindex.size(); i++)
    {
      n = coordindex[i];
      if (n==0) continue; // no check for this case yet
      st = symbol_number_and_type(i);
      if (n>0)
        {
          i2 = gens[n];
          if (i2!=i)
            {
              st2 = symbol_number_and_type(i2);
              // the next lines do not yet take into account alphas/sigmas a such that 2*a is integral
              // if(J(st.first) != st2.first)
              //   {
              //     cout<<"i="<<i<<", coordindex[i]="<<n<<" --> i2="<<i2<<", coordindex["<<i2<<"] = "<<coordindex[i2]<<endl;
              //     cout<<"st = ("<<st.first<<","<<st.second<<"), st2 = ("<<st2.first<<","<<st2.second<<")"<<endl;
              //     return 0;
              //   }
              if (st.second>0)
                {
                  if(alpha_flip[st.second] != st2.second)
                    {
                      cout<<"i="<<i<<", coordindex[i]="<<n<<" --> i2="<<i2<<", coordindex["<<i2<<"] = "<<coordindex[i2]<<endl;
                      cout<<"st = ("<<st.first<<","<<st.second<<"), st2 = ("<<st2.first<<","<<st2.second<<")"<<endl;
                      cout<<"alpha_flip["<<st.second<<"] = "<<alpha_flip[st.second]<<endl;
                      return 0;
                    }
                }
              else
                {
                  if(sigma_flip[-st.second] != -st2.second)
                    {
                      cout<<"i="<<i<<", coordindex[i]="<<n<<" --> i2="<<i2<<", coordindex["<<i2<<"] = "<<coordindex[i2]<<endl;
                      cout<<"st = ("<<st.first<<","<<st.second<<"), st2 = ("<<st2.first<<","<<st2.second<<")"<<endl;
                      cout<<"sigma_flip["<<-st.second<<"] = "<<sigma_flip[-st.second]<<endl;
                      return 0;
                    }
                }
            }
        }
      else // (n<0)
        {
          if (st.second<0)
            {
              cout<<"i="<<i<<", coordindex[i]="<<n<<"<0, but ";
              cout<<"st = ("<<st.first<<","<<st.second<<"), with t<0"<<endl;
              return 0;
            }
          i2 = gens[-n];
          st2 = symbol_number_and_type(i2);
          action M(P1, M_alphas[st.second]);
          // the next lines do not yet take into account alphas/sigmas a such that 2*a is integral
          if (!  (((alpha_inv[st.second] == st2.second) && (M(st2.first) == st.first))
                  || ((alpha_flip[alpha_inv[st.second]] == st2.second) && (M(J(st2.first)) == st.first))))
            {
              cout<<"i="<<i<<", coordindex[i]="<<n<<" --> i2="<<i2<<", coordindex["<<i2<<"] = "<<coordindex[i2]<<endl;
              cout<<"i ->  s t = ("<<st.first<<","<<st.second<<"), ";
              cout<<"i2 -> st2 = ("<<st2.first<<","<<st2.second<<"), \n";
              cout<<" alpha_inv["<<st.second<<"] = "<<alpha_inv[st.second];
              cout<<"; M("<<st2.first<<") = "<<M(st2.first)<<endl;
              cout<<" alpha_flip[alpha_inv["<<st.second<<"]] = "<< alpha_flip[alpha_inv[st.second]];
              cout<<"; MJ("<<st2.first<<") = "<<M(J(st2.first))<<endl;
              return 0;
            }
        }
    }
  if (verbose)
    cout<<"... finished checking edge relations..."<<endl;
  return 1;
}

void edge_relations::edge_relations_1()    // basic edge relations for alpha = 0
{
  Quad unit = fundunit, zero(0), one(1);
  int lenrel = Quad::nunits;
  if(!plusflag) {unit=fundunit*fundunit; lenrel/=2;}
  action eps(P1,unit,zero,zero,one);  assert (eps.det()==unit);
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
          if (triv && (characteristic!=2))
            for (k=0; k<lenrel; k++) coordindex[a[k]]=coordindex[b[k]]=0;
          else
            {
              ++ngens;
              gens.push_back(j);
              for(k=0; k<lenrel; k++)
                {
                  coordindex[a[k]] =  ngens;
                  if (!triv) coordindex[b[k]] = -ngens; // only in char.2
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

//#define DEBUG

// edge relations for one alpha, one sigma with denominator 2 when 2 is ramified d%4=1,2 (d>2)
void edge_relations::edge_relations_2_d12mod4()
{
  Quad w = Quad::w, a, b, zero(0), one(1);
  // alphas[1] = a/2, sigmas[1] = b/2
  if ((Quad::d)%4==1)   // for d=1(4), alpha_1=w/2, sigma_1=(w+1)/2
    {
      a = w;
      b = w+ONE;
    }
  else                 // for d=2(4), alpha_1=(w+1)/2, sigma_1=w/2
    {
      a = w+ONE;
      b = w;
    }

  action M(P1, M_alphas[1]);
  action L(P1, -one,a,zero,one); assert (L.det()==-one);

  // alpha_1 = a/2 = w/2 (d%4=1), (w+1)/2 (d%4=2)

  assert(check_rel({mat22::identity, M}, {1,1}));
  assert(check_rel({mat22::identity, L}, {1,1}, {1,-1}));

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

      int triv = (i==m);  // i==m iff l==k since M,L commute, both order 2
      if (triv && characteristic!=2)
        {
          coordindex[off + i] = 0;
          coordindex[off + l] = 0;
        }
      else
        {
          ++ngens;
          gens.push_back(off+i);
          coordindex[off + i] = ngens;
          if (!triv) coordindex[off + m] = -ngens;
#ifdef DEBUG
          cout<<" - increasing ngens to "<<ngens<<" and adding "<<off+i<< " to gens..."<<endl;
          cout<<" - setting coordindex["<<off+i<< "] to "<<ngens<<endl;
          cout<<" - setting coordindex["<<off+m<< "] to "<<-ngens<<endl;
#endif
          if ((i!=l) && (i!=k))
            {
              if (!plusflag)
                {
                  ++ngens;
                  gens.push_back(off+l);
#ifdef DEBUG
                  cout<<" - increasing ngens to "<<ngens<<" and adding "<<off+l<< " to gens..."<<endl;
#endif
                }
              coordindex[off + l] = ngens;
              if (!triv) coordindex[off + k] = -ngens;
#ifdef DEBUG
              cout<<" - setting coordindex["<<off+l<< "] to "<<ngens<<endl;
              cout<<" - setting coordindex["<<off+k<< "] to "<<-ngens<<endl;
#endif
            }
        }
    }

  // sigma_1 = b/2 = (1+w)/2  (d%4=1), w/2 (d%4=2)

  action K(P1, -one,b,zero,one); assert (K.det()==-one);
  assert (check_rel({mat22::identity, K}, {-1,-1}, {1,-1}));

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
#ifdef DEBUG
      cout<<" - increasing ngens to "<<ngens<<" and adding "<<off+i<< " to gens..."<<endl;
      cout<<" - setting coordindex["<<off+i<< "] to "<<ngens<<endl;
#endif
      if (!plusflag && (i!=k))
        {
          ++ngens;
          gens.push_back(off+k);
#ifdef DEBUG
          cout<<" - increasing ngens to "<<ngens<<" and adding "<<off+k<< " to gens..."<<endl;
#endif
        }
      coordindex[off + k] = ngens;
#ifdef DEBUG
      cout<<" - setting coordindex["<<off+k<< "] to "<<ngens<<endl;
#endif
    }
}

// edge relations for sigmas with denominator 2 whenever 2 is split, i.e. d%8=7 (d>7)

void edge_relations::edge_relations_2_d7mod8()
{
  Quad zero(0), one(1), two(2);
  // sigma_1 = w/2,  sigma_2 = (1-w)/2

  for (int t=1; t<3; t++) // types -t = -1, -2
    {
      Quad x;
      (two*sigmas[t]).is_integral(x);
      action L(P1, -one, x, zero,one); // fixes x/2 = sigma
      vector<int> done(nsymb, 0);
      long i, l, off = offset(-t);

      assert (check_rel({mat22::identity, L}, {-t,-t}, {1,-1}));

      for (i=0; i<nsymb; i++)
        {
          if (done[i])
            continue;
          l = L(i);
          done[i] = done[l] = 1;
          ++ngens;
          gens.push_back(off+i);
          coordindex[off + i] = ngens;
          if (l!=i)
            {
              if (!plusflag)
                {
                  ++ngens;
                  gens.push_back(off+l);
                }
              coordindex[off + l] = ngens;
            }
        }
    }
}

// edge relations for two alphas with denominator 2 whenever 2 is inert, i.e. d%8=3 (d>3)
void edge_relations::edge_relations_2_d3mod8()
{
  Quad w = Quad::w, zero(0), one(1);
  long j, k, l, m;

  // relevant alphas are  {1:w/2, 2:(w-1)/2}

  action M(P1, M_alphas[1]);
  action L(P1, -one,w-one,zero,one); assert (L.det()==-one);

  // M maps w/2 --> oo --> (w-1)/2, where M = [w-1,u;2,-w], det=1,  order 3
  // so (g)_(w-1/2) = {g((w-1)/2),g(oo)} = {gM(oo),gM(w/2)} = -(gM)_w/2.
  //
  // Also (g)_(w-1)/2 = (gL)_(w-1)/2      with L = [-1,w-1;0,1],  det=-1, order 2, if plus
  //                  = -(gLM)_w/2

  assert(check_rel({mat22::identity, M}, {2,1}, {1,1}));
  assert(check_rel({mat22::identity, M_alphas[2]}, {1,2}, {1,1}));
  assert(check_rel({mat22::identity, mat22(-one,w,zero,one)}, {1,1}, {1,-1}));
  assert(check_rel({mat22::identity, L}, {2,2}, {1,-1}));

  vector<int> done(nsymb, 0);
  long off1 = offset(1), off2 = offset(2);
  for (j=0; j<nsymb; j++)
    {
      if (!done[j])
        {
          k = M(j);
          l = L(j);
          m = M(l);
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

  assert (check_rel({mat22::identity, M_alphas[i]}, {i, i}));
  assert (check_rel({mat22::identity, M_alphas[i+1]}, {i+1, i+1}));
  assert (check_rel({mat22::identity, mat22::J}, {i, i+1}, {1,-1}));
  assert (check_rel({mat22::identity, mat22::J}, {i+1, i}, {1,-1}));

  for (j=0; j<nsymb; j++)
    {
      if (!done[j])
        {
          k = M(j);
          j2 = J(j);
          k2 = J(k);
          done[j] = done[k] = 1;
          int triv = (j==k); // symbol trivial except in char.2
          if (triv && characteristic!=2)
            {
              coordindex[off1+j] = 0;
              coordindex[off2+j2] = 0;
            }
          else
            {
              ++ngens;
              gens.push_back(off1+j);
              coordindex[off1+j] = ngens;
              if (!triv) coordindex[off1+k] = -ngens;
              if (!plusflag)
                {
                  ++ngens;
                  gens.push_back(off2+j2);
                }
              coordindex[off2+j2] = ngens;
              if (!triv) coordindex[off2+k2] = -ngens;
            }
        }
    }
}

// For use when alpha[i]=r/s with r^2=+1 (mod s)

void edge_relations::edge_pairing_plus(int i)
{
  long j, k, l, m;
  long off1 = offset(i), off2 = offset(i+1);
  action J(P1, mat22::J);
  action M(P1, M_alphas[i]);
  vector<int> done(nsymb, 0);

  assert (check_rel({mat22::identity, M_alphas[i]}, {i+1, i}));
  assert (check_rel({mat22::identity, M_alphas[i+1]}, {i, i+1}));
  assert (check_rel({mat22::identity, mat22::J}, {i, i+1}, {1,-1}));
  assert (check_rel({mat22::identity, mat22::J}, {i+1, i}, {1,-1}));

  for (j=0; j<nsymb; j++)
    {
      if (done[j])
        continue;
      k = M(j);
      l = J(j);
      m = J(k); assert (l==M(m));
      done[j] = done[m] = 1;

      int triv = (k==l && plusflag); // equivalently, j==m
      if (triv && characteristic!=2)
        {
          assert (j==m);
          coordindex[off1+k] = 0;
          coordindex[off2+j] = 0;
        }
      else
        {
          assert (j!=m || !plusflag || characteristic==2);
          ++ngens;
          gens.push_back(off1+l);
          coordindex[off1+l] = ngens;
          coordindex[off2+m] = -ngens;
          if (!plusflag && j!=m)
            {
              ++ngens;
              gens.push_back(off2+j);
            }
          coordindex[off2+j] = ngens;
          coordindex[off1+k] = -ngens;
        }
    }
}

// For use with alpha[i]=r1/s, alpha[i+1]=-r1/s, alpha[i+2]=r2/s, alpha[i+3]=-r2/s, where r1*r2=-1 (mod s)

void edge_relations::edge_pairing_double(int i)
{
  long off1 = offset(i), off2 = offset(i+1), off3 = offset(i+2), off4 = offset(i+3);

  // M has det 1, maps {alpha[i+2],oo} to {oo,  alpha[i]}
  action M(P1, M_alphas[i+2]);
  action J(P1, mat22::J);

  for (long j=0; j<nsymb; j++) // index of type i symbol
    {
      long k = M(j); // index of type i+2 symbol: (M)_{i+2} = - (I)_i
      long j2 = J(j); // index of type i+1 symbol: (I)_I = (J)_{i+1} if plusflag
      long k2 = J(k); // index of type i+3 symbol
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

void edge_relations::sigma_relations()          // for sigma with 2*sigma not integral
{
  action J(P1, mat22::J);
  int i, j, k, l, off1, off2;
  for (i=0; i<n_sigmas; i++)
    {
      j = sigma_flip[i];
      // if i==j,  dealt with in edge_relations_2()
      // if i>j we already dealt with this pair
      if (j<=i) continue;
      // otherwise the symbols of types -i, -j are identified in pairs when plusflag is true
      off1 = offset(-i);
      off2 = offset(-j);
      for (k=0; k<nsymb; k++)
        {
          ++ngens;
          gens.push_back(off1+k);
          coordindex[off1+k] = ngens;
          l = J(k);
          if (!plusflag)
            {
              ++ngens;
              gens.push_back(off2+l);
            }
          coordindex[off2+l] = ngens;
        }
    }
}
