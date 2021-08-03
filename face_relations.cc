// FILE FACE_RELATIONS.CC: Implemention of the face relations for class homspace

//#define TIME_RELATION_SOLVING

#include "mat22.h"
#include "ratquads.h"
#include "face_relations.h"
#include <assert.h>
#ifdef TIME_RELATION_SOLVING
#include <eclib/timer.h>
#endif

#ifdef USE_CRT
int liftmats_chinese(const smat& m1, scalar pr1, const smat& m2, scalar pr2, smat& m, scalar& dd);
#endif

// Each face relation is a signed sum of edges (M)_alpha = {M(alpha},
// M(oo)} for M in the list mats and alpha=alphas[t] (when t>=0) or
// sigmas[-t] (when t<0), for t in the list types.  Here we check that such a
// relation holds identically in H_3 (not just modulo the congruence
// subgroup!)

int check_face_rel(const vector<mat22>& mats, const vector<int>& types, const vector<int>& signs);

// Special case: all signs +1
int check_face_rel(const vector<mat22>& mats, const vector<int>& types)
{
  vector<int> signs(mats.size(), 1);
  return check_face_rel(mats, types, signs);
}

// General case:
int check_face_rel(const vector<mat22>& mats, const vector<int>& types, const vector<int>& signs)
{
  vector<mat22>::const_iterator mi;
  vector<int>::const_iterator ti, si;
  vector<RatQuad> alphas, betas;
  RatQuad a, b;
  mat22 M, M_alpha;
  for (mi=mats.begin(), ti=types.begin(), si=signs.begin(); ti!=types.end(); mi++, ti++, si++)
    {
      M = *mi;
      a = M(base_point(*ti));
      b = M(RatQuad(1,0));
      if (*si>0) // use {a,b} when sign is +1
        {
          alphas.push_back(a);
          betas.push_back(b);
        }
      else // use {b,a} when sign is -1
      {
        alphas.push_back(b);
        betas.push_back(a);
      }
    }

  // cout<<"Checking face relation "<<mats<<", "<<types<<endl;
  // cout<<"alphas: "<<alphas<<endl;
  // cout<<"betas:  "<<betas<<endl;

  vector<RatQuad>::const_iterator alpha, beta;
  int ok=1;
  for (alpha=alphas.begin()+1, beta=betas.begin(); beta!=betas.end() &&ok; alpha++, beta++)
    {
      RatQuad next_alpha = (alpha==alphas.end()? alphas[0]: *alpha);
      ok = ok && (*beta==next_alpha);
    }
  if (!ok)
    {
      cout<<"Bad face relation!\n";
      cout<<"alphas: "<<alphas<<endl;
      cout<<"betas:  "<<betas<<endl;
      exit(1);
    }
  // else
  //   {
  //     cout<<"Good face relation "<<endl;
  //   }
  return ok;
}

//
// face relations
//

face_relations::face_relations(edge_relations* er, int plus, int verb)
  :ER(er), plusflag(plus), verbose(verb)
{
  P1 = ER->P1;
  nsymb = P1->size();
  maxnumrel=2*(n_alphas+n_sigmas-1)*nsymb;
  ngens = ER->ngens;
  numrel = 0;
  hmod = 0;

#if(USE_SMATS)
  relmat = smat(maxnumrel,ngens);
#else
  relmat.init(maxnumrel,ngens);
#endif

  make_relations();

  if(verbose)
    {
      cout << "Finished making face relation matrix: ";
      cout << "number of relations = " << numrel;
      cout << " (bound was "<<maxnumrel<<")"<<endl;
    }

  solve_relations();

  if (verbose)
    {
      cout << "Finished solving face relation matrix: ";
      cout << "dimension (relative to cusps) = " << rk << endl;
      if (rk>0)
        {
          if (verbose>1) cout << "coord:" << coord;
          if (hmod)
            cout << "failed to lift, coord is only defined modulo "<<hmod<<endl;
          else
            cout << "lifted ok, denominator = " << denom << endl;
          cout << "pivots = " << pivs <<endl;
        }
    }
}

void face_relations::make_relations()
{
  long field = Quad::d;

  triangle_relation_0();

  if ((field==1)||(field==3))
    {
      if(!plusflag)
        triangle_relation_1_3();
      return;
    }

  if (field==2)
    {
      square_relation_2();
      return;
    }
  if (field==7)
    {
      rectangle_relation_7();
      return;
    }
  if (field==11)
    {
      hexagon_relation_11();
      return;
    }

  if (field==5)
    {
      square_relation_5();
      return;
    }

  // Now field = 19, 43, 67 or 163

  // additional triangle relations
  triangle_relation_2();
  if(verbose) cout<<"\nApplying "<<cyclic_triangles.size()<<" cyclic triangle relations"<<endl;
  for (vector<int>::const_iterator T = cyclic_triangles.begin(); T!=cyclic_triangles.end(); T++)
    {
      cyclic_triangle_relation(*T);
    }
  if(verbose) cout<<"\nApplying "<<triangles.size()<<" general triangle relations"<<endl;
  for (vector<vector<int>>::const_iterator T = triangles.begin(); T!=triangles.end(); T++)
    {
      general_triangle_relation(*T);
    }
  if(verbose) cout<<"\nApplying "<<squares.size()<<" general square relations"<<endl;
  for (vector<pair<vector<int>, vector<Quad>> >::const_iterator S = squares.begin(); S!=squares.end(); S++)
    {
      general_square_relation(S->first, S->second);
    }
}

// In add_face_rel(rel, types, signs):
//
//   rel is a list of (positive) (c:d)-symbol numbers i
//   types is a list of symbol types (negative for types involving a singular base point)
//   signs is a list of +1,-1
//
//   such that the signed sum of the corresponding (symbol,type) is 0 in homology.  We use
//   the map (i,t) -> j = ER.coords(i,t) to convert this to a vector of
//   dimension ngens, and append that as a new row to the relation
//   matrix relmat.
//

// Special case: all signs +1
void face_relations::add_face_rel(const vector<long>& rel, const vector<int>& types)
{
  vector<int> signs(rel.size(), 1);
  add_face_rel(rel, types, signs);
}

// General case:
void face_relations::add_face_rel(const vector<long>& rel, const vector<int>& types, const vector<int>& signs)
{
  vector<long>::const_iterator r;
  vector<int>::const_iterator t, s;
  if (verbose)
    {
      cout<<"Relation: ";
      Quad c, d;
      for (r = rel.begin(), t = types.begin(), s=signs.begin(); r!=rel.end(); r++, t++, s++)
        {
          P1->make_symb(*r, c, d);
          cout<<(*r)<<"_"<<(*t);
          cout<< ((*s)>0? " +": " -");
          cout<<"("<<c<<":"<<d<<") ";
        }
      cout <<" --> ";
    }
#ifdef USE_SMATS
  svec relation(ngens);
#else
  vec relation(ngens);
#endif
  for (r = rel.begin(), t = types.begin(), s=signs.begin(); r!=rel.end(); r++, t++, s++)
    {
      //      cout<<"Looking up edge coord of symbol "<<(*r)<<", type "<<(*t)<<"...";
      long c = (*s) * ER->coords(*r, *t);
      //      cout<<"c = "<<c<<endl;
      if(c)
#ifdef USE_SMATS
        relation.add(abs(c), sign(c));
#else
        relation[abs(c)] += sign(c);
#endif
    }
#ifdef USE_SMATS
  if(relation.size()==0)
#else
  if(trivial(relation))
#endif
    {
      if (verbose) cout<<relation<<endl;
      return;
    }
  numrel++;
  if(numrel<=maxnumrel)
    {
      make_primitive(relation);
      if (verbose) cout<<relation<<endl;
      relmat.setrow(numrel,relation);
    }
  else
    cerr<<"Too many face relations (numrel="<<numrel
        <<", maxnumrel="<<maxnumrel<<")"<<endl;
}


// triangle relation for all fields
void face_relations::triangle_relation_0()
{
  if(verbose)
    {
      cout << "Face relation type 1 (triangles):\n";
    }
  vector<long> rel(3);
  vector<int> types(3,0), done(nsymb, 0);
  long j, k;
  action TiS(P1, mat22::TiS);
  action R(P1, mat22::R);
  for (k=0; k<nsymb; k++)
    if (!done[k])
      {
        for(j=0; j<3; j++)
          {
            rel[j] = (j? TiS(rel[j-1]): k);
            done[rel[j]] = 1;
            if (plusflag)
              done[R(rel[j])] = 1;
          }
        add_face_rel(rel, types);
      }
  if(verbose)
    {
      cout << "After generic triangle relation, number of relations = " << numrel <<"\n";
    }
}

// extra triangle relation for fields 1, 3
void face_relations::triangle_relation_1_3()
{
  if(verbose)
    {
      cout << "Face relation type 2 (triangles):\n";
    }
  vector<long> rel(3);
  vector<int> types(3,0), done(nsymb, 0);
  long j, k;

  Quad w(0,1);
  long field = Quad::d;
  action X = (field==1? action(P1,w,1,1,0): action(P1,1,w,w-1,0));
  assert (X.det()==(field==1? -1: 1));

  for (k=0; k<nsymb; k++)
    if (!done[k])
      {
        for(j=0; j<3; j++)
          {
            rel[j] = (j? X(rel[j-1]): k);
            done[rel[j]] = 1;
          }
        add_face_rel(rel, types);
      }
  if(verbose)
    {
      cout << "After extra triangle relation, number of relations = " << numrel <<"\n";
    }
}

// extra square relation for field 2
void face_relations::square_relation_2()
{
  if(verbose)
    {
      cout << "Face relation type 2 (squares):\n";
    }
  vector<long> rel(4);
  vector<int> types(4,0), done(nsymb, 0);
  long j, k;

  Quad w(0,1);
  action U(P1,w,1,1,0);  assert (U.det()==-1);
  action S(P1, mat22::S);
  action J(P1, mat22::J);

  for (k=0; k<nsymb; k++)
    if (!done[k])
      {
        for(j=0; j<4; j++)
          {
            rel[j] = (j? U(rel[j-1]): k);
            if(plusflag) // NB det(U)=-1
              {
                done[rel[j]] = 1;
                done[S(rel[j])] = 1;
              }
            else
              {
                if(j%2==0)
                  done[rel[j]] = 1;
              }
          }
        if (!plusflag) // since det(U)=-1
          {
            rel[1] = J(rel[1]);
            rel[3] = J(rel[3]);
          }
        add_face_rel(rel, types);
      }
  if(verbose)
    {
      cout << "After face relation type 2, number of relations = " << numrel <<"\n";
    }
}

// extra rectangle relation for field 7
void face_relations::rectangle_relation_7()
{
  if(verbose)
    {
      cout << "Face relation type 2 (rectangles):\n";
    }
  vector<long> rel(4);
  vector<int> types(4,0), done(nsymb, 0);
  long j, k;
  Quad w(0,1);

  action Y(P1,1,-w,1-w,-1);  assert (Y.det()==1);
  action USof(P1,w,-1,1,0);  assert (USof.det()==1);
  action R(P1, mat22::R);

  for (k=0; k<nsymb; k++)
    if (!done[k])
      {
        for(j=0; j<4; j+=2) // j=0,2; j+1=1,3
          {
            rel[j] = (j? Y(rel[j-2]): k);
            done[rel[j]] = 1;
            rel[j+1] = USof(rel[j]);
            if (plusflag)
              done[R(rel[j+1])]=1;
          }
        add_face_rel(rel, types);
      }
}

// extra hexagon relation for field 11
void face_relations::hexagon_relation_11()
{
  if(verbose)
    {
      cout << "Face relation type 2 (hexagons):\n";
    }
  vector<long> rel(6);
  vector<int> types(6,0), done(nsymb, 0);
  long j, k;
  Quad w(0,1);

  //  action X(P1,1,-w,1-w,-2); // as in JC thesis (order 3)
  action X(P1,-2,w,w-1,1);      // its inverse, so the hexagon edges are in the right order
  assert (X.det()==1);
  action USof(P1,w,-1,1,0);
  assert (USof.det()==1);
  action R(P1, mat22::R);

  for (k=0; k<nsymb; k++)
    if (!done[k])
      {
        for(j=0; j<6; j+=2) // j=0,2,4; j+1=1,3,5
          {
            rel[j] = (j? X(rel[j-2]): k);
            done[rel[j]] = 1;
            rel[j+1] = USof(rel[j]);
            if(plusflag)
              done[R(rel[j+1])]=1;
          }
        add_face_rel(rel, types);
      }
  if(verbose)
    {
      cout << "After face relation type 2, number of relations = " << numrel <<"\n";
    }
}

// extra triangle relation(s) for fields 19+
// Triangles {oo, w/2, (w-1)/2} {oo, w/2, (w+1)/2}

void face_relations::triangle_relation_2()
{
  int field = Quad::d;
  Quad w = Quad::w;
  long j, k, u=(field-3)/8; // u=2, 5, 8, 20 for 19,43,67,163

  action K(P1, M_alphas[1]);  assert (K.det()==1); // oo --> (w-1)/2 --> w/2 --> oo
  action N(P1, 1+w,u-w,2,-w); assert (N.det()==1); // oo --> (w+1)/2 --> w/2 --> oo

  // N is the conjugate of K by [-1,w;0,1] which maps the first
  // triangle to the second with determinant -1.  Both have order 3 so
  // cycle the edges of each triangle around.

  // All symbols are type 1, i.e. images of {w/2,oo}.

  vector<long> rel(3);
  vector<int> types(3, 1), done(nsymb, 0);
  vector<mat22> mats = {mat22::identity, K, K*K};
  check_face_rel(mats, types);

  for (k=0; k<nsymb; k++)
    if (!done[k])
      {
        for(j=0; j<3; j++)
          {
            rel[j] = (j? K(rel[j-1]): k);
            done[rel[j]] = 1;
          }
        add_face_rel(rel, types);
      }
  if (!plusflag) // there's a second triangle (image of previous under L)
    {
      std::fill(done.begin(), done.end(), 0);
      for (k=0; k<nsymb; k++)
        if (!done[k])
          {
            for(j=0; j<3; j++)
              {
                rel[j] = (j? N(rel[j-1]): k);
                done[rel[j]] = 1;
              }
            add_face_rel(rel, types);
          }
    }
  if(verbose)
    {
      cout << "After type 2 triangle relations, number of relations = " << numrel <<"\n";
    }
}

// extra triangle relation(s) for fields 43+
//
// assume  -alpha[i] = alpha[i+1] if i is odd else alpha[i-1] for i>=3

int flip(int i) {return (i<3? i: (i&1? i+1: i-1));}

// Template for all other cyclic triangle relations, given M=M_alphas[i] of order 3
//
// The triangle has vertices [alpha_i, oo, alpha_i'] with M mapping alpha_i --> oo --> alpha_i' --> alpha_i,
// and edges
// (I)_i   = {alpha_i, oo},
// (M)_i   = {M(alpha_i), M(oo)} = {oo, alpha_i'}
// (M^2)_i = {M(oo), M(alpha_i')} = {alpha_i', alpha_i}
//

void face_relations::cyclic_triangle_relation(int i)
{
  if(verbose) cout << "Applying cyclic triangle relation "<<i<<"\n";

  long j, s;
  int Ji = flip(i);
  action M(P1, M_alphas[i]);
  action J(P1, mat22::J);

  vector<mat22> mats = {mat22::identity, M, M*M};
  vector<mat22> Jmats = {mat22::identity, J*M*J, J*M*M*J};
  vector<long> rel(3);
  vector<int> types(3, i), Jtypes(3, Ji), done(nsymb, 0);

  // if(verbose) cout<<"Checking cyclic triangle relation"<<endl;
  check_face_rel(mats, types);

  // if(verbose) cout<<"Checking cyclic triangle relation after applying J"<<endl;
  check_face_rel(Jmats, Jtypes);

  for (s=0; s<nsymb; s++)
    {
      if (done[s]) continue;
      for(j=0; j<3; j++)
        {
          rel[j] = (j? M(rel[j-1]): s);
          done[rel[j]] = 1;
        }
      add_face_rel(rel, types);
      if (!plusflag)
        {
          for (j=0; j<3; j++)
            {
              rel[j] = J(rel[j]);
            }
          add_face_rel(rel, Jtypes);
        }
    }

  if(verbose) cout << "After cyclic triangle relation "<<i<<", number of relations = " << numrel <<"\n";
}

// extra square relation for field 5
void face_relations::square_relation_5()
{
  if(verbose)
    cout << "Square relation for d=5:\n";

  vector<long> rel(4);
  vector<int> types(4, -1), // type -1 means {sigmas[1],oo}
    signs = {1,-1,1,-1},
    done(nsymb, 0);

  action M(P1, M_alphas[1]);
  action Ti(P1, mat22::Tmat(-1));
  vector<mat22> mats = {mat22::identity, Ti, M, M*Ti};

  long i, j, m, k;
  for (i=0; i<nsymb; i++)
      {
        if (done[i])
          continue;
        j = Ti(i);
        m = M(i);
        k = M(j);
        done[i] = done[m] = 1; // M has order 2
        rel = {i, j, m, k};
        check_face_rel(mats, types, signs);
        add_face_rel(rel, types, signs);
      }
  if(verbose)
    {
      cout << "After square relation, number of relations = " << numrel <<"\n";
    }
}



// Template for all other triangle relations, given M_alphas[i](alphas[j]) = x + alphas[k] with x integral

// The triangle has vertices [alpha_i, oo, alpha_j] and edges
// (I)_i = {alpha_i, oo},
// (M1)_j' = M_j' * {alpha_j',oo} = {oo, alpha_j},
// (M2)_k = M_i' * T^x * {alpha_k, oo} = M_i' * {x+alpha_k, oo} = {alpha_j, alpha_i}

// The general relation for a (c:d)-symbol s has symbols s, M1(s), M2(s).

// To see its image under J, recall that if s represents the matrix
// A=[a,b;c,d] (its coset modulo Gamma_0(N)), then J(s) represents
// JAJ.  So the image of (s)_alpha = (A)_alpha = {A(alpha), A(oo)}
// under J is {JA(alpha), JA(oo)} = {JAJ(-alpha)m JAJ(oo)} =
// (JAJ)_{-alpha} = (J(s))_{-alpha}.  Hence to apply J to a relation
// we need to both apply J to the M-symbols, and also 'flip' each
// alpha to -alpha.  This works for i>=3 since then -alpha[i] =
// alpha[flip(i)].  Hence this code needs adjusting for a general
// triangle whose types are not all at least 3.

// For field 43 the only general triangle relation is {3, 7, 4} which is OK.

// For field 67, three of the five general triangles are OK but we also
// have {0,19,24} and {1,22,17} with alpha_0=0 amd alpha_1=w/2.  The
// first is OK since we define flip(0)=0, but not the second.

void face_relations::general_triangle_relation(const vector<int>& tri)
{
  int i=tri[0], j=tri[1], k=tri[2], t;
  if(verbose)
    {
      cout << "Applying triangle relation "<<tri<<"\n";
    }
  RatQuad beta = M_alphas[j].inverse()(RatQuad(1,0));
  RatQuad gamma = M_alphas[k].inverse()(RatQuad(1,0));
  RatQuad x = M_alphas[i](beta)-gamma;
  assert (x.is_integral());

  action M1(P1, M_alphas[alpha_inv[j]]);
  action M2(P1, M_alphas[alpha_inv[i]]*mat22::Tmat(x.num()));
  action J(P1, mat22::J);
  vector<mat22> mats = {mat22::identity, M1, M2};

  int jd = alpha_inv[j];
  vector<int> types = {i,jd,k};
  // if(verbose) cout<<"Checking triangle relation"<<endl;
  check_face_rel(mats, types);

  vector<mat22> Jmats = {mat22::identity, J*M1*J, J*M2*J};
  vector<int> Jtypes = {flip(i), flip(jd), flip(k)};
  Quad w = Quad::w;
  action T1(P1, (mat22::Tmat(-w)));     // needed for type 1 since -alpha[1] = alpha[1]-w
  action T2(P1, (mat22::Tmat(1-w)));   // needed for type 2 since -alpha[2] = alpha[2] + 1-w
  for (t=0; t<3; t++)
    {
      if (types[t]==1) Jmats[t] = Jmats[t]*T1;
      if (types[t]==2) Jmats[t] = Jmats[t]*T2;
    }
  // if(verbose) cout<<"Checking triangle relation after applying J and adjustment"<<endl;
  check_face_rel(Jmats, Jtypes);

  vector<long> rel(3);
  for (long s=0; s<nsymb; s++)
    {
      rel[0] = s;
      rel[1] = M1(s);
      rel[2] = M2(s);
      add_face_rel(rel, types);
      if (!plusflag)
        {
          for (t=0; t<3; t++)
            {
              rel[t] = J(rel[t]);
              if (types[t]==1) rel[t] = T1(rel[t]);
              if (types[t]==2) rel[t] = T2(rel[t]);
            }
          add_face_rel(rel, Jtypes);
        }
    }

  if(verbose)
    {
      cout << "After triangle relation "<<tri<<", number of relations = " << numrel <<"\n";
    }
}


// generic square relation
// indices i,j,k,l and x,y,z such that beta = z + M_j(x+alpha_k') =  M_i'(y+alpha_l):

// Vertices [alpha_i, oo, z+alpha_j', beta]
// Edges
// {alpha_i, oo} = (I)_i
// {oo, z+alpha_j'} = (T^z*M_j)_j
// {z+alpha_j', beta} = (T^z*M_j*T^x*M_k)_k
// {beta, alpha_i} = (M_i'*T^y)_l

void face_relations::general_square_relation(const vector<int>& squ, const vector<Quad>& xyz)
{
  int i=squ[0], j=squ[1], k=squ[2], l=squ[3], t;
  Quad x = xyz[0], y=xyz[1], z=xyz[2];
  int symmetric = (i==k)&&(j==l);
  if(verbose)
    {
      cout << "Applying ";
      if (symmetric) cout<<"symmetric ";
      cout<<"square relation "<<squ<<" (x,y,z) = "<<xyz<<"\n";
    }

  action M1(P1, mat22::Tmat(z) * M_alphas[j]);
  action M2(P1, M1 * mat22::Tmat(x) * M_alphas[k]);
  action M3(P1, M_alphas[alpha_inv[i]] * mat22::Tmat(y));
  action J(P1, mat22::J);
  vector<mat22> mats = {mat22::identity, M1, M2, M3};

  vector<int> types = squ;
  // if(verbose) cout<<"Checking square relation"<<endl;
  check_face_rel(mats, types);

  vector<mat22> Jmats = {mat22::identity, J*M1*J, J*M2*J, J*M3*J};
  vector<int> Jtypes = {flip(i), flip(j), flip(k), flip(l)};
  Quad w = Quad::w;
  action T1(P1, (mat22::Tmat(-w)));     // needed for type 1 since -alpha[1] = alpha[1]-w
  action T2(P1, (mat22::Tmat(1-w)));   // needed for type 2 since -alpha[2] = alpha[2] + 1-w
  for (t=0; t<4; t++)
    {
      if (types[t]==1) Jmats[t] = Jmats[t]*T1;
      if (types[t]==2) Jmats[t] = Jmats[t]*T2;
    }
  // if(verbose) cout<<"Checking square relation after applying J and adjustment"<<endl;
  check_face_rel(Jmats, Jtypes);

  vector<long> rel(4);
  vector<int> done(nsymb, 0);
  for (long s=0; s<nsymb; s++)
    {
      if (done[s]) continue;
      rel[0] = s;
      rel[1] = M1(s);
      rel[2] = M2(s);
      rel[3] = M3(s);
      done[s] = 1;
      if (symmetric) done[rel[2]] = 1;
      add_face_rel(rel, types);
      if (!plusflag)
        {
          for (t=0; t<4; t++)
            {
              rel[t] = J(rel[t]);
              if (types[t]==1) rel[t] = T1(rel[t]);
              if (types[t]==2) rel[t] = T2(rel[t]);
            }
          add_face_rel(rel, Jtypes);
        }
    }

  if(verbose)
    {
      cout << "After square relation "<< squ <<", number of relations = " << numrel <<"\n";
    }
}

void face_relations::solve_relations()
{
  vec npivs; // pivs is a class attribute
  long i;
  if(verbose>1)
    {
      mat M = relmat.as_mat().slice(numrel,ngens);
      cout<<"relmat = "<<M<<endl;
      cout<<"rank(relmat) = "<<relmat.rank()<<", ngens = "<<ngens<<endl;
    }
#if(USE_SMATS)
#ifdef TIME_RELATION_SOLVING
  if (verbose)
    {
      cout<<"\nStarting to solve relation matrix of size "<<numrel<<" x "<<ngens<<"\n";
    }
   timer t;
   t.start("relation solver");
#endif
   smat_elim sme(relmat);
   int d1;
   smat ker = sme.kernel(npivs, pivs), sp;
#ifdef TIME_RELATION_SOLVING
   t.stopAll();
   if (verbose)
     {
       cout<<"Finished solving relation matrix in ";
       t.showAll();
       cout<<"\n";
     }
#endif
   int ok = liftmat(ker,MODULUS,sp,d1);
   if (!ok)
     {
       if(verbose)
         cout << "failed to lift modular kernel using modulus "
              << MODULUS << endl;
#ifdef USE_CRT
       int mod2 = 1073741783; // 2^30-41
       bigint mmod = to_ZZ(MODULUS)*to_ZZ(mod2);
       if(verbose)
         cout << "repeating kernel computation, modulo " << mod2 << endl;
       smat_elim sme2(relmat,mod2);
       vec pivs2, npivs2;
       smat ker2 = sme2.kernel(npivs2,pivs2), sp2;
       ok = (pivs==pivs2);
       if (!ok)
         {
           cout<<"pivs do not agree:\npivs  = "<<pivs<<"\npivs2 = "<<pivs2<<endl;
         }
       else
         {
           if(verbose) cout << " pivs agree" << endl;
           ok = liftmats_chinese(ker,MODULUS,ker2,mod2,sp,d1);
         }
       if (ok)
         {
           if(verbose)
             cout << "success using CRT, combined modulus = "<<mmod
                  <<", denominator= " << d1 << "\n";
         }
       else
         {
           if(verbose)
             cout << "CRT combination with combined modulus "<<mmod<<" failed\n" << endl;
         }
#endif
     }
   if (ok)
     {
       hmod = 0;
       denom = d1;
     }
   else
     {
       hmod = MODULUS;
       denom = 1;
     }
   relmat=smat(0,0); // clear space
   if(verbose>1)
     {
       cout << "kernel of relmat = " << sp.as_mat() << endl;
       cout << "pivots = "<<pivs << endl;
       cout << "denom = "<<d1 << endl;
     }
   rk = sp.ncols();
   coord.init(ngens,rk); // 0'th is unused
   for(i=1; i<=ngens; i++)
     coord.setrow(i,sp.row(i).as_vec());
   // if hmod>0, coord is only defined modulo hmod
   sp=smat(0,0); // clear space
#else
  relmat = relmat.slice(numrel,ngens);
  if(verbose)
    {
      cout << "relmat = "; relmat.output_pari(cout); cout << endl;
    }
  msubspace sp = kernel(mat_m(relmat),0);
  rk = dim(sp);
  if(verbose>2) cout<<"coord = "<<basis(sp)<<endl;
  coord = basis(sp).shorten((int)1);
  pivs = pivots(sp);
  denom = I2int(denom(sp));
  relmat.init(); sp.clear();
#endif
}

#ifdef USE_CRT

// Implementations of linear algebra utility function liftmats_chinese.

//#define DEBUG_CHINESE

int liftmats_chinese(const smat& m1, scalar pr1, const smat& m2, scalar pr2,
                     smat& m, scalar& dd)
{
  long modulus=(long)pr1*(long)pr2,n,d,mij;
  long nr,nc,u,v;
  float lim=floor(sqrt(modulus/2.0));

  dd = bezout(pr1,pr2,u,v); //==1
  if (dd!=1) return 0;

  // First time through: compute CRTs, common denominator and success flag
  m = m1; // NB We assume that m1 and m2 have nonzero entries in the same places
  for(nr=0; nr<m1.nro; nr++)
    for(nc=0; nc<m1.col[nr][0]; nc++)
      {
        mij = mod(v*m1.val[nr][nc],pr1)*pr2 + mod(u*m2.val[nr][nc],pr2)*pr1;
        mij = mod(mij,modulus);
#ifdef DEBUG_CHINESE
        if (((mij-m1.val[nr][nc])%pr1)||((mij-m2.val[nr][nc])%pr2))
          {
            cout<< "bad CRT(["<<m1.val[nr][nc]<<","<<m2.val[nr][nc]<<"],["<<pr1<<","<<pr2<<"]) = "<<mij<<endl;
          }
#endif
        m.val[nr][nc] = mij;
	if (modrat(mij,modulus,lim,n,d))
          dd=lcm(d,dd);
        else
          {
#ifdef DEBUG_CHINESE
            cout<<"CRT("<<m1.val[nr][nc]<<","<<m2.val[nr][nc]<<")="<<mij<<" (mod "<<modulus<<") fails to lift (lim="<<lim<<")\n";
            cout << "Problems encountered in chinese lifting of smat modulo "<<pr1<<" and "<<pr2<< endl;
#endif
            return 0;
          }
      }
  dd=abs(dd);
#ifdef DEBUG_CHINESE
  cout << "Common denominator = " << dd << "\n";
#endif
  // Second time through: rescale
  for(nr=0; nr<m.nro; nr++)
    for(nc=0; nc<m.col[nr][0]; nc++)
      {
        m.val[nr][nc] = mod(xmodmul((dd/d),(long)m.val[nr][nc],modulus),modulus);
      }
  return 1;
}

#endif

