// FILE FACE_RELATIONS.CC: Implemention of the face relations for class homspace

#include <eclib/msubspace.h>
#include <eclib/xmod.h>
#include "homspace.h"
#include "geometry.h"
#include <assert.h>

// Each face relation is a sum of edges (M)_alpha = {M(alpha}, M(oo)}
// for M in the list mats and alpha=alphas[i] for i in the list types.
// Here we check that such a relation holds identically in H_3 (not
// just modulo the congruence subgroup!)

int check_face_rel(const vector<mat22>& mats, const vector<int>& types)
{
  vector<mat22>::const_iterator mi;
  vector<int>::const_iterator ti;
  vector<RatQuad> alphas, betas;
  for (mi=mats.begin(), ti=types.begin(); ti!=types.end(); mi++, ti++)
    {
      mat22 M = *mi;
      mat22 M_alpha = M_alphas[*ti];
      alphas.push_back(M(RatQuad(-M_alpha.entry(1,1), M_alpha.entry(1,0))));
      betas.push_back(M(RatQuad(1,0)));
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
    }
  // else
  //   {
  //     cout<<"Good face relation "<<endl;
  //   }
  return ok;
}

// In add_face_rel(rel, types, check):
//
//   rel is a list of (positive) (c:d)-symbol numbers i
//   types is a list of symbol types
//
//   such that the corresponding (symbol,type) pairs add to 0 in homology.  We use
//   the map i -> j=coordindex[i] to convert this to a vector of
//   dimension ngens, and append that as a new row to the relation
//   matrix relmat.
//

void homspace::add_face_rel(const vector<int>& rel, const vector<int>& types)
{
  vector<int>::const_iterator r, t;
  long c;
  if (verbose)
    {
      cout<<"Relation: ";
      for (r = rel.begin(), t = types.begin(); r!=rel.end(); r++, t++)
        cout<<(*r)<<"_"<<(*t)<<" "<<symbol(*r)<<" ";
      cout <<" --> ";
    }
#ifdef USE_SMATS
  svec relation(ngens);
#else
  vec relation(ngens);
#endif
  for (r = rel.begin(), t = types.begin(); r!=rel.end(); r++, t++)
    {
      c = coordindex[*r+nsymb*(*t)];
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

//
// face relations
//

void homspace::face_relations()
{
  long field = Quad::d;

  maxnumrel=2*nsymbx;
  numrel = 0;

#if(USE_SMATS)
  relmat = smat(maxnumrel,ngens);
#else
  relmat.init(maxnumrel,ngens);
#endif
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

  // Now field = 19, 43, 67 or 163

  // additional triangle relations
  triangle_relation_2();
  for (vector<vector<int>>::const_iterator tri = triangles.begin(); tri!=triangles.end(); tri++)
    {
      general_triangle_relation(*tri);
    }
  for (vector<pair<vector<int>, vector<Quad>> >::const_iterator S = squares.begin(); S!=squares.end(); S++)
    {
      general_square_relation(S->first, S->second);
    }
  if(field==163)
    cerr<<"homspace not fully implemented for field "<<field<<endl;
}

// triangle relation for all fields
void homspace::triangle_relation_0()
{
  if(verbose)
    {
      cout << "Face relation type 1 (triangles):\n";
    }
  vector<int> rel(3), types(3,0), done(nsymb, 0);
  long j, k;
  symbop TiS(this, mat22::TiS);
  symbop R(this, mat22::R);
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
void homspace::triangle_relation_1_3()
{
  if(verbose)
    {
      cout << "Face relation type 2 (triangles):\n";
    }
  vector<int> rel(3), types(3,0), done(nsymb, 0);
  long j, k;

  Quad w(0,1);
  long field = Quad::d;
  symbop X = (field==1? symbop(this,w,1,1,0): symbop(this,1,w,w-1,0));
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
void homspace::square_relation_2()
{
  if(verbose)
    {
      cout << "Face relation type 2 (squares):\n";
    }
  vector<int> rel(4), types(4,0), done(nsymb, 0);
  long j, k;

  Quad w(0,1);
  symbop U(this,w,1,1,0);  assert (U.det()==-1);
  symbop S(this, mat22::S);
  symbop J(this, mat22::J);

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
void homspace::rectangle_relation_7()
{
  if(verbose)
    {
      cout << "Face relation type 2 (rectangles):\n";
    }
  vector<int> rel(4), types(4,0), done(nsymb, 0);
  long j, k;
  Quad w(0,1);

  symbop Y(this,1,-w,1-w,-1);  assert (Y.det()==1);
  symbop USof(this,w,-1,1,0);  assert (USof.det()==1);
  symbop R(this, mat22::R);

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
void homspace::hexagon_relation_11()
{
  if(verbose)
    {
      cout << "Face relation type 2 (hexagons):\n";
    }
  vector<int> rel(6), types(6,0), done(nsymb, 0);
  long j, k;
  Quad w(0,1);

  //  symbop X(this,1,-w,1-w,-2); // as in JC thesis (order 3)
  symbop X(this,-2,w,w-1,1);      // its inverse, so the hexagon edges are in the right order
  assert (X.det()==1);
  symbop USof(this,w,-1,1,0);
  assert (USof.det()==1);
  symbop R(this, mat22::R);

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

void homspace::triangle_relation_2()
{
  int field = Quad::d;
  Quad w = Quad::w;
  int j, k, u=(field-3)/8; // u=2, 5, 8, 20 for 19,43,67,163

  symbop K(this, M_alphas[1]);  assert (K.det()==1); // oo --> (w-1)/2 --> w/2 --> oo
  symbop N(this, 1+w,u-w,2,-w); assert (N.det()==1); // oo --> (w+1)/2 --> w/2 --> oo

  // N is the conjugate of K by [-1,w;0,1] which maps the first
  // triangle to the second with determinant -1.  Both have order 3 so
  // cycle the edges of each triangle around.

  // All symbols are type 1, i.e. images of {w/2,oo}.

  vector<int> rel(3), types(3, 1), done(nsymb, 0);
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

// extra square relation for field 19
// void homspace::square_relation_19()
// {
//   int j, k;
//   vector<int> rel(4), types(4), done(nsymb, 0);

//   symbop K(this, M_alphas[1]);
//   symbop S(this, M_alphas[0]);
//   types = {0,1,0,1};
//   vector<mat22> mats = {mat22::identity, K, K*S, S};
//   check_face_rel(mats, types);

//   for (j=0; j<nsymb; j++)
//     if (!done[j])
//       {
//         rel[0] = j;
//         rel[3] = S(j);
//         rel[1] = k = K(j);
//         rel[2] = S(k); // = KS(j)
//         done[j] = done[k] = 1;
//         add_face_rel(rel, types);
//       }
//   if(verbose)
//     {
//       cout << "After type 0011 square relations, number of relations = " << numrel <<"\n";
//     }
// }

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

void homspace::cyclic_triangle_relation(int i)
{
  if(verbose)
    {
      cout << "Applying cyclic triangle relation "<<i<<"\n";
    }
  int j,s;
  symbop M(this, M_alphas[i]);
  symbop J(this, mat22::J);

  vector<mat22> mats = {mat22::identity, M, M*M};
  vector<int> types = {i,i,i};
  // if(verbose) cout<<"Checking cyclic triangle relation"<<endl;
  check_face_rel(mats, types);

  vector<mat22> Jmats = {mat22::identity, J*M*J, J*M*M*J};
  int Ji = flip(i);
  vector<int> Jtypes = {Ji, Ji, Ji};
  // if(verbose) cout<<"Checking cyclic triangle relation after applying J"<<endl;
  check_face_rel(Jmats, Jtypes);

  vector<int> rel(3), done(3, 0);
  for (s=0; s<nsymb; s++)
    {
      if (!done[s])
        {
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
    }

  if(verbose)
    {
      cout << "After cyclic triangle relation "<<i<<", number of relations = " << numrel <<"\n";
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

void homspace::general_triangle_relation(const vector<int>& tri)
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

  symbop M1(this, M_alphas[alpha_inv[j]]);
  symbop M2(this, M_alphas[alpha_inv[i]]*mat22::Tmat(num(x)));
  symbop J(this, mat22::J);
  vector<mat22> mats = {mat22::identity, M1, M2};

  int jd = alpha_inv[j];
  vector<int> types = {i,jd,k};
  // if(verbose) cout<<"Checking triangle relation"<<endl;
  check_face_rel(mats, types);

  vector<mat22> Jmats = {mat22::identity, J*M1*J, J*M2*J};
  vector<int> Jtypes = {flip(i), flip(jd), flip(k)};
  Quad w = Quad::w;
  symbop T1(this, (mat22::Tmat(-w)));     // needed for type 1 since -alpha[1] = alpha[1]-w
  symbop T2(this, (mat22::Tmat(-1-w)));   // needed for type 2 since -alpha[2] = alpha[2]-1-w
  for (t=0; t<3; t++)
    {
      if (types[t]==1) Jmats[t] = Jmats[t]*T1;
      if (types[t]==2) Jmats[t] = Jmats[t]*T2;
    }
  // if(verbose) cout<<"Checking triangle relation after applying J and adjustment"<<endl;
  check_face_rel(Jmats, Jtypes);

  vector<int> rel(3);
  for (int s=0; s<nsymb; s++)
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

// extra square relations for field 43
// void homspace::square_relation_43()
// {
//   int verb=verbose; //verbose=1;
//   if(verbose) cout<<"Square relation 1"<<endl;
//   //  int field = Quad::d;
//   Quad w = Quad::w;
//   int j, k, t;
//   vector<int> rel(4), types(4), done(nsymb, 0);

//   // (1)
//   // vertices {0, oo, w/3, 3w/11}
//   // edge types 0 {0,oo}, 5 {(1-w)/3,oo}, 1 {w/2,oo}, 6 {(w-1)/3,oo} under I, M5, M, S
//   // -->  {0,oo}, {oo,w/3}, {w/3,3w/11}, {3w/11, 0}
//   // no symmetries

//   symbop M(this, M_alphas[5]*mat22::Tmat(1-w)*M_alphas[1]);
//   //  symbop M(this, -3,2*w,w-1,7); assert (M.det()==1);

//   symbop M5(this, M_alphas[5]);
//   symbop S(this, mat22::S);
//   symbop J(this, mat22::J);
//   types = {0,5,1,6};

//   vector<mat22> mats = {mat22::identity, M5, M, S};
//   check_face_rel(mats, types);

//   symbop T1(this, (mat22::Tmat(-w)));     // needed for type 1 since -alpha[1] = alpha[1]-w
//   vector<mat22> Jmats = {mat22::identity, J*M5*J, J*M*J*T1, J*S*J};
//   vector<int> Jtypes = {0, 6, 1, 5};
//   check_face_rel(Jmats, Jtypes);

//   for (j=0; j<nsymb; j++)
//     {
//       rel[0] = j;
//       rel[1] = M5(j);
//       rel[2] = M(j);
//       rel[3] = S(j);
//       add_face_rel(rel, types);
//       if (!plusflag)
//         {
//           for (t=0; t<4; t++) rel[t] = J(rel[t]);
//           rel[2] = T1(rel[2]);
//           add_face_rel(rel, Jtypes);
//         }
//     }

//   // (2)
//   // vertices {w/2,oo,(w+1)/3,(5*w+3)/13}
//   // edge types 1 {w/2,oo}, 7 {(w+1)/3,oo}, 1, 7 under I, M7, N, N*M7

//   if(verbose) cout<<"Square relation 2"<<endl;
//   symbop M7(this, M_alphas[7]);
//   //  symbop N(this, w-4,w+5,w+1,4-w);   assert (N.det()==1);
//   symbop N(this, M_alphas[7]*mat22::Tmat(1)*M_alphas[1]);
//   types = {1,7,1,7};
//   mats = {mat22::identity, M7, N, N*M7};
//   check_face_rel(mats, types);
//   Jmats = {T1, J*M7*J, J*N*J*T1, J*N*M7*J};
//   Jtypes = {1, 8, 1, 8};
//   for (j=0; j<nsymb; j++)
//     {
//       if (done[j]) continue;
//       rel[0] = j;
//       rel[1] = M7(j);
//       rel[2] = k = N(j);
//       rel[3] = M7(k);
//       done[j] = done[k] = 1;
//       add_face_rel(rel, types);
//       if (!plusflag)
//         {
//           for (int t=0; t<4; t++) rel[t] = J(rel[t]);
//           rel[0] = T1(rel[0]);
//           rel[2] = T1(rel[2]);
//           add_face_rel(rel, Jtypes);
//         }
//     }

//   if(verbose)
//     {
//       cout << "After extra square relations, number of relations = " << numrel <<"\n";
//     }
//   verbose=verb;
// }

// generic square relation
// indices i,j,k,l and x,y such that beta = M_j(x+alpha_k') =  M_i'(y+alpha_l):

// Vertices [alpha_i, oo, alpha_j', beta]
// Edges
// {alpha_i, oo} = (I)_i
// {oo, alpha_j'} = (M_j)_j
// {alpha_j', beta} = (M_j*T^x*M_k)_k
// {beta, alpha_i} = (M_i'*T^y)_l

void homspace::general_square_relation(const vector<int>& squ, const vector<Quad>& xy)
{
  int i=squ[0], j=squ[1], k=squ[2], l=squ[3], t;
  int symmetric = (i==k)&&(j==l);
  if(verbose)
    {
      cout << "Applying ";
      if (symmetric) cout<<"symmetric ";
      cout<<"square relation "<<squ<<"\n";
    }

  symbop M1(this, M_alphas[j]);
  symbop M2(this, M_alphas[j] * mat22::Tmat(xy[0]) * M_alphas[k]);
  symbop M3(this, M_alphas[alpha_inv[i]] * mat22::Tmat(xy[1]));
  symbop J(this, mat22::J);
  vector<mat22> mats = {mat22::identity, M1, M2, M3};

  vector<int> types = squ;
  // if(verbose) cout<<"Checking square relation"<<endl;
  check_face_rel(mats, types);

  vector<mat22> Jmats = {mat22::identity, J*M1*J, J*M2*J, J*M3*J};
  vector<int> Jtypes = {flip(i), flip(j), flip(k), flip(l)};
  Quad w = Quad::w;
  symbop T1(this, (mat22::Tmat(-w)));     // needed for type 1 since -alpha[1] = alpha[1]-w
  symbop T2(this, (mat22::Tmat(-1-w)));   // needed for type 2 since -alpha[2] = alpha[2]-1-w
  for (t=0; t<4; t++)
    {
      if (types[t]==1) Jmats[t] = Jmats[t]*T1;
      if (types[t]==2) Jmats[t] = Jmats[t]*T2;
    }
  // if(verbose) cout<<"Checking square relation after applying J and adjustment"<<endl;
  check_face_rel(Jmats, Jtypes);

  vector<int> rel(4), done(nsymb, 0);
  for (int s=0; s<nsymb; s++)
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
