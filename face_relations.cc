// FILE FACE_RELATIONS.CC: Implemention of the face relations for class homspace

#include <eclib/msubspace.h>
#include <eclib/xmod.h>
#include "homspace.h"
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
  vector<RatQuad>::const_iterator alpha, beta;
  int ok=1;
  for (alpha=alphas.begin()+1, beta=betas.begin(); alpha!=alphas.end(); alpha++, beta++)
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

  switch (field)
    {
    case 1:
    case 3:
      {
        if(!plusflag)
          triangle_relation_1_3();
        break;
      }
    case 2:
      {
        square_relation_2();
        break;
      }
    case 7:
      {
        rectangle_relation_7();
        break;
      }
    case 11:
      {
        hexagon_relation_11();
        break;
      }
    case 19:
      {
        triangle_relation_2();
        square_relation_19();
        break;
      }
    case 43:
      {
        triangle_relation_2();
        triangle_relation_3();
        square_relation_43();
        break;
      }
    default:
      {
        cerr<<"homspace not implemented for field "<<field<<endl;
      }
    }
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
  symbop TS(this,1,1,-1,0); assert (TS.det()==1);
  symbop R(this,0,1,1,0);   assert (R.det()==-1);
  for (k=0; k<nsymb; k++)
    if (!done[k])
      {
        for(j=0; j<3; j++)
          {
            rel[j] = (j? TS(rel[j-1]): k);
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
  symbop S(this,0,-1,1,0); assert (S.det()==1);
  symbop J(this,fundunit,0,0,1);  assert (J.det()==fundunit);

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
  symbop R(this,0,1,1,0);      assert (R.det()==-1);

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
  symbop R(this,0,1,1,0);
  assert (R.det()==-1);

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
  vector<mat22> mats = {mat22(1,0,0,1), K, K*K};
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
void homspace::square_relation_19()
{
  int j, k;
  vector<int> rel(4), types(4), done(nsymb, 0);

  symbop K(this, M_alphas[1]);
  symbop S(this, M_alphas[0]);
  types = {0,1,0,1};
  vector<mat22> mats = {mat22(1,0,0,1), K, K*S, S};
  check_face_rel(mats, types);

  for (j=0; j<nsymb; j++)
    if (!done[j])
      {
        rel[0] = j;
        rel[3] = S(j);
        rel[1] = k = K(j);
        rel[2] = S(k); // = KS(j)
        done[j] = done[k] = 1;
        add_face_rel(rel, types);
      }
  if(verbose)
    {
      cout << "After type 0011 square relations, number of relations = " << numrel <<"\n";
    }
}

// extra triangle relation(s) for fields 43+

// triangles {oo, w/3, (w+1)/3} and (its image under J) {oo, -w/3, -(w+1)/3}

void homspace::triangle_relation_3()
{
  int field = Quad::d;
  Quad w = Quad::w;
  int k, u=(field+5)/12; // u= 4, 6, 14 for 43, 67, 163

  vector<int> rel(3), types(3);

  // Triangle {oo, w/3, (w+1)/3}

  symbop M1(this, -w,u,-3,1-w); // maps {(1-w)/3,oo} to {oo, w/3}
  assert (M1.det()==1);
  symbop M2(this, w+1,u,3,2-w); // maps {(w-1)/3,oo} to {w/3, (w+1)/3}
  assert (M2.det()==1);
  types = {5,6,7};  // {(1-w)/3, oo}, {(w-1)/3, oo}, {(1+w)/3, oo}

  vector<mat22> mats = {M1, M2, mat22(1,0,0,1)};
  check_face_rel(mats, types);

  for (k=0; k<nsymb; k++)
    {
      rel[0] = M1(k);
      rel[1] = M2(k);
      rel[2] = k;
      add_face_rel(rel, types);
    }

  if (!plusflag) // there's a second triangle (image of previous under J)
    {
      // Triangle {oo, -w/3, -(w+1)/3}

      symbop M3(this, w,u,-3,w-1);    // = J*M1*J, maps {(w-1)/3,oo} to {oo, -w/3}
      assert (M3.det()==1);
      symbop M4(this, w+1,-u,-3,2-w); // = J*M2*J, maps {(1-w)/3,oo} to {-w/3, -(w+1)/3}
      assert (M4.det()==1);
      types = {6,5,8};  // {(w-1)/3, oo}, {(1-w)/3, oo}, {-(1+w)/3, oo}

      for (k=0; k<nsymb; k++)
        {
          rel[0] = M3(k);
          rel[1] = M4(k);
          rel[2] = k;
          add_face_rel(rel, types);
        }
    }
  if(verbose)
    {
      cout << "After type 3 triangle relations, number of relations = " << numrel <<"\n";
    }
}

// extra square relations for field 43
void homspace::square_relation_43()
{
  //  int field = Quad::d;
  Quad w = Quad::w;
  int j, k;
  vector<int> rel(4), types(4), done(nsymb, 0);

  // (1)
  // vertices {0, oo, w/3, 3w/11}
  // edge types 0 {0,oo}, 5 {(1-w)/3,oo}, 1 {w/2,oo}, 5 under I, V, N4, N4*V
  // no symmetries

  symbop N4L(this, -3,2*w,w-1,7); assert (N4L.det()==1);
  symbop V(this, w,-4,3,w-1);     assert (V.det()==1);
  symbop S(this, 0,-1,1,0);       assert (S.det()==1);
  types = {0,5,1,6};

  vector<mat22> mats = {mat22(1,0,0,1), V, N4L, S};
  check_face_rel(mats, types);

  for (j=0; j<nsymb; j++)
    {
      rel[0] = j;
      rel[1] = V(j);
      rel[2] = N4L(j);
      rel[3] = S(j);
      add_face_rel(rel, types);
    }

  // (2)
  // vertices {w/2,oo,(w+1)/3,(5*w+3)/13}
  // edge types 1 {w/2,oo}, 7 {(w+1)/3,oo}, 1, 7 under I, U, N5, N5*U

  symbop N5(this, w-4,w+5,w+1,4-w);   assert (N5.det()==1);
  symbop U(this, w+1,3-w,3,-1-w);     assert (U.det()==1);
  types = {1,7,1,7};
  mats = {mat22(1,0,0,1), U, N5, N5*U};
  check_face_rel(mats, types);

  for (j=0; j<nsymb; j++)
    {
      rel[0] = j;
      rel[1] = U(j);
      rel[2] = k = N5(j);
      rel[3] = U(k);
      add_face_rel(rel, types);
    }

  // (3) image of (2) under L=[-1,0;0,1] with det(L)=-1, only needed if !plusflag

  // vertices {w/2, oo, (2*w-1)/3, (8*w-3)/13}
  // edge types 1 {w/2,oo}, 8 {-(w+1)/3,oo}, 1, 7 under I, W1, W2, W2*W1

  if (!plusflag)
    {
      symbop W1(this, 2*w-1,w-8,3,w+1);   assert (W1.det()==1);
      symbop W2(this, w-7,4*w+5,w+1,7-w); assert (W2.det()==1);
      types = {1,8,1,8};
      mats = {mat22(1,0,0,1), W1, W2, W2*W1};
      check_face_rel(mats, types);

      for (j=0; j<nsymb; j++)
        {
          rel[0] = j;
          rel[1] = W1(j);
          rel[2] = k = W2(j);
          rel[3] = W1(k);
          add_face_rel(rel, types);
        }
    }

  if(verbose)
    {
      cout << "After type extra square relations, number of relations = " << numrel <<"\n";
    }
}

