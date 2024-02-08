// FILE FACE_RELATIONS.CC: Implemention of the face relations for class homspace

//#define TIME_RELATION_SOLVING
//#define USE_CRT // if using smats  mod MODULUS, try CRT-ing with another prime
                // NB this is experimental only


#include "mat22.h"
#include "ratquads.h"
#include "homspace.h"
#include "face_relations.h"
#include <assert.h>

#ifdef USE_CRT
int liftmats_chinese(const smat& m1, scalar pr1, const smat& m2, scalar pr2, smat& m, scalar& dd);
#endif

// Each relation is a signed sum of edges (M)_alpha = {M(alpha},
// M(oo)} for M in the list mats and alpha=alphas[t] (when t>=0) or
// sigmas[-t] (when t<0), for t in the list types.  Here we check that such a
// relation holds identically in H_3 (not just modulo the congruence
// subgroup!)

// Special case: all signs +1
int check_rel(const vector<mat22>& mats, const vector<int>& types)
{
  vector<int> signs(mats.size(), 1);
  return check_rel(mats, types, signs);
}

// General case:

//#define DEBUG_FACE_RELATION

int check_rel(const vector<mat22>& mats, const vector<int>& types, const vector<int>& signs)
{
#ifdef DEBUG_FACE_RELATION
  int n = mats.size();
  cout<<"  - Checking "<<(n==2? "edge": "face")<<" relation...\n";
  cout<<"    mats: "<<mats<<endl;
  cout<<"    types: "<<types<<endl;
  cout<<"    signs: "<<signs<<endl;
#endif
  auto mi = mats.begin();
  auto ti = types.begin(), si = signs.begin();
  vector<RatQuad> as, bs;
  while (ti!=types.end())
    {
      mat22 M = *mi++;
#ifdef DEBUG_FACE_RELATION
      cout<<"    M = "<<M<<" maps {"<<base_point(*ti)<<",oo} to ";
#endif
      RatQuad a = M(base_point(*ti++));
      RatQuad b = M.image_oo();
#ifdef DEBUG_FACE_RELATION
      cout<<"{"<<a<<","<<b<<"}"<<endl;
#endif
      if (*si>0) // use {a,b} when sign is +1
        {
          as.push_back(a);
          bs.push_back(b);
        }
      else // use {b,a} when sign is -1
        {
          as.push_back(b);
          bs.push_back(a);
        }
      si++;
    }

  auto ai = as.begin()+1, bi = bs.begin();
  int ok=1;
  while ( bi!=bs.end() &&ok)
    {
      RatQuad next_alpha = (ai==as.end()? as[0]: *ai++);
      ok = ok && (*bi++==next_alpha);
    }
  if (!ok)
    {
      int nsides = mats.size();
      cout<<"\n*************Bad "<< (nsides==2? "edge": "face") << " relation!\n";
      cout<<"alphas: "<<as<<endl;
      cout<<"betas:  "<<bs<<endl;
      exit(1);
    }
#ifdef DEBUG_FACE_RELATION
  else
    {
      cout<<"  - Good "<< (n==2? "edge": "face") << " relation!\n";
    }
#endif
  return ok;
}

//
// face relations
//

face_relations::face_relations(edge_relations* er, int plus, int verb, long ch)
  :ER(er), plusflag(plus), verbose(verb), characteristic(ch)
{
  P1 = ER->P1;
  nsymb = P1->size();
  maxnumrel = 2*(plusflag?1:2)*nsymb;
  if (verbose)
    cout<<"initial bound on #relations = "<<maxnumrel<<endl;

  if (!Quad::is_Euclidean)
    {
      if (verbose)
        {
          cout<<cyclic_triangles.size()<<" extra cyclic triangle relation types"<<endl;
          cout<<aaa_triangles.size()<<" extra aaa triangle relation types"<<endl;
          cout<<aas_triangles.size()<<" aas triangle relation types"<<endl;
          cout<<squares.size()<<" square relation types"<<endl;
          cout<<hexagons.size()<<" hexagon relation types"<<endl;
        }
      maxnumrel += (plusflag?1:2)*nsymb*(cyclic_triangles.size() + aaa_triangles.size() + aas_triangles.size() + squares.size() + hexagons.size());
      if (verbose)
        cout<<"final bound on #relations = "<<maxnumrel<<endl;
    }

  ngens = ER->ngens;
  numrel = 0;
  hmod = 0;

#ifdef USE_SMATS
  relmat = smat(maxnumrel,ngens);
  relmat_rowdata = vector<int>(10,0);
  maxrowsize = 0;
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
          if (hmod && characteristic==0)
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

  // additional triangle and square relations

  if (Quad::class_number==1)
    triangle_relation_2();

  if (!cyclic_triangles.empty())
    {
      if(verbose)
        cout<<"\nApplying "<<cyclic_triangles.size()<<" cyclic triangle relations"<<endl;
      for ( const auto& T : cyclic_triangles)
        cyclic_triangle_relation(T);
    }

  if (!aaa_triangles.empty())
    {
      if(verbose)
        cout<<"\nApplying "<<aaa_triangles.size()<<" general aaa-triangle relations"<<endl;
      for ( const auto& T : aaa_triangles)
        aaa_triangle_relation(T);
    }

  if (!squares.empty())
    {
      if(verbose)
        cout<<"\nApplying "<<squares.size()<<" general square relations"<<endl;
      for ( const auto& S : squares)
        general_square_relation(S);
    }

  if (!hexagons.empty())
    {
      if(verbose)
        cout<<"\nApplying "<<hexagons.size()<<" general hexgaon relations"<<endl;
      for ( const auto& H : hexagons)
        general_hexagon_relation(H);
    }

  if (Quad::class_number==1)
    return;

  if (!aas_triangles.empty())
    {
      if(verbose)
        cout<<"\nApplying "<<aas_triangles.size()<<" general aas-triangle relations"<<endl;
      for ( const auto& T : aas_triangles)
        aas_triangle_relation(T);
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
  if (verbose)
    {
      cout<<"Relation: ";
      Quad c, d;
      auto r = rel.begin();
      auto t = types.begin(), s=signs.begin();
      while (r!=rel.end())
        {
          P1->make_symb(*r, c, d);
          cout<< ((*s)>0? " +": " -");
          //cout<<"("<<c<<":"<<d<<")";
          cout<<"["<<(*r)<<";"<<(*t)<<"]";
          ++r, ++t, ++s;
        }
      cout <<" --> ";
    }
#ifdef USE_SMATS
  svec relation(ngens);
#else
  vec relation(ngens);
#endif
  auto r = rel.begin();
  auto t = types.begin(), s=signs.begin();
  while (r!=rel.end())
    {
      //      cout<<"Looking up edge coord of symbol "<<(*r)<<", type "<<(*t)<<"...";
      long c = (*s) * ER->coords(*r, *t);
      //      cout<<"c = "<<c<<endl;
      if(c)
        {
#ifdef USE_SMATS
          relation.add(abs(c), sign(c));
#else
          relation[abs(c)] += sign(c);
#endif
        }
      ++r, ++t, ++s;
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
  int n = relation.size();
  relmat_rowdata[n] +=1;
  if (n>maxrowsize) maxrowsize=n;
  numrel++;
  if(numrel<=maxnumrel)
    {
      if (characteristic==0)
        make_primitive(relation);
      else
        relation.reduce_mod_p(characteristic);
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
      cout << "Face relation type 1 (universal triangle):\n";
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
      cout << "After universal triangle relation, number of relations = " << numrel <<"\n";
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

  Quad w=Quad::w;
  long field = Quad::d;
  action X = (field==1? action(P1,w,ONE,ONE,ZERO): action(P1,ONE,w,w-ONE,ZERO));
  assert (X.det()==(field==1? -ONE: ONE));

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

  Quad w=Quad::w;
  action U(P1,w,ONE,ONE,ZERO);  assert (U.det()==-ONE);
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
  Quad w=Quad::w;

  action Y(P1,ONE,-w,ONE-w,-ONE); assert (Y.is_unimodular());
  action US(P1,w,-ONE,ONE,ZERO);   assert (US.is_unimodular());
  action R(P1, mat22::R);

  for (k=0; k<nsymb; k++)
    if (!done[k])
      {
        for(j=0; j<4; j+=2) // j=0,2; j+1=1,3
          {
            rel[j] = (j? Y(rel[j-2]): k);
            done[rel[j]] = 1;
            rel[j+1] = US(rel[j]);
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
  Quad w=Quad::w;

  //  action X(P1,ONE,-w,ONE-w,-2); // as in JC thesis (order 3)
  action X(P1,-TWO,w,w-ONE,ONE);      // its inverse, so the hexagon edges are in the right order
  assert (X.is_unimodular());
  action US(P1,w,-ONE,ONE,ZERO);
  assert (US.is_unimodular());
  action R(P1, mat22::R);

  for (k=0; k<nsymb; k++)
    if (!done[k])
      {
        for(j=0; j<6; j+=2) // j=0,2,4; j+1=1,3,5
          {
            rel[j] = (j? X(rel[j-2]): k);
            done[rel[j]] = 1;
            rel[j+1] = US(rel[j]);
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

// extra triangle relation(s)

// Triangles {oo, w/2, (w-1)/2} {oo, w/2, (w+1)/2}  for non-Euclidean fields of class number 1

void face_relations::triangle_relation_2()
{
  long field = Quad::d;
  Quad w=Quad::w;
  long j, k;
  Quad u(INT(field-3)/8); // u=2, 5, 8, 20 for 19,43,67,163

  action K(P1, M_alphas[1]);  assert (K.is_unimodular()); // oo --> (w-1)/2 --> w/2 --> oo
  action N(P1, ONE+w,u-w,TWO,-w); assert (N.is_unimodular()); // oo --> (w+1)/2 --> w/2 --> oo

  // N is the conjugate of K by [-1,w;0,1] which maps the first
  // triangle to the second with determinant -1.  Both have order 3 so
  // cycle the edges of each triangle around.

  // All symbols are type 1, i.e. images of {w/2,oo}.

  vector<long> rel(3);
  vector<int> types(3, 1), done(nsymb, 0);
  vector<mat22> mats = {mat22::identity, K, K*K};
  assert (check_rel(mats, types));

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

// extra relation(s) (not needed for Euclidean fields)

void face_relations::general_relation(const vector<action>& Mops,
                                      const vector<int>& types,
                                      const vector<int>& signs,
                                      int symmetry, int check)
{
  int len = types.size();
  vector<mat22> Mats(len);
  vector<int> sym(len, 0);
  for (int s=0; s<len; s++)
    {
      if (symmetry) sym[s] = (s%symmetry==0);
      Mats[s] = Mops[s];
    }
  mat22 J = mat22::J;
  vector<action> Jops(len, action(P1, J));
  vector<int> Jtypes(len);   // Jops, Jtypes used in applying J to a relation

  if (check)
    assert(check_rel(Mats, types, signs));

  // Adjustments will be needed on applying J when one of the alphas[t] has denominator 2.
  if (!plusflag)
    {
      vector<mat22> Jmats;  // Jmats only used in checking validity of relation
      for (int s=0; s<len; s++)
        {
          if (check)
            Jmats.push_back(J*Mats[s]*J);
          int t = types[s];
          RatQuad a;
          Quad x;
          if (t>=0)
            {
              Jtypes[s] = alpha_flip[t];
              a = TWO*alphas[t];
            }
          else
            {
              Jtypes[s] = -sigma_flip[-t];
              a = TWO*sigmas[-t];
            }
          if (a.is_integral(x))
            {
              mat22 T = mat22::Tmat(-x);
              Jops[s] = action(P1, J*T);
              if (check)
                Jmats[s] *= T;
            }
        }
      if (check)
        assert(check_rel(Jmats, Jtypes, signs));
    }

  vector<int> done(nsymb, 0);
  for (long j=0; j<nsymb; j++)
    {
      if (done[j]) continue;
      vector<long> rel(len);
      for (int s=0; s<len; s++)
        {
          long k = Mops[s](j);  // NB first matrix is NOT I for hexagons
          rel[s] = k;
          if (sym[s]) done[k]=1;
        }
      add_face_rel(rel, types, signs);
      if (!plusflag)
        {
          for (int s=0; s<len; s++)
            {
              rel[s] = Jops[s](rel[s]);
            }
          add_face_rel(rel, Jtypes, signs);
        }
    }
}

// Template for all other cyclic triangle relations, given M=M_alphas[i] of order 3
//
// The triangle has vertices [alpha_i, oo, alpha_i'] with M mapping alpha_i --> oo --> alpha_i' --> alpha_i,
// and edges
// (I)_i   = {alpha_i, oo},
// (M)_i   = {M(alpha_i), M(oo)} = {oo, alpha_i'}
// (M^2)_i = {M(oo), M(alpha_i')} = {alpha_i', alpha_i}
//

void face_relations::cyclic_triangle_relation(int i, int check)
{
  if(verbose) cout << "Applying cyclic triangle relation "<<i<<"\n";

  mat22 M = M_alphas[i];
  vector<int> types(3, i), signs(3, 1);
  vector<action> Mops = {action(P1, mat22::identity),
                         action(P1, M),
                         action(P1, M*M)};

  general_relation(Mops, types, signs, 1, check);

  if(verbose) cout << "After cyclic triangle relation "<<i<<", number of relations = " << numrel <<"\n";
}

// Template for other aaa-triangle relations, given M_alphas[i](alphas[j]+u) = x + alphas[k] with x integral

// The triangle has vertices [alpha_i, oo, alpha_j+u] and edges
// (I)_i = {alpha_i, oo},
// (M1)_j' = T^u * M_j' * {alpha_j',oo} = {oo, alpha_j +u},
// (M2)_k = M_i' * T^x * {alpha_k, oo} = M_i' * {x+alpha_k, oo} = {alpha_j +u, alpha_i}

// The general relation for a (c:d)-symbol s has symbols s, M1(s), M2(s).

// To see its image under J, recall that if s represents the matrix
// A=[a,b;c,d] (its coset modulo Gamma_0(N)), then J(s) represents
// JAJ.  So the image of (s)_alpha = (A)_alpha = {A(alpha), A(oo)}
// under J is {JA(alpha), JA(oo)} = {JAJ(-alpha)m JAJ(oo)} =
// (JAJ)_{-alpha} = (J(s))_{-alpha}.  Hence to apply J to a relation
// we need to both apply J to the M-symbols, and also 'flip' each
// alpha to -alpha.  This is easy when denom(alpha)!=2 since then
// -alpha[i] = alpha[alpha_flip[i]], but needs adjusting for a general
// triangle involving types with denominator 2.

// For field 43 the only general triangle relation is {3, 7, 4} which is OK.

// For field 67, three of the five general triangles are OK but we
// also have {0,19,24} and {1,22,17} with alpha_0=0 amd alpha_1=w/2.
// The first is OK since alpha_flip[0]=0, but not the second.

void face_relations::aaa_triangle_relation(const TRIANGLE& tri, int check)
{
  vector<int> T = tri.first;
  Quad u = tri.second[0];
  int i=T[0], j=T[1], k=T[2];
  if(verbose)
    cout << "Applying aaa-triangle relation ["<<T<<"; "<<u<<"]\n";

  RatQuad beta = M_alphas[j].preimage_oo();
  RatQuad gamma = M_alphas[k].preimage_oo();
  RatQuad x = M_alphas[i](beta+u) - gamma;
  Quad xnum;
  x.is_integral(xnum);

  vector<int> types = {i,alpha_inv[j],k}, signs(3, 1);
  vector<action> Mops = {action(P1, mat22::identity),
                         action(P1, mat22::Tmat(u) * M_alphas[alpha_inv[j]]),
                         action(P1, M_alphas[alpha_inv[i]]*mat22::Tmat(xnum))};

  general_relation(Mops, types, signs, 0, check);

  if(verbose)
    cout << "After aaa-triangle relation ["<<T<<"; "<<u<<"], number of relations = " << numrel <<"\n\n";
}

void face_relations::aas_triangle_relation(const TRIANGLE& tri, int check)
{
  vector<int> T = tri.first;
  Quad u = tri.second[0];
  int i=T[0], j=T[1], k=T[2];
  if(verbose)
    cout << "Applying aas-triangle relation ["<<T<<"; "<<u<<"]\n";

  RatQuad x = M_alphas[i](sigmas[j]+u) - sigmas[k];
  Quad xnum;
  x.is_integral(xnum);

  vector<int> types = {i,-j,-k}, signs = {+1,-1,+1};
  vector<action> Mops = {action(P1, mat22::identity),
                         action(P1, mat22::Tmat(u)),
                         action(P1, M_alphas[alpha_inv[i]]*mat22::Tmat(xnum))};

  general_relation(Mops, types, signs, 0, check);

  if(verbose)
    cout << "After aas-triangle relation ["<<T<<"; "<<u<<"], number of relations = " << numrel <<"\n\n";
}

// generic square relation
// indices i,j,k,l and x,y,z such that beta = z + M_j(x+alpha_k') =  M_i'(y+alpha_l):

// Vertices [alpha_i, oo, z+alpha_j', beta]
// Edges
// {alpha_i, oo} = (I)_i
// {oo, z+alpha_j'} = (T^z*M_j)_j
// {z+alpha_j', beta} = (T^z*M_j*T^x*M_k)_k
// {beta, alpha_i} = (M_i'*T^y)_l

void face_relations::general_square_relation(const POLYGON& S, int check)
{
  const vector<int>& squ = S.first;
  const vector<Quad>& xyz = S.second;
  int i=squ[0], j=squ[1], k=squ[2], l=squ[3];
  Quad x = xyz[0], y=xyz[1], z=xyz[2];
  int symmetry = (((i==k)&&(j==l))? 2: 0);
  if (Quad::d==5) symmetry=2; // not yet automatic
  if(verbose)
    {
      cout << "Applying ";
      if (symmetry) cout<<"symmetric ";
      cout<<"square relation "<<squ<<" (x,y,z) = "<<xyz<<"\n";
    }

  vector<int> types = squ, signs(4, 1);
  mat22 M1 = mat22::Tmat(z) * M_alphas[j];
  mat22 M2 = M1 * mat22::Tmat(x) * M_alphas[k];
  mat22 M3 = M_alphas[alpha_inv[i]] * mat22::Tmat(y);
  vector<action> Mops = {action(P1, mat22::identity),
                         action(P1, M1),
                         action(P1, M2),
                         action(P1, M3)};

  general_relation(Mops, types, signs, symmetry, check);

  if(verbose)
    cout << "After square relation "<< squ <<", number of relations = " << numrel <<"\n";
}

// generic hexagon relation: indices i,j,k,l,m,n and u,x1,y1,x2,y2 such that
// gamma = M_i'*T^x1*M_k'*T^y1(alpha[m]) = T^u*M_j'*T^x2*M_l'*T^y2(alpha[n]).

// The hexagon has vertices {beta_1, alpha_i, oo, u+alpha[j], beta_2, gamma}
// where beta1 = M_i'(x1+alpha[k]), beta2 = T^u*M_j'(x2+alpha[l]),

// Edges (3 in each orientation)

// +{gamma, beta1} = (M_i'*T^x1*M_k'*T^y1)_m
// +{beta1, alpha_i} = (M_i'*T^x1)_k
// +{alpha_i, oo} = (I)_i
// {oo, u+alpha_j} = - (T^u)_j
// {u+alpha_j, beta2} = - (T^u*M_j'*T^x2)_l
// {beta2, gamma} = - (T^u*M_j'*T^x2*M_l'*T^y2)_n

void face_relations::general_hexagon_relation(const POLYGON& H, int check)
{
  const vector<int>& hex = H.first;
  const vector<Quad>& ux1y1x2y2 = H.second;
  int i=hex[0], j=hex[1], k=hex[2], l=hex[3], m=hex[4], n=hex[5];
  Quad u = ux1y1x2y2[0], x1 = ux1y1x2y2[1], y1 = ux1y1x2y2[2], x2 = ux1y1x2y2[3], y2 = ux1y1x2y2[4];
  int symmetry = 0; // work this out later
  if(verbose)
    {
      cout << "Applying ";
      //      if (symmetry) cout<<"symmetric ";
      cout<<"hexagon relation "<<hex<<" (u,x1,y1,x2,y2) = "<<ux1y1x2y2<<"\n";
    }

  vector<int> types = {m,k,i,j,l,n}, signs = {1,1,1,-1,-1,-1};
  mat22 M2 = M_alphas[alpha_inv[i]] * mat22::Tmat(x1);
  mat22 M1 = M2 * M_alphas[alpha_inv[k]] * mat22::Tmat(y1);
  mat22 M3 = mat22::Tmat(u);
  mat22 M4 = M3 * M_alphas[alpha_inv[j]] * mat22::Tmat(x2);
  mat22 M5 = M4 * M_alphas[alpha_inv[l]] * mat22::Tmat(y2);
  vector<action> Mops = {action(P1, M1),
                         action(P1, M2),
                         action(P1, mat22::identity),
                         action(P1, M3),
                         action(P1, M4),
                         action(P1, M5)};

  general_relation(Mops, types, signs, symmetry, check);

  if(verbose)
    cout << "After hexagon relation "<< hex <<", number of relations = " << numrel <<"\n";
}

void face_relations::solve_relations()
{
  vec npivs; // pivs is a class attribute
  if(verbose>1)
    {
      mat M = relmat.as_mat().slice(numrel,ngens);
      cout<<"relmat = "<<M<<endl;
      if (characteristic)
        {
          vec pcols, npcols;
          long rk_modp, ny_modp;
          echmodp_uptri(M, pcols, npcols, rk_modp, ny_modp, characteristic);
          cout<<"rank_mod_p(relmat) = "<<rk_modp;
        }
      else
        {
          cout<<"rank(relmat) = "<<relmat.rank();
        }
      cout<<", ngens = "<<ngens<<endl;
    }
#ifdef USE_SMATS
#ifdef TIME_RELATION_SOLVING
  if (verbose)
    {
      cout<<"\nStarting to solve relation matrix of size "<<numrel<<" x "<<ngens<<", density "<<density(relmat)<<"\n"<<flush;
      cout<<"# nonzero entries in row\t\t#rows"<<endl;
      for (int i=1; i<=maxrowsize; i++)
        cout<<i<<"\t\t\t\t\t"<<relmat_rowdata[i]<<endl;
    }
   timer t;
   t.start("relation solver");
#endif
   scalar modulus = (characteristic==0? DEFAULT_MODULUS: characteristic);
   smat_elim sme(relmat, modulus);
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
   int ok = 1;
   if (characteristic==0)
     {
       ok = liftmat(ker,MODULUS,sp,d1);
     }
   else
     {
       d1 = 1;
       sp = ker;
     }
   if (!ok)
     {
       if(verbose)
         cout << "failed to lift modular kernel using modulus "
              << MODULUS << endl;
#ifdef USE_CRT
       int mod2 = 1073741783; // 2^30-41
       INT mmod(MODULUS); mmod*=mod2;
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
       hmod = characteristic;
       denom = (characteristic==0? d1: 1);
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
   for(long i=1; i<=ngens; i++)
     coord.setrow(i,sp.row(i).as_vec());
#ifdef USE_CRT
   long maxcoord =0;
   for(long i=1; i<=ngens; i++)
     for(long j=1; j<=rk; j++)
       {
         long cij = coord(i,j);
         if (abs(cij)>maxcoord) maxcoord=cij;
       }
   if (verbose)
     cout << "Max entry in coord is "<<maxcoord<<endl;
#endif
   // if hmod>0, coord is only defined modulo hmod
   sp=smat(0,0); // clear space
#else // not using smats; not yet implemented for characteristic>0
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

