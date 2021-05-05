// FILE HOMSPACE.CC: Implemention of class homspace

//#define USE_CRT // if using smats  mod MODULUS, try CRT-ing with another prime
                // NB this is experimental only

#include <eclib/msubspace.h>
#include <eclib/xmod.h>
#include "homspace.h"
#include <assert.h>
const string W_opname("W");
const string T_opname("T");

int liftmats_chinese(const smat& m1, scalar pr1, const smat& m2, scalar pr2,
                     smat& m, scalar& dd);

int sign(int a)
{
  if (a) return (a>0? 1: -1);
  return 0;
}

// In add_rel(rel), rel is a list of (positive) (c:d)-symbol numbers i
// such that the corresponding symbols add to 0 in homology.  We use
// the map i -> j=coordindex[i] to convert this to a vector of
// dimension ngens, and append that as a new row to the relation
// matrix relmat.

void homspace::add_rel(const vector<int>& rel, const vector<int>& types)
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

homspace::homspace(const Quad& n, int hp, int cuspid, int verb) :symbdata(n)
{
  verbose=verb;
  cuspidal=cuspid;
  hmod = 0;
  if (verbose) symbdata::display();
  plusflag=hp;                  // Sets static level::plusflag = hp
  ngens=0;
  coordindex.resize(nsymbx);
  gens.resize(nsymbx+1);
  //NB start of gens array is at 1 not 0

  edge_relations(); // sets coordindex

  if (verbose)
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
            cout << "Type " << j+1 << ", alpha = "<<alphalist[j]<<":\n";
          for (i=0; i<nsymb; i++)
            cout << i<<":\t"<<symbol(i)<<"\t"<<coordindex[i+offset] << "\n";
        }
      cout << endl;
    }

  face_relations(); // fills relmat with the relations

  if(verbose)
    {
      cout << "Finished face relations: ";
      cout << "number of relations = " << numrel;
      cout << " (bound was "<<maxnumrel<<")"<<endl;
    }

  solve_relations();

  if (verbose)
    {
      cout << "dimension (relative to cusps) = " << rk << endl;
    }

  kernel_delta();

  if(verbose)
    {
      cout << "number of cusps = " << ncusps << endl;
      if (cuspidal)
        cout << "dimension = " << dimension << endl;
    }

  make_freemods();

  if (verbose) cout << "Finished constructing homspace.\n";
}

// 2-term (edge) relations

void homspace::edge_relations()    // computes coordindex
{
  Quad unit = fundunit;
  long lenrel = Quad::nunits;
  if(!plusflag) {unit=fundunit*fundunit; lenrel/=2;}
  symbop eps(this,unit,0,0,1);  assert (eps.det()==unit);
  symbop sof(this,0,-1,1,0);  assert (sof.det()==1);
  vector<int> a(lenrel), b(lenrel);
  vector<int> check(nsymb, 0);
  int j, k, triv;
  if(verbose) cout << "About to start on 2-term (edge) relations.\n";
  if(verbose && n_alphas>1)
    cout<<"Generic edge relations for type 0 symbols\n";
  for (j=nsymb-1; j>=0; j--)
    {
      if (check[j]==0)
        { 
	  if(verbose>1) cout << "j = " << j << ":\t";
          a[0]=j; b[0]=sof(j); triv=(j==b[0]);
          for(k=1; k<lenrel; k++)
            {
              a[k]= eps(a[k-1]);
              b[k]= eps(b[k-1]);
              triv= triv | (j==b[k]);
            }
          for (k=0; k<lenrel; k++) check[a[k]]=check[b[k]]=1;
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
  int field = Quad::d;
  Quad w = Quad::w;
  int j, k, l, u=(field-3)/8; // u=2, 5, 8, 20 for 19,43,67,163

  // relevant alphas are  {1:w/2, 2:(w-1)/2}

  // (1) (g)_w/2 = -(gK)_(w-1)/2 with K = [w-1,u;2,-w], det=1,  order 3
  // (2) (g)_w/2 = (gL)_w/2      with L = [-1,w;0,1],  det=-1, order 2, if plus
  //             = -(gLK)_(w-1/2)
  symbop K(this, w-1,u,2,-w); assert (K.det()==1);
  symbop L(this, -1,w,0,1);   assert (L.det()==-1);

  vector<int> check(nsymb, 0);
  int offset1 = nsymb, offset2 = 2*nsymb;
  for (j=0; j<nsymb; j++)
    {
      if (check[j]==0)
        {
          k = K(j);
          check[j] = 1;
          gens[++ngens] = nsymb+j;
          coordindex[offset1 + j] = ngens;
          coordindex[offset2 + k] = -ngens;
          if (plusflag)
            {
              l = L(j);
              // if (verbose)
              //   cout<<symbol(j)<<" -L-> "<<symbol(l)<<endl;
              check[l] = 1;
              coordindex[offset1 + l] = ngens;
              coordindex[offset2 + K(l)] = -ngens;
            }
        }
    }
}

void homspace::edge_relations_3()    // extra edge relations for alphas with denominator 3
{
  int field = Quad::d;
  Quad w = Quad::w;
  int j, jj, k, l, u=-(field+5)/12; // u=-4, -6, -24 for 43,67,163

  // relevant alphas are  {3:w/3, 4:-w/3, 5:(1-w)/3, 6:(w-1)/3, 7:(1+w)/3, 8:(-1-w)/3}

  symbop J(this, -1,0,0,1);            assert (J.det()==-1);
  symbop L1(this, w,-u,3,1-w);         assert (L1.det()==-1);
  symbop K1(this, w,u,3,w-1);             assert (K1.det()==1);
  symbop K3(this, 1+w,-(w+u+1),3,-(1+w)); assert (K3.det()==1);
  symbop K4(this, -(1+w),-(w+u+1),3,1+w); assert (K4.det()==1);

  // (1) type 7, alpha=(1+w)/3, antisymmetric via K3:

  int offset7 = 7*nsymb;
  vector<int> check(nsymb, 0);
  for (j=0; j<nsymb; j++)
    {
      if (check[j]==0)
        {
          check[j] = 1;
          k = K3(j);
          if (j==k) // symbol trivial
            {
              coordindex[offset7+j] = 0;
            }
          else
            {
              check[k] = 1;
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
          coordindex[offset8+J(j)] = coordindex[offset7+j];
        }
    }
  else
    {
      std::fill(check.begin(), check.end(), 0);
      for (j=0; j<nsymb; j++)
        {
          if (check[j]==0)
            {
              check[j] = 1;
              k = K4(j);
              if (j==k) // symbol trivial
                {
                  coordindex[offset8+j] = 0;
                }
              else
                {
                  check[k] = 1;
                  gens[++ngens] = offset8+j;
                  coordindex[offset8+j] = ngens;
                  coordindex[offset8+k] = -ngens;
                }
            }
        }
    }

  // (3) types 3,4,5,6: identify in 4-tuples up to sign if (+), else in pairs if (0)

  int offset3 = 3*nsymb, offset4 = 4*nsymb, offset5 = 5*nsymb, offset6 = 6*nsymb;
  std::fill(check.begin(), check.end(), 0);
  for (j=0; j<nsymb; j++) // index of type 3 symbol
    {
      if (check[j]==0)
        {
          jj = J(j); // index of type 4 symbol
          k = K1(j); // index of type 5 symbol
          l = L1(j); // index of type 6 symbol
          check[j] = check[jj] = 1;
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
  vector<int> rel(3), types(3,0), check(nsymb, 0);
  long j, k;
  symbop TS(this,1,1,-1,0); assert (TS.det()==1);
  symbop R(this,0,1,1,0);   assert (R.det()==-1);
  for (k=0; k<nsymb; k++)
    if (check[k]==0)
      {
        for(j=0; j<3; j++)
          {
            rel[j] = (j? TS(rel[j-1]): k);
            check[rel[j]] = 1;
            if (plusflag)
              check[R(rel[j])] = 1;
          }
        add_rel(rel, types);
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
  vector<int> rel(3), types(3,0), check(nsymb, 0);
  long j, k;

  Quad w(0,1);
  long field = Quad::d;
  symbop X = (field==1? symbop(this,w,1,1,0): symbop(this,1,w,w-1,0));
  assert (X.det()==(field==1? -1: 1));

  for (k=0; k<nsymb; k++)
    if (check[k]==0)
      {
        for(j=0; j<3; j++)
          {
            rel[j] = (j? X(rel[j-1]): k);
            check[rel[j]] = 1;
          }
        add_rel(rel, types);
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
  vector<int> rel(4), types(4,0), check(nsymb, 0);
  long j, k;

  Quad w(0,1);
  symbop U(this,w,1,1,0);  assert (U.det()==-1);
  symbop S(this,0,-1,1,0); assert (S.det()==1);
  symbop J(this,fundunit,0,0,1);  assert (J.det()==fundunit);

  for (k=0; k<nsymb; k++)
    if (check[k]==0)
      {
        for(j=0; j<4; j++)
          {
            rel[j] = (j? U(rel[j-1]): k);
            if(plusflag) // NB det(U)=-1
              {
                check[rel[j]] = 1;
                check[S(rel[j])] = 1;
              }
            else
              {
                if(j%2==0)
                  check[rel[j]] = 1;
              }
          }
        if (!plusflag) // since det(U)=-1
          {
            rel[1] = J(rel[1]);
            rel[3] = J(rel[3]);
          }
        add_rel(rel, types);
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
  vector<int> rel(4), types(4,0), check(nsymb, 0);
  long j, k;
  Quad w(0,1);

  symbop Y(this,1,-w,1-w,-1);  assert (Y.det()==1);
  symbop USof(this,w,-1,1,0);  assert (USof.det()==1);
  symbop R(this,0,1,1,0);      assert (R.det()==-1);

  for (k=0; k<nsymb; k++)
    if (check[k]==0)
      {
        for(j=0; j<4; j+=2) // j=0,2; j+1=1,3
          {
            rel[j] = (j? Y(rel[j-2]): k);
            check[rel[j]] = 1;
            rel[j+1] = USof(rel[j]);
            if (plusflag)
              check[R(rel[j+1])]=1;
          }
        add_rel(rel, types);
      }
}

// extra hexagon relation for field 11
void homspace::hexagon_relation_11()
{
  if(verbose)
    {
      cout << "Face relation type 2 (hexagons):\n";
    }
  vector<int> rel(6), types(6,0), check(nsymb, 0);
  long j, k;
  Quad w(0,1);

  symbop X(this,1,-w,1-w,-2);
  assert (X.det()==1);
  symbop USof(this,w,-1,1,0);
  assert (USof.det()==1);
  symbop R(this,0,1,1,0);
  assert (R.det()==-1);

  for (k=0; k<nsymb; k++)
    if (check[k]==0)
      {
        for(j=0; j<6; j+=2) // j=0,2,4; j+1=1,3,5
          {
            rel[j] = (j? X(rel[j-2]): k);
            check[rel[j]] = 1;
            rel[j+1] = USof(rel[j]);
            if(plusflag)
              check[R(rel[j+1])]=1;
          }
        add_rel(rel, types);
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

  symbop K(this, w-1,u,2,-w);   assert (K.det()==1); // oo --> (w-1)/2 --> w/2 --> oo
  symbop N(this, 1+w,u-w,2,-w); assert (N.det()==1); // oo --> (w+1)/2 --> w/2 --> oo

  // N is the conjugate of K by [-1,w;0,1] which maps the first
  // triangle to the second with determinant -1.  Both have order 3 so
  // cycle the edges of each triangle around.

  vector<int> rel(3), types(3, 2), check(nsymb, 0);
  for (k=0; k<nsymb; k++)
    if (check[k]==0)
      {
        for(j=0; j<3; j++)
          {
            rel[j] = (j? K(rel[j-1]): k);
            check[rel[j]] = 1;
          }
        add_rel(rel, types);
      }
  if (!plusflag) // there's a second triangle (image of previous under L)
    {
      std::fill(check.begin(), check.end(), 0);
      for (k=0; k<nsymb; k++)
        if (check[k]==0)
          {
            for(j=0; j<3; j++)
              {
                rel[j] = (j? N(rel[j-1]): k);
                check[rel[j]] = 1;
              }
            add_rel(rel, types);
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
  int field = Quad::d;
  Quad w = Quad::w;
  int j, k, u=(field-3)/8; // u=2, 5, 8, 20 for 19,43,67,163

  symbop K(this, w-1,u,2,-w);   assert (K.det()==1);
  symbop S(this, 0,-1,1,0);     assert (S.det()==1);

  vector<int> rel(4), types(4), check(nsymb, 0);
  types[0] = types[1] = 0;
  types[2] = types[3] = 1;
  for (j=0; j<nsymb; j++)
    if (check[j]==0)
      {
        rel[0] = j;
        rel[2] = S(j);
        rel[3] = k = K(j);
        rel[1] = S(k); // = KS(j)
        check[j] = check[k] = 1;
        add_rel(rel, types);
      }
  if(verbose)
    {
      cout << "After type 0022 square relations, number of relations = " << numrel <<"\n";
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
  symbop M2(this, w+1,u,3,2-w); // maps {(w-1)/3,oo} to {w/3, (w+1)/3}
  types = {5,6,7};  // {(1-w)/3, oo}, {(w-1)/3, oo}, {(1+w)/3, oo}

  for (k=0; k<nsymb; k++)
    {
      rel[0] = M1(k);
      rel[1] = M2(k);
      rel[2] = k;
      add_rel(rel, types);
    }

  if (!plusflag) // there's a second triangle (image of previous under J)
    {
      // Triangle {oo, -w/3, -(w+1)/3}

      symbop M3(this, w,u,-3,w-1);    // = J*M1*J, maps {(w-1)/3,oo} to {oo, -w/3}
      symbop M4(this, w+1,-u,-3,2-w); // = J*M2*J, maps {(1-w)/3,oo} to {-w/3, -(w+1)/3}
      types = {6,5,8};  // {(w-1)/3, oo}, {(1-w)/3, oo}, {-(1+w)/3, oo}

      for (k=0; k<nsymb; k++)
        {
          rel[0] = M3(k);
          rel[1] = M4(k);
          rel[2] = k;
          add_rel(rel, types);
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
  vector<int> rel(4), types(4);

  // (1)
  // vertices {0, oo, w/3, 3w/11}
  // edge types 0 {0,oo}, 5 {(1-w)/3,oo}, 1 {w/2,oo}, 5 under I, V, N4, N4*V

  symbop N4(this, -3,w,w-1,4);   assert (N4.det()==-1);
  symbop V(this, w,-4,3,w-1);     assert (V.det()==1);
  types = {0,5,1,5};

  types[0] = 0; types[2] = 1;
  types[1] = types[3] = 5;
  for (j=0; j<nsymb; j++)
    {
      rel[0] = j;
      rel[1] = V(j);
      rel[2] = k = N4(j);
      rel[3] = V(k);
      add_rel(rel, types);
    }

  // (2)
  // vertices {w/2,oo,(w+1)/3,(5*w+3)/13}
  // edge types 1 {w/2,oo}, 7 {(w+1)/3,oo}, 1, 7 under I, U, N5, N5*U

  symbop N5(this, w-4,w+5,w+1,4-w);   assert (N5.det()==1);
  symbop U(this, w+1,3-w,3,-1-w);     assert (U.det()==1);
  types = {1,7,1,7};

  for (j=0; j<nsymb; j++)
    {
      rel[0] = j;
      rel[1] = U(j);
      rel[2] = k = N5(j);
      rel[3] = U(k);
      add_rel(rel, types);
    }

  // (3) image of (2) under L=[-1,0;0,1] with det(L)=-1, only if !posflag

  // vertices {w/2, oo, (2*w-1)/3, (8*w-3)/13}
  // edge types 1 {w/2,oo}, 8 {-(w+1)/3,oo}, 1, 7 under I, W1, W2, W2*W1

  symbop W1(this, 2*w-1,w-8,3,w+1);   assert (W1.det()==1);
  symbop W2(this, w-7,4*w+5,w+1,7-w); assert (W2.det()==1);
  types = {1,8,1,8};

  for (j=0; j<nsymb; j++)
    {
      rel[0] = j;
      rel[1] = W1(j);
      rel[2] = k = W2(j);
      rel[3] = W1(k);
      add_rel(rel, types);
    }

  if(verbose)
    {
      cout << "After type extra square relations, number of relations = " << numrel <<"\n";
    }
}

void homspace::solve_relations()
{
   vec pivs, npivs;
   int i;
#if(USE_SMATS)
   smat_elim sme(relmat);
   int d1;
   smat ker = sme.kernel(npivs,pivs), sp;
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
       denom1 = d1;
     }
   else
     {
       hmod = MODULUS;
       denom1 = 1;
     }
   relmat=smat(0,0); // clear space
   if(verbose>1)
     {
       cout << "kernel of relmat = " << sp.as_mat() << endl;
       cout << "pivots = "<<pivs << endl;
       cout << "denom = "<<d1 << endl;
     }
   rk = sp.ncols();
   coord.init(ngens+1,rk); // 0'th is unused
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
  denom1 = I2int(denom(sp));
  relmat.init(); sp.clear(); 
#endif
  if (verbose)
    {
      cout << "rk = " << rk << endl;
      cout << "coord:" << coord;
      if (hmod)
	cout << "failed to lift, coord is only defined modulo "<<hmod<<endl;
      else
        cout << "lifted ok, denominator = " << denom1 << endl;
      cout << "pivots = " << pivs <<endl;
    }
  if (rk>0)
    {
      freegens.resize(rk);
      for (i=0; i<rk; i++) freegens[i] = gens[pivs[i+1]];
      if (verbose)
        {
          cout << "freegens: ";
          for (i=0; i<rk; i++) cout << freegens[i] << " ";
          cout << endl;
        }
    }
}

void homspace::kernel_delta()
{
  if (verbose)
    cout<<"Computing boundary map"<<endl;
  cusplist cusps(2*rk,this);
  mat deltamat(2*rk,rk);
  int i, j, s, t=0;
  modsym m;
  for (i=0; i<rk; i++)
    {
      s = j = freegens[i];
      if (n_alphas>1)
        {
          //      cout<<"j = "<<j<<": ";
          std::div_t st = div(j, nsymb);
          s = st.rem;  // remainder gives (c:d) symbol number
          t = st.quot; // quotient gives symbol type
          //      cout<<"(s,t) = ("<<s<<","<<t<<")\n";
        }
      m = modsym(symbol(s), t);
      deltamat(cusps.index(m.beta())+1, i+1) += 1;  // N.B. offset of 1
      deltamat(cusps.index(m.alpha())+1, i+1) -= 1;
    }
  ncusps=cusps.count();

  kern = kernel(smat(deltamat));
  vec pivs, npivs;
  int d2;
  smat sk;
  int ok = liftmat(smat_elim(deltamat).kernel(npivs,pivs),MODULUS,sk,d2);
  if (!ok)
    cout << "**!!!** failed to lift modular kernel\n" << endl;

  tkernbas = transpose(kern.bas());         // dim(kern) x rank
  if(verbose>1)
    cout<<"tkernbas = "<<tkernbas.as_mat()<<endl;

  dimension = (cuspidal? dim(kern): rk);
  denom2 = d2;
  denom3 = denom1 * denom2;
}

void homspace::make_freemods()
{
  const smat& basiskern = basis(kern);
  if (verbose)  cout << "Freemods:\n";
  int i,j,s,t=0;
  modsym m;
  for (i=0; i<rk; i++)
    {
      s = j = freegens[i];
      if (n_alphas>1)
        {
          //      cout<<"j = "<<j<<": ";
          std::div_t st = div(j, nsymb);
          s = st.rem;  // remainder gives (c:d) symbol number
          t = st.quot; // quotient gives symbol type
          //      cout<<"(s,t) = ("<<s<<","<<t<<")\n";
        }
      m = modsym(symbol(s), t);
      freemods.push_back(m);
      int n = (cuspidal? ! trivial(basiskern.row(i+1).as_vec()) : 1);
      needed.push_back(n);
      if (verbose)
        {
          cout << i << ": " << m;
          if (!n) cout << " (not needed)";
          cout << endl;
        }
    }
  if (verbose)
    {
      cout << "Basis of ker(delta):\n";
      cout << basiskern;
      cout << "pivots: " << pivots(kern) << endl;
    }
}


vec homspace::chain(const symb& s) const  //=old getcoord
{
 long i= coordindex[index(s)];
 vec ans(coord.ncols());
 if (i) ans = sign(i)*(coord.row(abs(i)));
 return ans;
}

vec homspace::chaincd(const Quad& c, const Quad& d) const //=old getcoord2
{
 long i= coordindex[index2(c,d)];
 vec ans(coord.ncols());
 if (i) ans = sign(i)*(coord.row(abs(i)));
 return ans;
}

vec homspace::projchaincd(const Quad& c, const Quad& d) const 
{
 long i= coordindex[index2(c,d)];
 vec ans(projcoord.ncols());
 if (i) ans = sign(i)*(projcoord.row(abs(i)));
 return ans;
}

vec homspace::chain(const Quad& nn, const Quad& dd) const //=old qtovec2
{
   vec ans = chaincd(0,1), part;
   Quad c=0, d=1, e, a=nn, b=dd, q, f;
   while (b!=0)
   { q=a/b;
     f=a; a=-b; b= f-q*b;
     e=d; d= c; c= e+q*c; c=-c;
     part = chaincd(c,d);
     if(hmod)
       ans.addmodp(part,hmod);
     else
       ans += part;
   }
   if(hmod) ans=reduce_modp(ans,hmod);
   return ans;
}

vec homspace::projcycle(const Quad& nn, const Quad& dd) const  //
{
  vec ans = projchaincd(0,1), part;
   Quad c=0, d=1, e, a=nn, b=dd, q, f;
   while (b!=0)
   { q=a/b;
     f=a; a=-b; b= f-q*b;
     e=d; d= c; c= e+q*c; c=-c;
     part = projchaincd(c,d);
     if(hmod)
       ans.addmodp(part,hmod);
     else
       ans += part;
   }
   if(hmod) ans=reduce_modp(ans,hmod);
   return ans;
}

vec homspace::applyop(const matop& mlist, const RatQuad& q) const
{ vec ans(rk), part;
  long i=mlist.length();
  while (i--)
    {
      part = chain(mlist[i](q));
      if(hmod)
        ans.addmodp(part,hmod);
      else
        ans += part;
    }
  if(hmod) ans=reduce_modp(ans,hmod);
  return ans;
}

mat homspace::calcop(const string opname, const Quad& p, const matop& mlist, int dual, int display) const
{
  mat m(rk,rk);
  for (long j=0; j<rk; j++) if (needed[j])
     { vec colj = applyop(mlist,freemods[j]);
       if(hmod) colj=reduce_modp(colj,hmod);
       m.setcol(j+1,colj);
     }
  if(cuspidal) m = restrict_mat(smat(m),kern).as_mat();
  if(dual) m = transpose(m);
  if (display) cout << "Matrix of " << opname << "(" << p << ") = " << m;
  if (display && (dimension>1)) cout << endl;
  return m;
}

vec homspace::calcop_col(const string opname, const Quad& p, const matop& mlist, int j, int display) const
{
  vec colj = applyop(mlist,freemods[j-1]);
  if(hmod) colj=reduce_modp(colj,hmod);
  return colj;
}

mat homspace::calcop_cols(const string opname, const Quad& p, const matop& mlist, const vec& jlist, int display) const
{
  int i, j, d = dim(jlist);
  mat m(d,rk);
  for (i=1; i<=d; i++)
    {
      j = jlist[i];
      vec colj = applyop(mlist,freemods[j-1]);
      if(hmod) colj=reduce_modp(colj,hmod);
      m.setcol(i,colj);
     }
  return m;
}
 
smat homspace::s_calcop_cols(const string opname, const Quad& p, const matop& mlist, const vec& jlist, int display) const
{
  int i, j, d = dim(jlist);
  smat m(d,rk);
  for (i=1; i<=d; i++)
    {
      j = jlist[i];
      svec colj = applyop(mlist,freemods[j-1]);
      if(hmod) colj.reduce_mod_p(hmod);
      m.setrow(i,colj);
     }
  return m;
}
 
smat homspace::s_calcop(const string  opname, const Quad& p, const matop& mlist, 
			int dual, int display) const
{
  smat m(rk,rk);
  for (long j=0; j<rk; j++) if (needed[j])
     { svec colj = applyop(mlist,freemods[j]);
       if(hmod) colj.reduce_mod_p(hmod);
       m.setrow(j+1,colj);
     }
  if(cuspidal)
    {
      m = restrict_mat(transpose(m),kern);
      if(dual) m = transpose(m);
    }
  else if(!dual) {m=transpose(m);}
  if (display) 
    {
      cout << "Matrix of " << opname << "(" << p << ") = ";
      if (dimension>1) cout << "\n";
      cout<<m.as_mat();
    }
  return m;
}

mat homspace::calcop_restricted(const string opname, const Quad& p, const matop& mlist, const subspace& s, int dual, int display) const
{
  long d=dim(s);
  mat m(d,rk);
  for (long j=0; j<d; j++)
     {
       long jj = pivots(s)[j+1]-1;
       vec colj = applyop(mlist,freemods[jj]);
       if(hmod) colj=reduce_modp(colj,hmod);
       m.setrow(j+1,colj);
     }
  if(hmod)
    m = matmulmodp(m,basis(s),hmod);
  else
    m = m*basis(s);
  if(!dual) m=transpose(m); // dual is default for restricted ops
  if (display) cout << "Matrix of " << opname << "(" << p << ") = " << m;
  if (display && (dimension>1)) cout << endl;
  return m;
}
 
smat homspace::s_calcop_restricted(const string opname, const Quad& p, const matop& mlist, const ssubspace& s, int dual, int display) const
{
  long d=dim(s);
  smat m(d,rk);
  for (long j=1; j<=d; j++)
     { 
       long jj = pivots(s)[j];
       svec colj = applyop(mlist,freemods[jj-1]);
       if(hmod) colj.reduce_mod_p(hmod);
       m.setrow(j,colj);
     }
  m = mult_mod_p(m,basis(s),MODULUS);
  m.reduce_mod_p();
  if(!dual) m=transpose(m); // dual is default for restricted ops
  if (display)
    {
      cout << "Matrix of " << opname << "(" << p << ") = " << m.as_mat();
      if (dimension>1) cout << endl;
    }
  return m;
}
 
mat homspace::heckeop(const Quad& p, int dual, int display) const
{
 matop matlist(p,modulus);
 string name = (div(p,modulus)) ? W_opname : T_opname;
 return calcop(name,p,matlist,dual,display);
}
 
vec homspace::heckeop_col(const Quad& p, int j, int display) const
{
 matop matlist(p,modulus);
 string name = (div(p,modulus)) ? W_opname : T_opname;
 return calcop_col(name,p,matlist,j,display);
}
 
mat homspace::heckeop_cols(const Quad& p, const vec& jlist, int display) const
{
 matop matlist(p,modulus);
 string name = (div(p,modulus)) ? W_opname : T_opname;
 return calcop_cols(name,p,matlist,jlist,display);
}
 
smat homspace::s_heckeop_cols(const Quad& p, const vec& jlist, int display) const
{
 matop matlist(p,modulus);
 string name = (div(p,modulus)) ? W_opname : T_opname;
 return s_calcop_cols(name,p,matlist,jlist,display);
}
 
smat homspace::s_heckeop(const Quad& p, int dual, int display) const
{
 matop matlist(p,modulus);
 string name = (div(p,modulus)) ? W_opname : T_opname;
 return s_calcop(name,p,matlist,dual,display);
}
 
mat homspace::heckeop_restricted(const Quad& p, const subspace& s, int dual, int display) const
{
 matop matlist(p,modulus);
 string name = (div(p,modulus)) ? W_opname : T_opname;
 return calcop_restricted(name,p,matlist,s,dual,display);
}
 
smat homspace::s_heckeop_restricted(const Quad& p, const ssubspace& s, int dual, int display) const
{
 matop matlist(p,modulus);
 string name = (div(p,modulus)) ? W_opname : T_opname;
 return s_calcop_restricted(name,p,matlist,s,dual,display);
}
 
mat homspace::wop(const Quad& q, int dual, int display) const
{
 matop matlist(q,modulus);
 return calcop(W_opname,q,matlist,dual,display);
}
 
mat homspace::fricke(int dual, int display) const
{
 matop frickelist(modulus,modulus);
 return calcop(W_opname,modulus,frickelist,dual,display);
}

mat homspace::opmat(int i, int dual, int verb)
{
  if((i<0)||(i>=nap)) return mat(dimension);  // shouldn't happen
  Quad p = primelist[i];
  if(verbose) 
      cout<<"Computing " << ((div(p,modulus)) ? W_opname : T_opname) <<"("<<p<<")...";
  return heckeop(p,dual,verb); // Automatically chooses W or T
}

vec homspace::opmat_col(int i, int j, int verb)
{
  if((i<0)||(i>=nap)) return vec(dimension);  // shouldn't happen
  Quad p = primelist[i];
  if(verbose) 
      cout<<"Computing " << ((div(p,modulus)) ? W_opname : T_opname) <<"("<<p<<")...";
  return heckeop_col(p,j,verb); // Automatically chooses W or T
}

mat homspace::opmat_cols(int i, const vec& jlist, int verb)
{
  if((i<0)||(i>=nap)) return mat(dimension);  // shouldn't happen
  Quad p = primelist[i];
  if(verbose) 
      cout<<"Computing " << ((div(p,modulus)) ? W_opname : T_opname) <<"("<<p<<")...";
  return heckeop_cols(p,jlist,verb); // Automatically chooses W or T
}

smat homspace::s_opmat_cols(int i, const vec& jlist, int verb)
{
  if((i<0)||(i>=nap)) return smat(dimension);  // shouldn't happen
  Quad p = primelist[i];
  if(verbose)
    cout<<"Computing " << ((div(p,modulus)) ? W_opname : T_opname) <<"("<<p<<")..."<<flush;
  return s_heckeop_cols(p,jlist,verb); // Automatically chooses W or T
}

mat homspace::opmat_restricted(int i, const subspace& s, int dual, int verb)
{
  if((i<0)||(i>=nap)) return mat(dim(s));  // shouldn't happen
  Quad p = primelist[i];
  if(verbose) 
      cout<<"Computing " << ((div(p,modulus)) ? W_opname : T_opname) <<"("<<p
	  <<") restricted to subspace of dimension "<<dim(s)<<" ..."<<flush;
  return heckeop_restricted(p,s,dual,verb); // Automatically chooses W or T
}

smat homspace::s_opmat(int i, int dual, int v)
{
  //  if(i==-1) return s_conj(dual,v);
  if((i<0)||(i>=nap)) 
    {
      return smat(dimension);  // shouldn't happen
    }
  Quad p = primelist[i];
  if(v) 
    {
      cout<<"Computing " << ((div(p,modulus)) ? W_opname : T_opname) <<"("<<p<<")...";
      smat ans = s_heckeop(p,dual,0); // Automatically chooses W or T
      cout<<"done."<<endl;
      return ans;
    }
  else return s_heckeop(p,dual,0); // Automatically chooses W or T
}

smat homspace::s_opmat_restricted(int i, const ssubspace& s, int dual, int v)
{
  if((i<0)||(i>=nap)) 
    {
      return smat(dim(s));  // shouldn't happen
    }
  Quad p = primelist[i];
  if(v) 
    {
      cout<<"Computing " << ((div(p,modulus)) ? W_opname : T_opname) <<"("<<p
	  <<") restricted to subspace of dimension "<<dim(s)<<" ..."<<flush;
      smat ans = s_heckeop_restricted(p,s,dual,0); // Automatically chooses W or T
      cout<<"done."<<endl;
      return ans;
    }
  else return s_heckeop_restricted(p,s,dual,0); // Automatically chooses W or T
}

vector<long> homspace::eigrange(long i)  // implementing virtal function in matmaker
{
  vector<long> ans;
  if((i<0)||(i>=nap)) return ans;  // shouldn't happen
  Quad p = primelist[i];
  long normp = quadnorm(p);
  if (verbose) 
    cout << "eigrange for p = " << p << ":\t";
  if(div(p,modulus))
    {
      vector<long> ans(2);
      ans[0]=-1;
      ans[1]=1;
      if (verbose) 
	cout << ans << endl;
      return ans;
    }
  else
    {
      long aplim=2;
      while (aplim*aplim<=4*normp) aplim++;
      aplim--;
      if(verbose)
	cout << "|ap| up to "<<aplim<<":\t";
      long ap, l = 2*aplim+1;
      vector<long> ans(l);
      ans[0]=0;
      for(ap=-aplim; ap<=aplim; ap++)
	ans[ap+aplim] = ap;
      if (verbose) 
	cout << ans << endl;
      return ans;
    }
}

vec homspace::maninvector(const Quad& p) const
{
  vector<Quad> resmodp=residues(p);
  vec ans = chain(0,p), part;             // =0, but sets the right length.
  vector<Quad>::const_iterator res=resmodp.begin();
  while(res!=resmodp.end())
    {
      part = chain(*res++,p);
      if(hmod)
        ans.addmodp(part,hmod);
      else
        ans += part;
    }
  if(hmod) ans=reduce_modp(ans,hmod);
  return ans;
}

vec homspace::manintwist(const Quad& lambda, const vector<Quad>& res, vector<int> chitable) const
{
  vec ans = chain(0,lambda), part;          // =0, but sets the right length.
  vector<int>::const_iterator chi=chitable.begin();
  vector<Quad>::const_iterator r=res.begin();
  while(r!=res.end())
   {
     part = (*chi++)*chain(*r++,lambda);
      if(hmod)
        ans.addmodp(part,hmod);
      else
        ans += part;
   }
  if(hmod) ans=reduce_modp(ans,hmod);
 return ans;
}

vec homspace::projmaninvector(const Quad& p) const    // Will only work after "proj"
{
  vector<Quad> resmodp=residues(p);
  vec ans = projcycle(0,p), part;         // =0, but sets the right length.
  vector<Quad>::const_iterator res=resmodp.begin();
  while(res!=resmodp.end())
    {
      part = projcycle(*res++,p);
      if(hmod)
        ans.addmodp(part,hmod);
      else
        ans += part;
    }
  if(hmod) ans=reduce_modp(ans,hmod);
  return ans;
}

vec homspace::newhecke(const Quad& p, const Quad& n, const Quad& d) const
                                     // Will only work after "proj"
{ 
  vec ans = projcycle(p*n,d), part;
  vector<Quad> resmodp=residues(p);  Quad dp = d*p;
  vector<Quad>::const_iterator res=resmodp.begin();
  while(res!=resmodp.end())
    {
      part = projcycle(n+d*(*res++),dp);
      if(hmod)
        ans.addmodp(part,hmod);
      else
        ans += part;
    }
  if(hmod) ans=reduce_modp(ans,hmod);
  return ans;
}

#if(USE_SMATS)
void mergeposval(long* pos, long* val, long& npos, long f)
{
  if(f==0) return;
  long af=abs(f), sf=sign(f), i=0,j;
  while((pos[i]<af)&&(i<npos)) i++;
  if(i==npos) //run off end
    {
      pos[npos]=af; val[npos]=sf; 
      npos++;
      return;
    }
  if(pos[i]==af) // change existing entry
    {
      val[i]+=sf;
      if(val[i]!=0) return;  // else close up
      for(j=i; j<npos; j++) {pos[j]=pos[j+1]; val[j]=val[j+1];}
      npos--;
      return;
    }
  // else pos[i-1] < af < pos[i], so insert
  for(j=npos; j>i; j--) {pos[j]=pos[j-1]; val[j]=val[j-1];}
  npos++;
  pos[i]=af; val[i]=sf;
  return;
}
#endif


#define DEBUG_CHINESE

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

vec reduce_modp(const vec& v, const scalar& p)
{
  int i, d=dim(v);
  scalar ai, p2 = p>>1;
  vec ans(d);
  for(i=1; i<=dim(v); i++)
    {
      ai = v[i]%p;
      while( ai>p2) ai-=p;
      while(-ai>p2) ai+=p;
      ans[i] = ai;
    }
  return ans;
}

mat reduce_modp(const mat& m, const scalar& p)
{
  int i, j, nr=m.nrows(), nc=m.ncols();
  scalar aij, p2 = p>>1;
  mat ans(nr,nc);
  for(i=1; i<=nr; i++)
    for(j=1; j<=nc; j++)
      {
        aij = m(i,j)%p;
        while( aij>p2) aij-=p;
        while(-aij>p2) aij+=p;
        ans(i,j) = aij;
      }
  return ans;
}
