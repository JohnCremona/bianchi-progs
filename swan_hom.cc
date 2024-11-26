// FILE SWAN_HOM.CC: implementation of functions for computing integral 1-homology

#include "swan_utils.h"
#include "swan_tess.h"
#include "swan_hom.h"
#include "pari_snf.h"
#include "flint_snf.h"

// Define this to use FLINT instead of PARI for HNF,SNF:
//#define INVARIANTS_VIA_FLINT

//#define CHECK_INVARIANTS

ostream& operator<<(ostream& os, const vector<vector<int>>& M)
{
  for (auto Mi : M) os << Mi << "\n";
  return os;
}

// Given alphas (and pluspairs, minuspairs, fours), sigmas, faces,
// return the invariants of H_1 as a Z-module in either the SL2 or GL2
// cases or both

// group=1 for GL2 only, 2 for SL2 only, 3 for both

vector<vector<INT>> integral_homology(const vector<CuspList>& faces,
                                      const CuspList& alphas, const CuspList& sigmas,
                                      const vector<vector<Quad>>& pluspairs,
                                      const vector<vector<Quad>>& minuspairs,
                                      const vector<vector<Quad>>& fours,
                                      int group, int debug)
{
  if (debug)
    {
      cout<<"alphas: "<<alphas<<endl;
      cout<<"sigmas: "<<sigmas<<endl;
    }
  vector<vector<int>> M10 = edge_boundary_matrix(alphas, sigmas);
  if (debug)
    {
      cout << "edge boundary matrix M10 has size " << M10.size() << " x " << M10[0].size() << endl;
      if (debug>1) cout << "M10 = \n" << M10 << endl;
    }
  vector<vector<int>> M21G, M21S;
  if (group&1)
    {
      M21G = face_boundary_matrix(faces, alphas, sigmas, pluspairs, minuspairs, fours, 1);
      if (debug)
        {
          cout << "GL2 face boundary matrix M21 has size " << M21G.size() << " x " << M21G[0].size() << endl;
        }
    }
  if (group&2)
    {
      M21S = face_boundary_matrix(faces, alphas, sigmas, pluspairs, minuspairs, fours, 0);
      if (debug)
        {
          cout << "SL2 face boundary matrix M21 has size " << M21S.size() << " x " << M21S[0].size() << endl;
        }
    }
  vector<vector<INT>> invs;
  if (group&1)
    invs.push_back(homology_invariants(M10, M21G, debug>1));
  if (group&2)
    invs.push_back(homology_invariants(M10, M21S, debug>1));
  return invs;
}

// Return the edge boundary matrix M10 (matrix of delta: 1-chains -> 0-chains)

// The edge basis consists of first the [a,oo] for a in alphas (whose
// boundary is trivial), then the [s,oo] for s in sigmas[1:].  The
// cusp basis is indexed by ideal classes. Same for SL2 and GL2.
vector<vector<int>> edge_boundary_matrix(const CuspList& alphas, const CuspList& sigmas)
{
  int na = alphas.size(), ns = sigmas.size()-1; // omit oo
  int h = Quad::class_number;
  // Each row as  h columns, and there are na+ns rows
  vector<vector<int>> M;
  M.reserve(na+ns);
  for (int i=0; i<na; i++) // first na rows are 0 (boundary of principal edge)
    {
      vector<int> row(h, 0);
      M.push_back(row);
    }
  for (int i=0; i<ns; i++) // next ns rows are boundary of edge {oo,s}
    {
      vector<int> row(h, 0);
      row[0] = -1;
      row[sigmas[i+1].ideal_class()] = +1;
      M.push_back(row);
    }
  return M;
}

// Return the index of an edge {a,b} in the range 0..#alphas+#sigmas-2
// For e={a,b} with a,b principal return i (>=0), where [oo,alphas[i]]
// is SL(2,Ok)-equivalent to e.

// If a is principal and b singular, return i+len(alphas)-1 (>=1),
// where [oo,sigmas[i]] is SL(2,Ok)-equivalent to e.

// If a is singular and b principal, return -(i+len(alphas)-1) (<0),
// where [sigmas[i],oo] is SL(2,Ok)-equivalent to e.

// Raise an error if both a and b are singular.

// Note that we assume sigmas[0]=oo which is not a singular point;
// the number of singular points is len(sigmas)-1.

int edge_index(const EDGE& e, const CuspList& alphas, const CuspList& sigmas)
{
  RatQuad a = e.alpha(), b = e.beta();
  Quad temp;
  int i;

  auto is_finite_singular = [alphas, sigmas, &temp](const RatQuad& c)
  {
    int j = cusp_index_with_translation(c, sigmas, temp);
    return (j>0? j : 0);
  };

  if (is_finite_singular(b))
    {
      assert (!is_finite_singular(a));
      mat22 M = Malpha(a);  // M(a)=oo
      i = cusp_index_with_translation(M(b), sigmas, temp);
      assert (i>0);
      return alphas.size() + i-1;
    }
  // now b is principal
  if (is_finite_singular(a))
    {
      return - edge_index(EDGE(b,a), alphas, sigmas);
    }
  // now a and b are principal
  mat22 M = Malpha(a);  // M(a)=oo
  i = cusp_index_with_translation(M(b), alphas, temp);
  assert (i>=0);
  return i;
}

// Return the image under delta of the face, as a vector of length #alphas+#sigmas-1
vector<int> face_boundary_vector(const CuspList& face, const CuspList& alphas, const CuspList& sigmas)
{
  vector<int> v(alphas.size()+sigmas.size()-1,0);
  int n=face.size();
  for (int i=0; i<n; i++)
    {
      int j = edge_index(EDGE(face[i], face[(i+1)%n]), alphas, sigmas);
      if (j<0) v[-j]-=1; else v[j]+=1;
    }
  return v;
}

// Return the first half of the face boundary matrix M21 (matrix of
// delta: 2-chains -> 1-chains), giving the edge identifications.
// This depends on whether we want to GL2 or SL2 homology, as the edge
// identifications change.  For the edge identifications we need
// either (first version) the pairing data in plus_pairs, minus_pairs
// and fours, or the same encoded in globals edge_pairs_plus,
// edge_pairs_minus, edge_fours

// This block of rows encodes edge identifications.  The number of
// columns is #alphas+#sigmas-1.

vector<vector<int>> edge_pairings(const CuspList& alphas, const CuspList& sigmas,
                                        const vector<vector<Quad>>& pluspairs,
                                        const vector<vector<Quad>>& minuspairs,
                                        const vector<vector<Quad>>& fours,
                                        int GL2)
{
  int
    nalphas = alphas.size(), nsigmas = sigmas.size()-1,
    nplus = pluspairs.size(), nminus = minuspairs.size(), nfours = fours.size();
  int
    ncols = nalphas + nsigmas, // size of edge basis
    nrows = nplus + nminus + nfours + (GL2? nalphas + nsigmas : nfours);
  vector<vector<int>> M;
  M.reserve(nrows);

  int i, j, n = 0; // will count number of rows pushed onto M
  Quad temp, s, r, r1, r2;

  // edge identifications
  // (0) GL2 only: {a,oo}={-a,oo} for a in alphas and {s,oo}={-s,oo} for s in sigmas
  // (1) {a,oo}+{-a,oo}=0     for a=r/s, (r,s) in pluspairs
  // (2) {a,oo}=0            for a=r/s, (r,s) in minuspairs
  // (3) {a1,oo}+{a2,oo}=0    for a1=r1/s, a2=r2/s, (s,r1,r2) in fours
  //   and if not GL2: {-a1,oo}+{-a2,oo}=0 for a1=r1/s, a2=r2/s, (s,r1,r2) in fours

  if (GL2) // type (0)
    {
      for (i=0; i<nalphas; i++)
        {
          vector<int> row(ncols, 0);
          j = cusp_index_with_translation(-alphas[i], alphas, temp);
          assert ((j>=0)&&(j<nalphas));
          row[i] +=1;
          row[j] -=1;
          // cout<<"i="<<i<<": alphas[i]="<<alphas[i]<<"; j="<<j<<": alphas[j]="<<alphas[j]<<" so row is "<<row<<endl;
          M.push_back(row);
          n++;
        }
      for (i=0; i<nsigmas; i++)
        {
          vector<int> row(ncols, 0);
          j = cusp_index_with_translation(-sigmas[i+1], sigmas, temp) -1;
          assert ((j>=0)&&(j<nsigmas));
          row[nalphas+i] +=1;
          row[nalphas+j] -=1;
          M.push_back(row);
          n++;
        }
    }
  // type (1)
  for (i=0; i<nplus; i++)
    {
      vector<int> row(ncols, 0);
      Quad ri = pluspairs[i][0], si = pluspairs[i][1];
      j = cusp_index_with_translation(RatQuad(ri,si), alphas, temp);
      assert ((j>=0)&&(j<nalphas));
      row[j] +=1;
      j = cusp_index_with_translation(RatQuad(-ri,si), alphas, temp);
      assert ((j>=0)&&(j<nalphas));
      row[j] +=1;
      // cout<<"pluspairs["<<i<<"]="<<pluspairs[i]<<" so row is "<<row<<endl;
      M.push_back(row);
      n++;
    }
  // type (2)
  for (i=0; i<nminus; i++)
    {
      vector<int> row(ncols, 0);
      Quad ri = minuspairs[i][0], si = minuspairs[i][1];
      j = cusp_index_with_translation(RatQuad(ri,si), alphas, temp);
      assert ((j>=0)&&(j<nalphas));
      // row[j] +=2; // no: if g(e)=-e then e=0, not 2e=0
      row[j] +=1;
      // cout<<"minuspairs["<<i<<"]="<<minuspairs[i]<<" so row is "<<row<<endl;
      M.push_back(row);
      n++;
    }
  // types (3) and (3')
  for (i=0; i<nfours; i++)
    {
      vector<int> row(ncols, 0);
      Quad si = fours[i][0], r1i = fours[i][1], r2i = fours[i][2];
      j = cusp_index_with_translation(RatQuad(r1i,si), alphas, temp);
      assert ((j>=0)&&(j<nalphas));
      row[j] +=1;
      j = cusp_index_with_translation(RatQuad(r2i,si), alphas, temp);
      assert ((j>=0)&&(j<nalphas));
      row[j] +=1;
      // cout<<"fours["<<i<<"]="<<fours[i]<<" so row is "<<row<<endl;
      M.push_back(row);
      n++;
      if (!GL2)
        {
          vector<int> row2(ncols, 0);
          j = cusp_index_with_translation(RatQuad(-r1i,si), alphas, temp);
          assert ((j>=0)&&(j<nalphas));
          row2[j] +=1;
          j = cusp_index_with_translation(RatQuad(-r2i,si), alphas, temp);
          assert ((j>=0)&&(j<nalphas));
          row2[j] +=1;
          M.push_back(row2);
          n++;
        }
    }
  assert (n==nrows && n==(int)M.size());
  return M;
}

vector<vector<int>> face_boundaries(const vector<CuspList>& faces,
                                    const CuspList& alphas, const CuspList& sigmas,
                                    int GL2, int debug)
{
  if (debug) cout<<"In face_boundaries("<<GL2<<")"<<endl;
  int nrows = faces.size() * (GL2? 1 : 2);
  if (debug) cout<<"nrows = "<<nrows<<endl;
  vector<vector<int>> M;
  M.reserve(nrows);
  int n=0;
  for (const auto& face : faces)
    {
      M.push_back(face_boundary_vector(face, alphas, sigmas));
      n++;
      if (!GL2)
        {
          M.push_back(face_boundary_vector(negate_polygon(face), alphas, sigmas));
          n++;
        }
    }
  if (debug) cout << "Matrix has size "<<faces.size()<<" x "<< M[0].size() << ":\n" << M << endl;
  assert (n==nrows);
  return M;
}

// Return the face boundary matrix M21 (matrix of delta: 2-chains ->
// 1-chains).  This depends on whether we want to GL2 or SL2 homology,
// as the edge identifications change.  For the edge identifications
// we need to pairing data in plus_pairs, minus_pairs and fours.

// The first block of rows encodes edge identifications, the second
// block face boundaries.  The number of columns is #alphas+#sigmas-1.
// For both GL2 and SL2 we have #faces+#plus_pairs+#minus_pairs+#fours
// rows; then for GL2 we have an extra #alphas+#sigmas-1 rows, while
// for SL2 we have an extra #faces+#fours rows.

vector<vector<int>> face_boundary_matrix(const vector<CuspList>& faces,
                                         const CuspList& alphas, const CuspList& sigmas,
                                         const vector<vector<Quad>>& pluspairs,
                                         const vector<vector<Quad>>& minuspairs,
                                         const vector<vector<Quad>>& fours,
                                         int GL2)
{
  vector<vector<int>> M = edge_pairings(alphas, sigmas, pluspairs, minuspairs, fours, GL2);
  // cout << "edge pairing matrix has size " << M.size() << " x " << M[0].size() << endl;
  // cout << M << endl;
  vector<vector<int>> M2 = face_boundaries(faces, alphas, sigmas, GL2);
  // cout << "face boundaries matrix has size " << M2.size() << " x " << M2[0].size() << endl;
  // cout << M2 << endl;
  M.insert(M.end(), M2.begin(), M2.end());
  return M;
}

int is_product_zero(const vector<vector<int>>& M1, const vector<vector<int>>& M2)
{
  // M1 represents a n1xn0 matrix and M2 a n2xn1, with M2*M1=0
  long n0 = M1[0].size(), n1 = M1.size(), n2 = M2.size();
  assert (n1==(long)M2[0].size());
  // convert from vector<vector<int>> to fmpz_mats:
  fmpz_mat_t A10, A21, Z;
  fmpz_mat_init(A10, n1, n0);
  fmpz_mat_init(A21, n2, n1);
  fmpz_mat_init(Z, n2, n0);

  make_mat(A10, M1); // size n1xn0
  make_mat(A21, M2); // size n2xn1

  // Check that A21*A10=0:
  fmpz_mat_init(Z, n2, n0);
  fmpz_mat_mul(Z, A21, A10);

  int res = fmpz_mat_is_zero(Z);
  fmpz_mat_clear(Z);
  fmpz_mat_clear(A10);
  fmpz_mat_clear(A21);
  return res;
}

// Given integer matrices (encoded as vector<vector<int>>) of the boundary maps
// M10: 1-chains -> 0-chains (as from edge_boundary_matrix())
// M21: 2-chains -> 1-chains (as from face_boundary_matrix())
// return the invariants of the integral 1-homology

// NB Both matrices are formed by rows, and act on row-vectors on the right

vector<INT> homology_invariants(const vector<vector<int>>& M10, const vector<vector<int>>& M21, int debug)
{
  vector<INT> invs =
#ifdef INVARIANTS_VIA_FLINT
    homology_invariants_via_flint(M10, M21, debug);
#else
    homology_invariants_via_pari(M10, M21, debug);
#ifdef CHECK_INVARIANTS
    vector<INT> invs2 = homology_invariants_via_flint(M10, M21, debug);
    assert (invs==invs2);
#endif
#endif
  return invs;
}

void show_matrix_data(const vector<vector<int>>& M)
{
  long nrows = M.size(), ncols = M[0].size();
  cout << "Matrix with nrows="<<nrows<<", ncols="<<ncols<<endl;
  map<int,int> entry_dist;
  for (auto Mi : M)
    for (auto Mij : Mi)
      entry_dist[Mij] +=1;
  cout << "Distribution of entries:\n";
  for (auto x : entry_dist)
    {
      cout<<x.first<<":\t"<<x.second<<endl;
    }
}

void show_invariants(const vector<INT>& v, int pretty)
{
  map<INT,int> mults;
  for (const auto& vi : v)
    mults[vi]++;
  if (pretty)
    {
      int i=0;
      for (const auto& mi : mults)
        {
          INT d=mi.first;
          int e=mi.second;
          if (i)
            cout<<" x ";
          if (e>1 && d!=0) cout << "(";
          if (d>0)
            cout << "Z/"<<d<<"Z";
          else
            cout << "Z";
          if (e>1)
            {
              if (d!=0) cout << ")";
              cout << "^" << e;
            }
          i++;
        }
      if (i==0)
        cout << "trivial";
    }
  else
    {
      cout << mults[INT(0)] << " (";
      int i=0;
      for (const auto& mi : mults)
        {
          INT d=mi.first;
          if (d==0)
            continue;
          int e=mi.second;
          if (i) cout << ",";
          cout << d;
          if (e>1) cout << "^" << e;
          i++;
        }
      cout << ")";
    }
}
