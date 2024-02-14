// FILE SWAN_HOM.CC: implementation of functions for computing integral 1-homology

#include "swan_utils.h"
#include "swan_tess.h"
#include "swan_hom.h"

#include <flint/fmpz_mat.h>

// Given alphas (and pluspairs, minuspairs, fours), sigmas, faces,
// return the invariants of H_1 as a Z-module in either the SL2 or GL2
// cases

vector<int> integral_homology(const vector<CuspList>& faces,
                              const CuspList& alphas, const CuspList& sigmas,
                              const vector<vector<Quad>>& pluspairs,
                              const vector<vector<Quad>>& minuspairs,
                              const vector<vector<Quad>>& fours,
                              int GL2)
{
  vector<vector<int>> M10 = edge_boundary_matrix(alphas, sigmas);
  cout<<"Computed M10 = \n";
  for (const auto& row: M10)
    cout << row <<endl;
  vector<vector<int>> M21 = face_boundary_matrix(faces, alphas, sigmas, pluspairs, minuspairs, fours, GL2);
  cout<<"Computed M21 = \n";
  for (const auto& row: M21)
    cout << row <<endl;
  return homology_invariants(M10, M21);
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
  int i, j, n=face.size();
  for (i=0; i<n; i++)
    {
      j = edge_index(EDGE(face[i], face[(i+1)%n]), alphas, sigmas);
      if (j<0) v[-j]-=1; else v[j]+=1;
    }
  return v;
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
  int nalphas = alphas.size(), nsigmas = sigmas.size()-1,
    nplus = pluspairs.size(), nminus = minuspairs.size(), nfours = fours.size(),
    nfaces = faces.size();
  int ncols = nalphas + nsigmas; // size of edge basis
  int nrows = nfaces + nplus + nminus + nfours + (GL2? nalphas + nsigmas : nfaces + nfours);
  vector<vector<int>> M;
  M.reserve(nrows);

  cout<<"nalphas="<<nalphas<<endl;
  cout<<"nsigmas="<<nsigmas<<endl;
  cout<<"nplus="<<nplus<<endl;
  cout<<"nminus="<<nminus<<endl;
  cout<<"nfours="<<nfours<<endl;
  cout<<"nfaces="<<nfaces<<endl;

  cout<<"M should have "<<nrows<<" rows"<<endl;

  int i, j;
  int n = 0; // will count number of rows pushed onto M
  Quad temp, s, r, r1, r2;

  // edge identifications
  // (0) GL2 only: {a,oo}={-a,oo} for a in alphas and {s,oo}={-s,oo} for s in sigmas
  // (1) {a,oo}+{-a,oo}=0     for a=r/s, (r,s) in pluspairs
  // (2) 2{a,oo}=0            for a=r/s, (r,s) in minuspairs
  // (3) {a1,oo}+{a2,oo}=0    for a1=r1/s, a2=r2/s, (s,r1,r2) in fours
  //   and if not GP2: {-a1,oo}+{-a2,oo}=0 for a1=r1/s, a2=r2/s, (s,r1,r2) in fours

  if (GL2) // type (0)
    {
      for (i=0; i<nalphas; i++)
        {
          vector<int> row(ncols, 0);
          j = cusp_index_with_translation(-alphas[i], alphas, temp);
          row[i] +=1;
          row[j] -=1;
          M.push_back(row);
          n++;
        }
      for (i=0; i<nsigmas; i++)
        {
          vector<int> row(ncols, 0);
          j = cusp_index_with_translation(-sigmas[i+1], sigmas, temp) -1;
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
      Quad r = pluspairs[i][0], s = pluspairs[i][1];
      j = cusp_index_with_translation(RatQuad(r,s), alphas, temp);
      row[j] +=1;
      j = cusp_index_with_translation(RatQuad(-r,s), alphas, temp);
      row[j] +=1;
      M.push_back(row);
      n++;
    }
  // type (2)
  for (i=0; i<nminus; i++)
    {
      vector<int> row(ncols, 0);
      Quad r = minuspairs[i][0], s = minuspairs[i][1];
      j = cusp_index_with_translation(RatQuad(r,s), alphas, temp);
      row[j] +=2;
      M.push_back(row);
      n++;
    }
  // types (3) and (3')
  for (i=0; i<nfours; i++)
    {
      vector<int> row(ncols, 0);
      Quad s = fours[i][0], r1 = fours[i][1], r2 = fours[i][2];
      j = cusp_index_with_translation(RatQuad(r1,s), alphas, temp);
      row[j] +=1;
      j = cusp_index_with_translation(RatQuad(r2,s), alphas, temp);
      row[j] +=1;
      M.push_back(row);
      n++;
      if (!GL2)
        {
          vector<int> row2(ncols, 0);
          j = cusp_index_with_translation(RatQuad(-r1,s), alphas, temp);
          row2[j] +=1;
          j = cusp_index_with_translation(RatQuad(-r2,s), alphas, temp);
          row2[j] +=1;
          M.push_back(row2);
          n++;
        }
    }

  // face boundary relations

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
  assert (n==nrows);
  return M;
}

// Given a matrix as vector<vector<int>>, construct a FLINT fmpz_mat

// NB The caller must have called call fmpz_mat_init() on the supplied
// fmpz_mat_t, and must call fmpz_mat_clear() when finished with it.

void make_mat( fmpz_mat_t A, const vector<vector<int>>& M)
{
  long nrows = M.size(), ncols = M[0].size();
  for (long i=0; i<nrows; i++)
    {
      const vector<int>& row = M[i];
      for (long j=0; j<ncols; j++)
        {
          fmpz_set_si(fmpz_mat_entry(A, i, j), row[j]);
        }
    }
}

// Given integer matrices (encoded as vector<vector<int>>) of the boundary maps
// M10: 1-chains -> 0-chains (as from edge_boundary_matrix())
// M21: 2-chains -> 1-chains (as from face_boundary_matrix())
// return the invariants of the integral 1-homology

// NB Both matrices are formed by rows, and act on row-vectors on the right

vector<int> homology_invariants(const vector<vector<int>>& M10, const vector<vector<int>>& M21)
{
  // M10 represents a n1xn0 matrix and M21 a n2xn1, with M21*M10=0
  long n0 = M10[0].size(), n1 = M10.size(), n2 = M21.size();
  assert (n1==(long)M21[0].size());

  // convert from vector<vector<int>> to fmpz_mats:
  fmpz_mat_t A10, A21, Z, U, H, M, S;
  fmpz_mat_init(A10, n1, n0);
  fmpz_mat_init(A21, n2, n1);
  make_mat(A10, M10); // size n1xn0
  make_mat(A21, M21); // size n2xn1
  cout << "M10 as a FLINT matrix:\n";
  fmpz_mat_print_pretty(A10);
  cout<<endl;
  cout << "M21 as a FLINT matrix:\n";
  fmpz_mat_print_pretty(A21);
  cout<<endl;

  // Check that A21*A10=0:
  fmpz_mat_init(Z, n2, n0);
  fmpz_mat_mul(Z, A21, A10);
  assert (fmpz_mat_is_zero(Z));

  // find H = HNF of A10, with transform matrix U s.t. H=U*A10:
  fmpz_mat_init(H, n1, n2);
  fmpz_mat_init(U, n1, n1);
  fmpz_mat_hnf_transform(H, U, A10);
  long r = fmpz_mat_rank(H);

  // replace U by U^-1:
  fmpz_t den;
  fmpz_init(den);
  int ok = fmpz_mat_inv(U, den, U);
  assert (ok); // means U was invertible over Q, i.e. det(U) nonzero
  cout<<"After inverting U, denom is ";
  fmpz_print(den);
  cout<<endl;
  assert (fmpz_equal_si(den, 1) || fmpz_equal_si(den, -1)); // denominator is +-1

  // multiply A21 by U^-1:
  fmpz_mat_mul(A21, A21, U);

  // drop first r rows of this:
  cout<<"A21 has size "<<n2<<" x "<<n1<<", and r="<<r<<endl;
  fmpz_mat_window_init(M, A21, 0, r, n2, n1);
  cout<<"The window has size "<<fmpz_mat_nrows(M)<<" x "<<fmpz_mat_ncols(M)<<endl;
  assert (fmpz_mat_nrows(M)==n2);
  assert (fmpz_mat_ncols(M)==n1-r);

  // Compute Smith Normal Form of that:
  fmpz_mat_init_set(S, M); // to set to the right size
  fmpz_mat_snf(S, M);

  // Extract the diagonal entries of S:
  long n = min(n1, n2-r);
  vector<int> v(n);
  for (long i=0; i<n; i++)
    v[i] = fmpz_get_si(fmpz_mat_entry(S, i, i));

  fmpz_mat_clear(A10);
  fmpz_mat_clear(A21);
  fmpz_mat_clear(Z);
  fmpz_mat_clear(H);
  fmpz_mat_clear(M);
  fmpz_mat_clear(S);
  fmpz_mat_clear(U);

  return v;
}

