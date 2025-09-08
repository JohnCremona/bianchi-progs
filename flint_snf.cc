#include "eclib.h"
#include "flint_snf.h"

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
        fmpz_set_si(fmpz_mat_entry(A, i, j), row[j]);
    }
}

void make_mat( fmpz_mat_t A, const vector<vector<INT>>& M)
{
  long nrows = M.size(), ncols = M[0].size();
  for (long i=0; i<nrows; i++)
    {
      const vector<INT>& row = M[i];
      for (long j=0; j<ncols; j++)
        fmpz_set(fmpz_mat_entry(A, i, j), row[j].z);
    }
}

// Inversely, given a FLINT fmpz_mat, construct a matrix as vector<vector<int>>

void unmake_mat( fmpz_mat_t A, vector<vector<int>>& M)
{
  long nrows = fmpz_mat_nrows(A), ncols = fmpz_mat_ncols(A);
  for (long i=0; i<nrows; i++)
    {
      vector<int> row(ncols);
      for (long j=0; j<ncols; j++)
        row[j] = fmpz_get_si(fmpz_mat_entry(A, i, j));
      M.push_back(row);
    }
}

void unmake_mat( fmpz_mat_t A, vector<vector<INT>>& M)
{
  long nrows = fmpz_mat_nrows(A), ncols = fmpz_mat_ncols(A);
  for (long i=0; i<nrows; i++)
    {
      vector<INT> row(ncols);
      for (long j=0; j<ncols; j++)
        row[j] = INT(fmpz_mat_entry(A, i, j));
      M.push_back(row);
    }
}

long rank(const vector<vector<int>>& M)
{
  fmpz_mat_t A;
  fmpz_mat_init(A, M.size(), M[0].size());
  make_mat(A, M);
  long r = fmpz_mat_rank(A);
  fmpz_mat_clear(A);
  return r;
}

long rank(const vector<vector<INT>>& M)
{
  fmpz_mat_t A;
  fmpz_mat_init(A, M.size(), M[0].size());
  make_mat(A, M);
  long r = fmpz_mat_rank(A);
  fmpz_mat_clear(A);
  return r;
}

vector<vector<int>> HNF(const vector<vector<int>>& M)
{
  fmpz_mat_t A;
  fmpz_mat_init(A, M.size(), M[0].size());
  make_mat(A, M);

  cout << "Before HNF:\n";
  fmpz_mat_print_pretty(A);
  cout<<endl;

  fmpz_mat_hnf(A, A);

  cout << "After HNF:\n";
  fmpz_mat_print_pretty(A);
  cout<<endl;

  vector<vector<int>> H;
  unmake_mat(A, H);
  fmpz_mat_clear(A);
  return H;
}

vector<vector<INT>> HNF(const vector<vector<INT>>& M)
{
  fmpz_mat_t A, At;
  fmpz_mat_init(A, M.size(), M[0].size());
  fmpz_mat_init(At, M[0].size(), M.size());

  make_mat(A, M);

  // cout << "Before HNF:\n";
  // fmpz_mat_print_pretty(A);
  // cout<<endl;

  fmpz_mat_transpose(At, A);
  fmpz_mat_hnf(At, At); // transpose and hnf
  fmpz_mat_transpose(A, At);

  // cout << "After HNF:\n";
  // fmpz_mat_print_pretty(A);
  // cout<<endl;

  vector<vector<INT>> H;
  unmake_mat(A, H);
  fmpz_mat_clear(A);
  fmpz_mat_clear(At);
  return H;
}

// Return a list of the pivotal columns of the HNF of a matrix
// (encoded as vector<vector<int>>) for which the pivots are =1
vector<int> HNF_pivots(const vector<vector<int>>& M)
{
  vector<int> ans;
  auto H = HNF(M);
  for (const auto& row : H)
    {
      // find first nonzero entry (if any)
      auto search = std::find_if(row.begin(), row.end(), [](int x){return (x!=0);});
      if (search==row.end())
        continue;
      if (*search >1)
        continue;
      int j = search-row.begin();
      // cout << "pivot=1 in column "<<j<<" in row "<<row<<endl;
      ans.push_back(j);
    }
  return ans;
}

void SNF(fmpz_mat_t& S, fmpz_mat_t& A)
{
  long nrows = fmpz_mat_nrows(A), ncols = fmpz_mat_ncols(A);
  // cout<<"In SNF(A) with "<<nrows<<" rows and "<<ncols<<" columns"<<endl;
  fmpz_mat_t H, H0, Ht;
  fmpz_mat_init(H, nrows, ncols);
  fmpz_mat_init(H0, nrows, ncols);
  fmpz_mat_init(Ht, ncols, nrows);
  fmpz_mat_set(H, A);
  // fmpz_mat_hnf(H, H);
  // cout<<" - initial HNF step done"<<endl;
  int ok = 0;
  int nsteps = 0;
  while (!ok)
    {
      nsteps++;
      fmpz_mat_set(H0, H);  // keep current H in H0
      fmpz_mat_transpose(Ht, H);
      fmpz_mat_hnf(Ht, Ht); // transpose and hnf
      fmpz_mat_transpose(H, Ht);
      fmpz_mat_hnf(H, H);   // transpose back and hnf again
      ok = fmpz_mat_equal(H, H0); // see if anything changed
      // cout<<" - step "<<nsteps<<" done"<<endl;
    }
  // cout<<"stabilised after "<<nsteps<<" steps. Now doing final snf.\n";
  fmpz_mat_snf(S, H);
  // cout<<"SNF finished"<<endl;
  fmpz_mat_clear(H);
  fmpz_mat_clear(H0);
  fmpz_mat_clear(Ht);
}

vector<INT> homology_invariants_via_flint(const vector<vector<int>>& M10, const vector<vector<int>>& M21, int debug)
{
  // M10 represents a n1xn0 matrix and M21 a n2xn1, with M21*M10=0
  long n0 = M10[0].size(), n1 = M10.size(), n2 = M21.size();
  assert (n1==(long)M21[0].size());

  // convert from vector<vector<int>> to fmpz_mats:
  fmpz_mat_t A10, A21, Z, U, H, M;
  fmpz_mat_init(A10, n1, n0);
  fmpz_mat_init(A21, n2, n1);
  make_mat(A10, M10); // size n1xn0
  make_mat(A21, M21); // size n2xn1
  if (debug>2)
    {
      cout << "M10 as a FLINT matrix:\n";
      fmpz_mat_print_pretty(A10);
      cout<<endl;
      cout << "M21 as a FLINT matrix:\n";
      fmpz_mat_print_pretty(A21);
      cout<<endl;
    }

  // Check that A21*A10=0:
  fmpz_mat_init(Z, n2, n0);
  fmpz_mat_mul(Z, A21, A10);
  assert (fmpz_mat_is_zero(Z));

  // find H = HNF of A10, with transform matrix U s.t. H=U*A10:
  if (debug)
    {
      cout<<"Computing HNF of A10 which has size "<<n1<<" x "<<n0<<endl;
    }
  fmpz_mat_init(H, n1, n0);
  fmpz_mat_init(U, n1, n1);
  fmpz_mat_hnf_transform(H, U, A10);
  long r = fmpz_mat_rank(H);
  if (debug)
    {
      cout<<"H has size "<<n1<<" x "<<n0<<", and r="<<r<<endl;
      cout<<"U has size "<<n1<<" x "<<n1<<", now inverting U..."<<endl;
    }
  // replace U by U^-1:
  fmpz_t den;
  fmpz_init(den);
  int ok = fmpz_mat_inv(U, den, U);
  assert (ok); // means U was invertible over Q, i.e. det(U) nonzero
  if (debug)
    {
      cout<<"...finished inverting U";
      if (debug>1)
        {
          cout<<"; denom is ";
          fmpz_print(den);
        }
      cout<<endl;
    }
  assert (fmpz_equal_si(den, 1) || fmpz_equal_si(den, -1)); // denominator is +-1

  // multiply A21 by U^-1:
  if (debug)
    {
      cout<<"Multiplying A21 by U^-1...";
    }
  fmpz_mat_mul(A21, A21, U);
  if (debug)
    {
      cout<<"...done, A21*U^{-1} has size "<<n2<<" x "<<n1<<", now dropping first r columns..."<<endl;
      if (debug>2)
        {
          fmpz_mat_print_pretty(A21);
          cout<<endl;
        }
    }
  // drop first r cols of this:
  fmpz_mat_window_init(M, A21, 0, r, n2, n1);
  int homrank = n1 - r - fmpz_mat_rank(M); // ==nullity(M)
  if (debug>1)
    {
      cout<<"The window has size "<<fmpz_mat_nrows(M)<<" x "<<fmpz_mat_ncols(M)<<endl;
      if (debug>2)
        {
          fmpz_mat_print_pretty(M);
          cout<<endl;
        }
    }
  if (debug)
    cout << "Homology rank = " << homrank << "; ";
  assert (fmpz_mat_nrows(M)==n2);
  assert (fmpz_mat_ncols(M)==n1-r);

  fmpz_mat_t S;
  // Compute Smith Normal Form of that:
  fmpz_mat_init_set(S, M); // to set to the right size
  if (debug)
    {
      cout<<" (about to compute SNF of M with "<<fmpz_mat_nrows(M)<<" rows, "<<fmpz_mat_ncols(M)<<" columns)" <<endl;
      cout<<"=============="<<endl;
      // fmpz_mat_print(S);
      // cout<<"\n=============="<<endl;
    }

  SNF(S, M);

  // Extract the diagonal entries of S (omitting any 1s):
  long n = min(n2, n1-r);
  vector<INT> v;
  for (long i=0; i<n; i++)
    {
      INT m(fmpz_mat_entry(S, i, i));
      if (debug) cout<<" S["<<i<<","<<i<<"] =  "<<m<<endl;
      if (m!=1)
        v.push_back(m);
    }
  cout << "non-trivial invariants: "<<v<<endl;
  fmpz_mat_clear(S);
  fmpz_mat_window_clear(M);
  fmpz_mat_clear(A10);
  fmpz_mat_clear(A21);
  fmpz_mat_clear(Z);
  fmpz_mat_clear(H);
  fmpz_mat_clear(U);

  return v;
}
