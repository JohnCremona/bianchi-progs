#include "arith_extras.h"
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

