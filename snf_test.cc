// Compile with gcc -std=c++11 snf_test.cc -o snf_test -lflint -lgmp -lstdc++

#include <iostream>
#include <flint/fmpz.h>
#include <flint/fmpz_mat.h>

int main ()
{
  fmpz_mat_t M, S;
  int r, c;
  std::cin >> r >> c;
  fmpz_mat_init(M,r,c);
  fmpz_mat_read(M);
  std::cout << "Input matrix (" << r<<"x"<<c<<"):\n";
  fmpz_mat_print_pretty(M);
  std::cout << "\n";
  fmpz_mat_init_set(S, M); // to set to the right size
  fmpz_mat_hnf(M, M);
  std::cout << "HNF matrix:\n";
  fmpz_mat_print_pretty(M);
  std::cout << "\n";
  fmpz_mat_snf(S, M);
  std::cout << "SNF matrix:\n";
  fmpz_mat_print_pretty(S);
  std::cout << "\n";
  fmpz_mat_clear(M);
  fmpz_mat_clear(S);
}
