// FILE P1N.CC: implementation of class for P^1(O_K/N) for an arbitrary ideal N

#include <iostream>

#include "P1N.h"

// utilities for a standard bijection between [0,1,...,n-1] and the
// product of [0,1,...,n_i-1] where n is the product of the n_i

long merge_indices(const vector<long>& nlist, const vector<long>& klist)
// return ((k1*n2+k2)*n3+k3)*n4+k4 (etc.)
// NB n1 is not used except implicitly
{
  long tot=0;
  vector<long>::const_iterator n=nlist.begin()+1, k=klist.begin();
  while (n!=nlist.end())
    {
      tot += *k++;
      tot *= *n++;
    }
  tot += *k;
  return tot;
}

vector<long> split_indices(const vector<long>& nlist, long k)
{
  long tot=k;
  vector<long>::const_reverse_iterator n=nlist.rbegin();
  vector<long> klist;
  while (n!=nlist.rend())
    {
      std::ldiv_t qr = ldiv(tot, *n++);
      tot = qr.quot;
      klist.push_back(qr.rem);
    }
  std::reverse(klist.begin(), klist.end());
  return klist;
}

P1N::P1N(const Qideal& I)
  : N(I)
{
  Factorization F = N.factorization(); //  factorization is cached in N
  np = F.size();
  nrm = I2long(N.norm());
  phi = psi = 1;
  for (int i=0; i<np; i++)
    {
      long normp = I2long(F.prime(i).norm());
      long e = F.exponent(i);
      phi *= (normp-1);
      psi *= (normp+1);
      e--;
      while (e--)
        {
          phi *= normp;
          psi *= normp;
        }
    }
  switch (np) {
  case 0:
    {
      break; // nothing left to do
    }
  case 1:
    {
      residue_codes.reserve(nrm);
      noninvertible_residue_indices.reserve(nrm-phi);
      for(long i=0; i<nrm; i++)
        {
          Quad r = N.resnum(i), s;
          if (N.is_coprime_to(r, s)) // r invertible with r*s = 1 (N)
            {
              residue_codes.push_back(N.numres(s));
            }
          else // r not invertible
            {
              residue_codes.push_back(-noninvertible_residue_indices.size());
              noninvertible_residue_indices.push_back(i);
            }
        }
      // cout << "P1N constructor, prime power case" << endl;
      // cout << "Residues: " << N.residues() << endl;
      // cout << "residue_codes: " << residue_codes << endl;
      // cout << "noninvertible_residue_indices: " << noninvertible_residue_indices << endl;
      break;
    }
  default: // not a prime power
    {
      P1PP.reserve(np);
      psilist.reserve(np);
      for (int i=0; i<np; i++)
        {
          P1N P1PPi(F.prime_power(i));
          P1PP.push_back(P1PPi);
          psilist.push_back(P1PPi.psi);
        }
    }
  }
}

void P1N::make_symb(long i, Quad& c, Quad& d) // the i'th (c:d) symbol
{
  if (np==0) // at level (1), the only symbol is (0:1), lifting to the identity matrix
    {
      c = Quad::zero;
      d = Quad::one;
      return;
    }
  if (np==1)
    {
      assert (i>=0);
      assert (i<psi);
      if (i<nrm)
        {
          c = N.resnum(i);
          d = Quad::one;
        }
      else
        {
          c = Quad::one;
          d = N.resnum(noninvertible_residue_indices[i-nrm]);
        }
      return;
    }
  // not a prime power, use CRT
  vector<long> ilist = split_indices(i);
  vector<Quad> clist, dlist;
  clist.reserve(np);
  dlist.reserve(np);
  Quad ck, dk;
  for (int k=0; k<np; k++)
    {
      P1PP[k].make_symb(ilist[k], ck, dk);
      clist.push_back(ck);
      dlist.push_back(dk);
    }
  c = N.factorization().solve_CRT(clist);
  d = N.factorization().solve_CRT(dlist);
  reduce(c, d);
  return;
}

void P1N::reduce(Quad& c, Quad& d)   // simplify c,d without changing (c:d)
{
  // cout<<"Reducing (c:d)=("<<c<<":"<<d<<") mod "<<N<<"..."<<endl;
  Quad u;
  if (N.is_coprime_to(d,u)) // scale to (c:1)
    {
      c = N.reduce(c*u);
      d = Quad::one;
      // cout<<"  -- returning ("<<c<<":"<<d<<")"<<endl;
      return;
    }
  if (N.is_coprime_to(c,u)) // scale to (1:d)
    {
      d = N.reduce(d*u);
      c = Quad::one;
      // cout<<"  -- returning ("<<c<<":"<<d<<")"<<endl;
      return;
    }

  // divide by gcd(c,d) if it is principal
  u = quadgcd(c,d); // if (c,d) is principal this is a generator, else 0
  if (u.norm()>1)
    {
      c /= u;
      d /= u;
      // cout<<"  -- after dividing by "<<u<<":  ("<<c<<":"<<d<<")"<<endl;
    }
  else // at least divide out by the gcd of the contents
    {
      QUINT g = gcd(c.content(), d.content());
      if (g>1)
        {
          c /= g;
          d /= g;
          // cout<<"  -- after dividing by "<<u<<":  ("<<c<<":"<<d<<")"<<endl;
        }
    }
  return;
}

long P1N::index(const Quad& c, const Quad& d) // index i of (c:d)
{
  if (np==0) return 0;
  if (np==1)
    {
      // cout<<"Finding index of ("<<c<<":"<<d<<") mod "<<N<<"..."<<endl;
      long t = residue_codes[N.numres(d)];
      if (t>0) // d is invertible
        return N.numres(c*N.resnum(t));
      // now d is not invertible so c must be
      t = residue_codes[N.numres(c)];
      assert (t>0);
      t = residue_codes[N.numres(d*N.resnum(t))];
      assert (t<=0);
      return nrm - t;
    }
  // general case, np>1 so not a prime power
  vector<long> ilist;
  ilist.reserve(np);
  for (int i=0; i<np; i++)
    ilist.push_back(P1PP[i].index(c,d));
  return merge_indices(ilist);
}

// each M in GL2 permutes the (c:d) symbols by right multiplcation:
long P1N::apply(const mat22& M, long i)
{
  if (M.is_scalar()) return i;
  if (np==0) return 0;
  if (np==1)  // works in general but much simpler in this case (no CRT in constructing c,d)
    {
      Quad c, d;
      make_symb(i, c, d);
      M.apply_right(c,d);
      return index(c,d);
    }
  // general case, np>1 so not a prime power
  vector<long> jlist, ilist = split_indices(i);
  jlist.reserve(np);
  for (int k=0; k<np; k++)
    jlist.push_back(P1PP[k].apply(M,ilist[k]));
  return merge_indices(jlist);
}

// return a matrix [a, b; c, d] with det=1 and (c:d)= i'th symbol in P^1(N)
mat22 P1N::lift_to_SL2(long ind)
{
  Quad c, d;
  make_symb(ind, c, d); // this is reduced, in particular is (c:1) or (1:d) if possible
  mat22 M = ::lift_to_SL2(N, c, d);
  assert (M.is_unimodular());
  assert (index(M.entry(1,0), M.entry(1,1))==ind);
  return M;
}

// return a matrix [a, b; c, d] with det=1, c in M and (c:d) equal
// to the i'th symbol (requires M,N coprime). If (u,v)!=(0,0) they
// should satisfy u+v=1 with u in N, v in M, otherwise such u,v will
// be computed and returned.
mat22 P1N::lift_to_Gamma_0(long i, Qideal M, const Quad& u, const Quad& v)
{
  Quad uu(u), vv(v);
  if (uu.is_zero() && vv.is_zero())
    M.is_coprime_to(N, uu, vv); // u+v=1, u in M, v in N
  Quad c, d;
  make_symb(i, c, d); // this is reduced, in particular is (c:1) or (1:d) if possible
  mat22 G = ::lift_to_Gamma_0(M, N, c, d, uu, vv);
  assert (G.is_unimodular());
  assert (index(G.entry(1,0), G.entry(1,1))==i);
  assert (M.contains(G.entry(1,0)));
  return G;
}

// test functions
void P1N::check(int verbose)       // checks indexing
{
  if (verbose)
    {
      cout << "Testing P1(N) for N = " << N << ": ";
      cout << "psi(N) = "<<psi<<endl;
    }
  long i;
  Quad c, d;
  for (i=0; i<psi; i++)
    {
      make_symb(i, c, d);
      long j = index(c, d);
      if (verbose)
        {
          cout << i << " --> ("<<c<<":"<<d << ") --> " << j << endl;
          assert(i==j);
        }
      else
        {
          if (i!=j)
            {
              cout << i << " --> ("<<c<<":"<<d << ") --> " << j << endl;
            }
        }
    }
}

void P1N::check_lifts(int verbose) // checks lifts to SL2 and Gamma_0(P) for a P coprime to N
{
  long i;
  Quad c, d, u, v;

  Quadprime P; // the first prime not dividing the level
  for (vector<Quadprime>::const_iterator pr = Quadprimes::list.begin(); pr!=Quadprimes::list.end(); ++pr)
    {
      P = *pr;
      if (!P.divides(N))
        break;
    }

  if (verbose)
    {
      cout << "Testing lifts from P1(N) to SL(2,O_K) and to Gamma_0(P) for N = " << N << " and P= "<<P<<":\n";
    }
  for (i=0; i<psi; i++)
    {
      make_symb(i, c, d);
      mat22 M = lift_to_SL2(i);
      if (verbose)
        {
          cout << i << " --> ("<<c<<":"<<d << ") lifts to " << M;
        }
      M = lift_to_Gamma_0(i, P, u, v);
      if (verbose)
        {
          cout << " --> "<< M << " (in Gamma_0(P))" << endl;
        }
      // NB The lift_to_SL2() method already tests these assertions:
      // assert(M.is_unimodular());
      // assert (i==index(M.entry(1,0), M.entry(1,1)));
    }
}
