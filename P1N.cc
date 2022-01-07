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

//#define DEBUG_LIFT

// compute a matrix M = [a, b; c, d] with det=1 lifting (c:d) in P^1(N)
mat22 lift_to_SL2(Qideal& N, const Quad& cc, const Quad& dd)
{
  Quad a, b, c(cc), d(dd), inv, x, y, z, h;
#ifdef DEBUG_LIFT
  cout<<"Lifting symbol (c:d)=("<<c<<":"<<d<<") mod "<<N<<" to SL2"<<endl;
#endif
  // Special cases (1): (c:1), (1:d) need no work:
  Quad one=Quad::one, zero=Quad::zero;
  if (d==one) return mat22(one,zero,c,one);
  if (c==one) return mat22(zero,-one,one,d);

#ifdef DEBUG_LIFT
  cout<<"Neither c nor d is invertible modulo "<<N<<": testing whether ideal (c,d) is principal"<<endl;
#endif
  // General case: neither c nor d is invertible.

  int ok=0;
  for (long i=1; i<N.norm() && !ok; i++)
    {
      Quad t = N.resnum(i);
      if (!N.is_coprime_to(t)) continue;
      Quad ct=N.reduce(c*t), dt=N.reduce(d*t);
      h = quadbezout(ct,dt, x, y);
      ok = (!h.is_zero());
      if (ok) // then is principal with ct*x+dt*y=h
        {
#ifdef DEBUG_LIFT
          cout<<"ideal (c,d)=("<<h<<"), success"<<endl;
#endif
          a = y;
          b = -x;
          c = ct / h;
          d = dt / h;
        }
    }

  if (ok)
    {
#ifdef DEBUG_LIFT
      cout<<" replacing c by "<<c<<" and d by "<<d<<", which are coprime"<<endl;
#endif
      assert (a*d-b*c==Quad::one);
      mat22 M(a,b,c,d);
#ifdef DEBUG_LIFT
      cout<<" returning  "<< M <<endl;
#endif
      return M;
    }

  // Test if (c,d)=(h), principal:
  h = quadbezout(c,d, x, y);
  if (!h.is_zero()) // then is principal with c*x+d*y=h, and h=1 since reduced
    {
#ifdef DEBUG_LIFT
      cout<<"ideal (c,d)=("<<h<<"), success"<<endl;
#endif
      a = y;
      b = -x;
      if (h.norm()>1) // should not happen as (c:d) was reduced
        {
          c /= h;
          d /= h;
        }
    }
  else
    {  // Now we must work harder.
#ifdef DEBUG_LIFT
      cout<<" (c,d) not principal, working harder..."<<endl;
#endif
      int t = N.is_coprime_to(c, d, x, y, 1);   // c*x+d*y = 1 mod N with y invertible
      assert (t==1);
#ifdef DEBUG_LIFT
      cout<<" c*x+d*y=1 mod N with x = "<<x<<" and y = "<<y<<endl;
#endif
      t = N.is_coprime_to(y, z);            // y*z = 1 mod N
      assert (t==1);
#ifdef DEBUG_LIFT
      cout<<" inverse of y mod N is z = "<<z<<" with y*z="<<y*z<<endl;
#endif
      z *= y;  // now z=0 (mod y), z=1 (mod N)
      a = y;
      b =  one - (x+one)*z;  // = 1 (mod y), -x (mod N)
      c = -one + (c+one)*z;  // =-1 (mod y),  c (mod N)
      d = (one+b*c)/a;     // 1+b*c = 0 (mod y), 1-c*x = d*y = a*d (mod N)
    }
#ifdef DEBUG_LIFT
  cout<<" replacing c by "<<c<<" and d by "<<d<<", which are coprime"<<endl;
#endif
  assert (a*d-b*c==one);
  mat22 M(a,b,c,d);
#ifdef DEBUG_LIFT
  cout<<" returning  "<< M <<endl;
#endif
  return M;
}

// compute a matrix M = [a, b; c, d] with det=1 lifting (c:d)
mat22 P1N::lift_to_SL2(long ind)
{
  Quad c, d;
  make_symb(ind, c, d); // this is reduced, in particular is (c:1) or (1:d) if possible
  mat22 M = ::lift_to_SL2(N, c, d);
  assert (M.is_unimodular());
  assert (index(M.entry(1,0), M.entry(1,1))==ind);
  return M;
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

void P1N::check_lifts(int verbose) // checks lifts to SL2
{
  if (verbose)
    {
      cout << "Testing lifts from P1(N) to SL(2,O_K) for N = " << N << ": ";
    }
  long i;
  Quad c, d;
  for (i=0; i<psi; i++)
    {
      make_symb(i, c, d);
      mat22 M = lift_to_SL2(i);
      if (verbose)
        {
          cout << i << " --> ("<<c<<":"<<d << ") lifts to " << M << endl;
        }
      // NB The lift_to_SL2() method already tests these assertions:
      // assert(M.is_unimodular());
      // assert (i==index(M.entry(1,0), M.entry(1,1)));
    }
}
