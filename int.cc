// FILE INT.CC -- implementation of some of the interface functions in int.h

#include <assert.h>
#include "int.h"

#ifdef FLINT // else none of these are needed

std::ostream& operator<<(std::ostream& s, const INT& a)
{
  s << std::string(fmpz_get_str(NULL, 10, a.z));
  return s;
}

std::istream& operator>>(std::istream& s, INT& x)
{
  std::string n;
  s>>n;
  fmpz_set_str(x.z, n.c_str(), 10);
  return s;
}

int divrem(const INT& a, const INT& b, INT& quo, INT& rem)
{
  fmpz_fdiv_qr(quo.z, rem.z, a.z, b.z);
  return is_zero(rem);
}

INT rounded_division(const INT& a, const INT& b)
{
  INT q, r;
  divrem(a,b,q,r);
  INT r2 = 2*r;
  // We want -b <= r2 < +b
  if (r2<-b)
    q-=1;
  else
    if (r2>=b)
      q+=1;
  return q;
}

std::vector<INT> pdivs(const INT& a)
{
  fmpz_factor_t f;
  fmpz_factor_init(f);
  fmpz_factor(f, a.z);
  std::vector<INT> ans;
  for (int i =0; i< f->num; i++)
    {
      INT p;
      fmpz_init_set(p.z, f->p + i);
      ans.push_back(p);
    }
  return ans;
}

std::vector<INT> sqdivs(const INT& a)
{
  fmpz_factor_t f;
  fmpz_factor_init(f);
  fmpz_factor(f, a.z);
  std::vector<INT> plist;
  std::vector<int> elist;
  int nd = 1;
  for (slong i =0; i< f->num; i++)
    {
      INT p;
      fmpz_init_set(p.z, f->p + i);
      plist.push_back(p);
      int e = (f->exp[i])/2;
      elist.push_back(e);
      nd *= (1+e);
    }
  std::vector<INT> dlist(1, INT(1));
  dlist.resize(nd);
  nd = 1;
  auto pr = plist.begin();
  auto ei = elist.begin();
  while(pr!=plist.end())
    {
      INT p=*pr++;
      int e=*ei++;
      for (int j=0; j<e; j++)
        for (int k=0; k<nd; k++)
          dlist[nd*(j+1)+k] = p*dlist[nd*j+k];
      nd*=(e+1);
    }
  return dlist;
}

INT invmod(const INT&a, const INT& p)
{
  INT b;
  fmpz_invmod(b.z, a.z, p.z);
  return b;
}

long invmod(const INT&a, long p)
{
  fmpz_t b, P;
  fmpz_init_set_si(P,p);
  fmpz_invmod(b, a.z, P);
  return (long)fmpz_get_si(b);
}

long I2long(const INT& a)
{
  if (a.is_long())
    return (long)fmpz_get_si(a.z);
  assert(0 && "INT does not fint into an slong");
}

INT sqrt_mod_p(const INT& a, const INT& p)
{
  INT b;
  sqrt_mod_p(b, a, p);
  return b;
}

#endif