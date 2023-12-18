#include "looper.h"

// Set bmin and bmax for this value of n, and initialise the b-loop
void Quadlooper::setblims()
{
  n4 = (d%4==3? 4*n: n);
  INT zero(0);
  switch (d)
    {
    case 1:
      bmin = db2 = zero;
      bmax = include_conjugates? isqrt(n-1) : isqrt(n/2);
      break;
    case 3:
      bmin = db2 = zero;
      bmax = include_conjugates? isqrt(n-1) : isqrt(n/3);
      break;
    default:
      bmax = isqrt(n4/d);
      bmin = include_conjugates? -bmax : zero;
      db2 = d*bmin*bmin;
      if (db2==n4)
        {
          bmin+=1;
          db2 = d*bmin*bmin;
        }
    }
  // cout<<"n="<<n<<": Setting bmin="<<bmin<<", bmax="<<bmax<<endl;
  b = bmin;
  db2 = d*b*b;
  while(!testb())
    bstep();
}

// test if current b is valid, setting val if so
int Quadlooper::testb()
{
  INT a, a2 = (Quad::t? 4*n-db2 : n-db2);
  if(!isqrt(a2, a))
    return 0;
  if(Quad::t)
    a = (a-b)/2;
  val = Quad(a,b);
  return 1;
}

// increment b if b<bmax, otherwise, increment n
void Quadlooper::bstep()
{
  if (b<bmax)
    {
      b+=1;      // cout<<"Increasing b to "<<b<<endl;
      db2 = d*b*b;
    }
  else
    nstep(); // calls setblims which sets b and db2
}

// increment n if possible, skipping values certainly not norms
void Quadlooper::nstep()
{
  n+=1;
  if (!ok())
    return;
  while(kronecker(disc,n)==-1)
    n+=1;
  setblims();
}

// increment b until next valid value
void Quadlooper::operator++()
{
  bstep();
  while(!testb())
    bstep();
}


vector<Quad> Quadlooper::values_with_current_norm()
{
  // must make a copy
  INT m = n;
  return values_with_norm_up_to(m);
}

vector<Quad> Quadlooper::values_with_norm_up_to(const INT& m)
{
  vector<Quad> values;
  if (n>m)
    return values;
  values.push_back(val);
  operator++();
  while (n <= m)
    {
      values.push_back(val);
      operator++();
    }
  return values;
}


// Lists of elements of norm in ranges (up to units, excluding 0)
vector<Quad> elements_of_norm_between(const INT& n1, const INT& n2)
{
  vector<Quad> ans;
  for(Quadlooper alpha(I2long(n1), I2long(n2), 1); alpha.ok(); ++alpha)
    ans.push_back(alpha);
  return ans;
}
