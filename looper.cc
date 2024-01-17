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
  if(Quad::t) // then d is odd so a=b(mod 2)
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

// Lists of quads of one norm or a range (up to units)

// return 1 iff a is *not* the chosen one of a conjugate pair
int other_conj(const Quad& a)
{
  if (Quad::nunits >2)
    return a.re()<a.im();
  else
    return a.im()<0;
}

vector<Quad> quads_of_norm(const INT& n, int conj)
{
  vector<Quad> ans;
  long d=Quad::d, t=Quad::t;

  INT n4 = (t? 4*n: n);
  INT bmax = isqrt(n4/d), a, b;
  for ( b=0; b<=bmax; b+=1)
    {
      if(isqrt(n4-d*b*b, a))
        {
          if (d==1 && a==0)
            continue;
          if(t) // then d is odd so a=b(mod 2)
            a = (a-b)/2;
          Quad z(a,b);
          // Keep this unless conj==0 and other_conj(z), i.e. this is not the chosen conjugate
          if (conj || !other_conj(z))
            ans.push_back(z);
        }
    }
  if (!conj)
    ans.erase(std::remove_if(ans.begin(), ans.end(), [](Quad val) { return other_conj(val); }),
              ans.end());
  return ans;
}

//#define DEBUG_QUADS_OF_NORM_BETWEEN

vector<Quad> quads_of_norm_between(const INT& n1, const INT& n2, int conj, int sorted)
{
#ifdef DEBUG_QUADS_OF_NORM_BETWEEN
  cout<<"Finding quads with norm between "<<n1<<" and "<<n2<<endl;
#endif
  vector<Quad> ans;
  long d=Quad::d, t=Quad::t;
  int extra_units = (d==1 || d==3);

  INT n1x = (t ? 4*n1 : n1), n2x = (t ? 4*n2 : n2);
  INT b, bmax = isqrt(n2x/d);
#ifdef DEBUG_QUADS_OF_NORM_BETWEEN
  cout<<"bmax = "<<bmax<<endl;
#endif
  for (b=0; b<=bmax; b+=1)
    {
#ifdef DEBUG_QUADS_OF_NORM_BETWEEN
      cout<<"b="<<b<<endl;
#endif
      INT db2 = d*b*b;
      INT aminsq = n1x-db2, amax = isqrt(n2x-db2);
#ifdef DEBUG_QUADS_OF_NORM_BETWEEN
      cout << "aminsq = "<<aminsq<<endl;
#endif
      INT amin = 0;
      if (aminsq.sign()>0)
        {
          amin = isqrt(aminsq);
          if (amin*amin!=aminsq) // isqrt rounds down but we want to round up
            amin+=1;
        }
#ifdef DEBUG_QUADS_OF_NORM_BETWEEN
      cout << " then amin = "<<amin<<", amax = "<<amax<<endl;
#endif
      if (extra_units && amin<b)
        {
          amin = b;
        }
      if (t)
        {
          if ((amin-b)%2) // ensure a=b (mod 2)
            amin +=1;
          amin = (amin-b)/2; // so this is rounded up
          amax = (amax-b)/2; // rounded down
        }

#ifdef DEBUG_QUADS_OF_NORM_BETWEEN
      cout<<"amin = "<<amin<<", amax = "<<amax<<endl;
#endif
      for (INT a = amin; a<=amax; a+=1)
        {
          if (extra_units && a==0)
            continue;
          Quad val = Quad(a, b);
#ifdef DEBUG_QUADS_OF_NORM_BETWEEN
          cout<<" val = "<<val<<" with norm "<<val.norm()<<endl;
#endif
          assert (val.norm()>=n1 && val.norm()<=n2 && pos(val));
          ans.push_back(val);
          if (b!=0 && conj)
            {
              val = makepos(val.conj());
#ifdef DEBUG_QUADS_OF_NORM_BETWEEN
              cout<<" conj val = "<<val<<" with norm "<<val.norm()<<endl;
#endif
              assert (val.norm()>=n1 && val.norm()<=n2 && pos(val));
              ans.push_back(val);
            }
        }
    }
#ifdef DEBUG_QUADS_OF_NORM_BETWEEN
  cout<<"initial list: "<<ans<<endl;
#endif

  if (sorted)
    std::sort(ans.begin(), ans.end(), Quad_cmp);

  // This works but is inefficient unless n1==n2
  // for(Quadlooper alpha(I2long(n1), I2long(n2), 1); alpha.ok(); ++alpha)
  //   ans.push_back(alpha);

#ifdef DEBUG_QUADS_OF_NORM_BETWEEN
  cout<<"final list:   "<<ans<<endl;
#endif

  return ans;
}
