#include "looper.h"

// Set bmin and bmax for this value of n, and initialise the b-loop
void Quadlooper::setblims()
{
  n4 = (d%4==3? 4*n: n);
  INT zero(0);
  switch (d)
    {
    case 1:
      b = bmin = db2 = zero;
      bmax = include_conjugates? isqrt(n-1) : isqrt(n/2);
      break;
    case 3:
      b = bmin = db2 = zero;
      bmax = include_conjugates? isqrt(n-1) : isqrt(n/3);
      break;
    default:
      bmax = isqrt(n4/d);
      b = bmin = include_conjugates? -bmax : zero;
      db2 = d*bmin*bmin;
      if (db2==n4)
        {
          bmin+=1;
          db2 = d*bmin*bmin;
        }
    }
  // cout<<"n="<<n<<": Setting bmin="<<bmin<<", bmax="<<bmax<<endl;
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
  val = makepos(Quad(a,b));
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


vector<Quad> Quadlooper::values_with_current_norm(int sorted)
{
  // must make a copy
  INT m = n;
  auto values = values_with_norm_up_to(m);
  if (sorted)
    std::sort(values.begin(), values.end(), Quad_cmp);
  return values;
}

vector<Quad> Quadlooper::values_with_norm_up_to(const INT& m, int sorted)
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
  if (sorted)
    std::sort(values.begin(), values.end(), Quad_cmp);
  return values;
}

// Lists of quads of one norm or a range (up to units)

// we solve a^2+db^2=n4 where n4=n (t=0),  or =4n (t=1), t=Quad::t
// with 0<=b<=sqrt(n4/d) and a>=0 with additional restrictions
// d=1: a>0 if conj, a>=b if not conj;
// d=3: a>b if conj, a>=3b if not conj.
// Then use a+b*w (t=0) or (a-b)/2 + b*w (t=1).
// if d not 1,3, a>0, and conj, use also (-a,b) --> -a+b*w or -(a+b)/2+b*w

vector<Quad> quads_of_norm(const INT& n, int conj, int sorted)
{
  vector<Quad> ans;
  long d=Quad::d, t=Quad::t;

  INT n4 = (t? 4*n: n);
  INT bmax = isqrt(n4/d), a, b, a1;
  for ( b=0; b<=bmax; b+=1)
    {
      if(isqrt(n4-d*b*b, a))
        {
          // now we have a^2+db^2=n4 with a,b>=0
          if ((d==1) && (a < (conj? ONE : b)))
            continue;
          if ((d==3) && (a < (conj? b+1 : 3*b)))
            continue;
          a1 = (t? (a-b)/2 : a);
          ans.push_back(Quad(a1,b));
          if ((d!=1)&&(d!=3)&&(a>0)&&(b>0)&&conj)
            {
              a1 = (t? (-a-b)/2 : -a);
              ans.push_back(Quad(a1,b));
            }
        }
    }
  if (sorted)
    std::sort(ans.begin(), ans.end(), Quad_cmp);

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

  INT n1x = (t ? 4*n1 : n1), n2x = (t ? 4*n2 : n2);
  INT b, bmax = isqrt(n2x/d), db2, aminsq, amax, amin;
#ifdef DEBUG_QUADS_OF_NORM_BETWEEN
  cout<<"bmax = "<<bmax<<endl;
#endif
  for (b=0; b<=bmax; b+=1)
    {
#ifdef DEBUG_QUADS_OF_NORM_BETWEEN
      cout<<"b="<<b<<endl;
#endif
      db2 = d*b*b;
      aminsq = n1x-db2;
      amax = isqrt(n2x-db2);
#ifdef DEBUG_QUADS_OF_NORM_BETWEEN
      cout << "aminsq = "<<aminsq<<endl;
#endif
      amin = 0;
      if (aminsq.sign()>0)
        {
          amin = isqrt(aminsq);
          if (amin*amin!=aminsq) // isqrt rounds down but we want to round up
            amin+=1;
        }
      // special cases d=1,3
      // d=1: a>0 if conj, a>=b if not conj;
      // d=3: a>b if conj, a>=3b if not conj.
      if (d==1)
        amin = max(amin, (conj? ONE : b));
      if (d==3)
        amin = max(amin, (conj? b+1 : 3*b));
#ifdef DEBUG_QUADS_OF_NORM_BETWEEN
      cout << " then amin = "<<amin<<", amax = "<<amax<<endl;
#endif
      if (t)
        {
          if ((amin-b)%2)  // ensure a=b (mod 2)
            amin +=1;
          if ((amax-b)%2)  // ensure a=b (mod 2)
            amax -=1;
          amin = (amin-b)/2; // so this is rounded up
          amax = (amax-b)/2; // and this rounded down
        }

#ifdef DEBUG_QUADS_OF_NORM_BETWEEN
      cout<<"amin = "<<amin<<", amax = "<<amax<<endl;
#endif
      for (INT a = amin; a<=amax; a+=1)
        {
          Quad value(a, b);
#ifdef DEBUG_QUADS_OF_NORM_BETWEEN
          cout<<" value = "<<value<<" with norm "<<value.norm()<<endl;
#endif
          if (!(value.norm()>=n1 && value.norm()<=n2 && pos(value)))
            {
              cout<<"In quads_of_norm_between() with n1="<<n1<<", n2="<<n2<<", conj="<<conj<<", sorted="<<sorted<<endl;
              cout<<"bmax="<<bmax<<", amin="<<amin<<", amax="<<amax<<endl;
              cout<<"a="<<a<<", b="<<b<<endl;
              cout << "value = "<<value<<" has norm "<<value.norm()<<endl;
              assert (value.norm()>=n1 && value.norm()<=n2 && pos(value));
            }
          ans.push_back(value);
          INT a1 = (t? 2*a+b : a);
          if ((d!=1)&&(d!=3)&&(a1>0)&&(b>0)&&conj)
            {
              value = -value.conj();
#ifdef DEBUG_QUADS_OF_NORM_BETWEEN
              cout<<" conj value = "<<value<<" with norm "<<value.norm()<<endl;
#endif
              assert (value.norm()>=n1 && value.norm()<=n2 && pos(value));
              ans.push_back(value);
            }
        }
    }
#ifdef DEBUG_QUADS_OF_NORM_BETWEEN
  cout<<"unsorted list: "<<ans<<endl;
#endif

  if (sorted)
    {
      std::sort(ans.begin(), ans.end(), Quad_cmp);
#ifdef DEBUG_QUADS_OF_NORM_BETWEEN
      cout<<"sorted list:   "<<ans<<endl;
#endif
    }
  return ans;
}
