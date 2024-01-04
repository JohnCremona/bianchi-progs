// FILE MAT22.CC

#include "primes.h"
#include "mat22.h"

RatQuad mat22::operator()(const RatQuad& q)const
{
  Quad r = q.num(), s = q.den();
  apply_left(r, s);
  return RatQuad(r,s);
}

RatQuad mat22::image_oo() const {return RatQuad(a,c);}
RatQuad mat22::preimage_oo() const {return RatQuad(-d,c);}
RatQuad mat22::image_0() const {return RatQuad(b,d);}
RatQuad mat22::preimage_0() const {return RatQuad(b,-a);}

//#define DEBUG_LIFT

// return a matrix [a, b; c, d] with det=1 and (c:d)=(cc:dd) in P^1(N)
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
      a = Quad::one;
      b = N.reduce(-x*z);
      c = N.reduce(c*y); // so b*c == -x*c = d*y-1
      d = a + b*c; // = d*y mod N
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

// return a matrix [a, b; c, d] with det=1 and c in M and (c:d)=(cc:dd) in P^1(N)
// If (u,v)!=(0,0) they should satisfy u+v=1 with u in N, v in M,
// otherwise such u,v will be computed and returned.
mat22 lift_to_Gamma_0(Qideal& M, Qideal& N, const Quad& cc, const Quad& dd, const Quad& u, const Quad& v)
{
  // CRT: lift (cc:dd) in P^1(N) to (c:d) in P^1(MN) which also lifts (0:1) in P^1(M), then lift that
  Quad uu(u), vv(v);
  if (uu.is_zero() && vv.is_zero())
    M.is_coprime_to(N, uu, vv); // u+v=1, u in M, v in N
  assert (M.contains(uu));
  assert (N.contains(vv));
  assert (uu+vv==Quad::one);
  Qideal MN = M*N;
  return lift_to_SL2(MN, cc*uu, dd*uu+vv);
}
