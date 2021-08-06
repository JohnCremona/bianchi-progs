// FILE MAT22.CC

#include "mat22.h"
#include "qideal.h"

matop::matop(const Quad& p, const Quad& n)
{
 if (p==n)
   {
     mats.resize(1, mat22(0,-1,n,0));
   }
 else
 if (div(p,n))   // W involution, 1 term
   {
      Quad u,v,a,b;
      for (u=1, v=n; div(p,v); v/=p, u*=p) ;
      quadbezout(u,v,a,b);
      mats.resize(1, mat22(u*a,-b,n,u));
   }
else                 // Hecke operator, p+1 terms
  {
    vector<Quad> resmodp = residues(p);
    vector<Quad>::const_iterator r=resmodp.begin();
    while(r!=resmodp.end())
      mats.push_back(mat22(1,*r++,0,p));
    mats.push_back(mat22(p,0,0,1));
  }
}

// partial implementation of Hecke and Atkin-Lehner matrices, assuming
// all ideals principal:

matop::matop(Qideal& P, Qideal& N)
{
  // cout<<"In matop constructor with P="<<P<<", N="<<N<<" = ("<<N.gen()<<")..."<<flush;
  if (P==N)
   {
     assert (N.is_principal());
     mat22 W(0,-1,N.gen(),0);
     mats.push_back(W);
   }
 else
   if (P.divides(N))   // W involution, 1 term
     {
       Qideal Q(1), M(N);
       while (P.divides(M))
         {
           M /= P;
           Q *= P;
         }
       assert (Q.is_principal()); // has side-effect of setting their gens!
       assert (M.is_principal());
       Quad u = Q.gen(), v = M.gen(), a,b;
       quadbezout(u,v,a,b);
       mat22 W(u*a,-b,u*v,u);
       assert (W.det()==u);
       mats.push_back(W);
     }
   else                 // Hecke operator, p+1 terms
     {
       assert (P.is_principal());
       Quad p = P.gen();
       vector<Quad> resmodp = P.residues();
       vector<Quad>::const_iterator r=resmodp.begin();
       while(r!=resmodp.end())
         mats.push_back(mat22(1,*r++,0,p));
       mats.push_back(mat22(p,0,0,1));
     }
  // cout<<mats<<endl;
}


RatQuad mat22::operator()(const RatQuad& q)const
{
  Quad r = q.num(), s = q.den();
  apply_left(r, s);
  return RatQuad(r,s, 1);
}
