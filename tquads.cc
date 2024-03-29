#include "looper.h"
#include "geometry.h"
#include "primes.h"

int main ()
{
  long d, max;
 cout << "Enter field: " << flush;  cin >> d;
 cout << "Enter max. norm for primes: " << flush;  cin >> max;
 Quad::field(d,max);
 cout << "The field is "; Quad::displayfield(cout); cout << endl;
 Quad w = Quad::w;
 cout << "w   = " << w << endl;
 cout << "w*w = " << (w*w) << endl;
 Quad a,b,c,p;

 cout << "Enter a Quad, a: ";
 cin >> a;
 cout << "a = " << a << endl;
 cout << "real(a) = " << real(a) << endl;
 cout << "imag(a) = " << imag(a) << endl;
 cout << "quadconj(a) = " << quadconj(a) << endl;
 cout << "quadnorm(a) = " << quadnorm(a) << endl;
 cout << "HNF(a) = " << HNF(a) << endl;
 cout << "ideal_label(a) = " << ideal_label(a) << endl;

  if (Quad::class_number==1)
   {
     cout << "Here are the first 20 primes out of "<<quadprimes.size()
          << " (one per line):\n";
     auto pr = quadprimes.begin();
     while((pr-quadprimes.begin() < 21) && (pr!=quadprimes.end()))
       {
         p = *pr++;
         cout << "("<<p << ") has norm "<<p.norm()<<" and label "<<ideal_label(p)<<endl;
       }
   }
 else
   {
     cout << "---------------------------------------------------------------" << endl;
     Quadprimes::display();
     cout << "Here are the first 20 prime ideals:"<<endl;
     for (auto Pi = Quadprimes::list.begin();
          (Pi-Quadprimes::list.begin())<20 && (Pi != Quadprimes::list.end()); ++Pi)
       {
         Quadprime P = *Pi;
         vector<Quad> gg = P.gens();
         cout << P << " = " << ideal_label(P) << " = " << (Qideal)P << " = (" << gg[0] <<","<<gg[1] << ")";
         if (P.is_principal())
           cout << " = ("<<gg[0]<<") (principal)";
         else
           cout << " (not principal)";
         cout<<endl;
       }
     cout << "---------------------------------------------------------------" << endl;
   }

 b=a*a;
 cout << "b=a*a=" << b << endl;
 cout << "val(a,b) = " << val(a,b) << endl;
 cout << "div(a,b) = " << div(a,b) << endl;
 cout << "div(b,a) = " << div(b,a) << endl;
 c=b/a;
 cout << "b/a (=a?) = " << c;
 if (c==a) {cout<<" \t--OK!";}
 else      {cout<<" \t--NO!";}
 cout<<endl;
 c=a/THREE;
 cout << "a/3 = " << c << endl;
 c=a+b;
 cout << "c = a+b = " << c << endl;
 c=c-b;
 cout << "c-b (=a?)  = " << c;
 if (c==a) {cout<<" \t--OK!";}
 else      {cout<<" \t--NO!";}
 cout<<endl;
 cout << "pos(a) = " << pos(a) << endl;

//Test of residues
 vector<Quad> resmoda = residues(a);  INT norma=a.norm();
 cout << norma << " residues modulo a: "<<resmoda<<endl;

//Test of primes and divisor functions
 cout<<"Testing prime divisor functions.\nEnter a: "; cin>>a;

 if (Quad::class_number==1)
   {
     Quad pda = primdiv(a);
     cout << "A prime divisor of a: " << pda << endl;
     vector<Quad> plist = pdivs(a);
     cout << "The list of all prime divisors of a: " << plist << endl;
     vector<int> elist;
     for (const auto& p : plist)
       elist.push_back(val(p,a));
     cout << "Exponents: " << elist << endl;
     vector<Quad> dlist = posdivs(a);
     cout << "The list of "<<dlist.size()<<" divisors of a (up to units): \n";
     cout << dlist << endl;
     dlist = alldivs(a);
     cout << "The list of all "<<dlist.size()<<" divisors of a: \n";
     cout << dlist << endl;
     dlist = sqdivs(a);
     cout << "The list of "<<dlist.size()<<" divisors of a whose square divides: \n";
     cout << dlist << endl;
     dlist = sqfreedivs(a);
     cout << "The list of "<<dlist.size()<<" square-free divisors of a: \n";
     cout << dlist << endl;
   }
 else
   {
     Qideal A(a);
     Factorization F = A.factorization();
     cout << "Prime ideal factorization of A = (a):\t" << F << endl;
     cout<<endl;
     cout << "Ideals with the same norm as A:    \t"
          << Qideal_lists::ideals_with_norm(a.norm()) <<endl;
     cout << "Prime ideal divisors of A:         \t"
          << pdivs(A) <<endl;
     cout << "All ideal divisors of A:           \t"
          << alldivs(A) <<endl;
     cout << "All ideals whose square divides A: \t"
          << sqdivs(A) <<endl;
     cout << "All squarefree ideal divisors of A:\t"
          << sqfreedivs(A) <<endl;
   }

 if (!check_field(d))
   {
     cout << "Skipping gcd and bezout tests as this field not yet fully implemented"<<endl;
     exit(0);
   }

 //Test of gcd and bezout
 Quad g,h,x,y;
 cout << "Testing gcd and bezout."<<endl;
 while(cout<<"Enter Quads a and b (a=0 to stop): ", cin >> a >> b, !a.is_zero())
   {
     g = quadgcd(a,b);
     if (!g.is_zero())
       cout << "gcd("<<a<<","<<b<<") = "<< g <<endl;
     if (Quad::class_number==1)
       {
         g=quadbezout(a,b,x,y);
         cout << "bezout returns x="<<x<<", y="<<y<<", g="<<g<<".  ";
         if(g==a*x+b*y)cout<<"OK!"; else cout<<"WRONG!";
         cout<<endl;
       }
   }
 cout << "Systematic testing of bezout..."<<flush;
 if (Quad::class_number==1)
   {
     for (auto ap=quadprimes.begin(); ap-quadprimes.begin()<10 && ap!=quadprimes.end(); ++ap)
       for (auto bp=quadprimes.begin(); bp-quadprimes.begin()<10 && bp!=quadprimes.end(); ++bp)
         {
           a = *ap;
           b = *bp;
           cout<<"quadbezout("<<a<<","<<b<<") = " << quadbezout(a,b,x,y) << endl;
         }
   }
 else
   {
     Quad::setup_geometry();
     long minn(1), maxn(10);
     int s;
     mat22 M;
     for(Quadlooper alpha(minn,maxn,1); alpha.ok(); ++alpha)
       for(Quadlooper beta(minn,maxn,1); beta.ok(); ++beta)
         {
           a = alpha;
           b = beta;
           M = generalised_extended_euclid(a, b, s);
           cout<<"generalised_extended_euclid("<<a<<","<<b<<") = " << M << "; ";
           g = a; h = b;
           M.apply_left(g,h);
           while (!pos(h))
             {
               g*=fundunit;
               h*=fundunit;
             }
           if (h.is_zero())
             while (!pos(g))
               {
                 g*=fundunit;
               }
           RatQuad z(g,h);
           cout<<"("<<a<<","<<b<<") is ";
           if (s!=0)
             cout << "not principal, maps to singular point " << sigmas[s] << endl;
           else
             cout << "principal, maps to ("<<g<<")/("<<h<<") = " <<z<<endl;
           assert (sigmas[s]==z);
         }
   }
 cout<<"done.\n";
}
