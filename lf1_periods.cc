#include <fstream>
#include "newforms.h"
#include "lf1.h"
//#define LOOPER
#ifdef LOOPER
#include "looper.h"
#endif

#define RECOMPUTE_RATIOS

scalar modulus = default_modulus<scalar>();

int main ()
{
 cout.precision(10);
 int d,maxpnorm=10000;
 cout << "Enter field: " << flush;
 cin >> d;
 Quad::field(d,maxpnorm);
 Quad n; int verbose=0;
 cout << "Verbose? ";
 cin>>verbose;
 long denom_norm_bound = 0;
 cout << "Bound on denominator norm for extra periods (0 for none)? ";
 cin >> denom_norm_bound;
#ifdef LOOPER
 long firstn, lastn;
 cout<<"Enter first and last norm for Quads: ";
 cin >> firstn >> lastn;
 cerr<<endl;
 int both_conj=0;
 Qidealooper loop(firstn, lastn, both_conj, 1); // sorted within norm
 while( loop.not_finished() )
   {
     Qideal N = loop.next();
#else
     Qideal N;
     while(cerr<<"Enter level (ideal label or generator): ", cin>>N, !N.is_zero())
#endif
   {
     cout << ">>>> Level " << label(N) <<" = "<<gens_string(N)<<", norm = "<<N.norm()<<" <<<<" << endl;
     newforms nf(N, modulus, 1);
     nf.read_from_file();
#ifdef RECOMPUTE_RATIOS
     nf.makebases();
#endif
     if (verbose)
       nf.display();
     for(int i=0; i<nf.n1ds; i++)
       {
         cout << "\nForm number " << i+1 << ": " << endl << endl;
         Quad lambda = nf.nflist[i].lambda;
         REAL abs_lambda = realnorm(lambda);
         int trivial_twist = (lambda==Quad(1));
         string chi_string = (trivial_twist? "": "chi,");
         int sign = nf.nflist[i].sfe;
         int deriv = (sign>0? 0 : 1);
         if (deriv)
           {
             cout << "Computing L'(F,1)" << endl;
           }
         else
           {
             cout << "Computing L(F,"<<chi_string<<"1) and period" << endl;
             if(!trivial_twist)
               {
                 cout<<"Using twisting prime lambda = "<<lambda<<endl;
               }
           }
         lf1 summer(&(nf.nflist[i]), deriv, verbose);

         if (deriv)
           {
             cout << "L'(F,1) = " << summer.get_ldash1() << endl;
             continue;
           }

         REAL lf1chi = summer.get_lf1chi();
         cout << "L(F," << chi_string << "1) = " << lf1chi << endl;
         REAL lf1chi_abs_lambda = lf1chi*abs_lambda;;

         REAL P_from_L = summer.get_period();
         rational ratio = nf.nflist[i].loverp;
         if (num(ratio)==0)
           ratio = rational(I2long(nf.nflist[i].lambdadot), I2long(nf.nflist[i].cuspidalfactor));
         cout << "Period = " << P_from_L
              << " (via L(F,chi,1), using L/P ratio = " << ratio << ")"<< endl << endl;

         cout<<"Finding period by direct integration, using stored matrix and scaling factor:"<<endl;
         period_direct summer2(&(nf.nflist[i]), verbose);
         REAL P0 = summer2.compute_base_period();
         Quad b0=nf.nflist[i].b, d0=nf.nflist[i].d;
         scalar matdot0 = nf.nflist[i].matdot;
         cout << "Base period P0 = " << P0 << " = I_F({0,"<<RatQuad(b0,d0)<<"}) / "<<matdot0<<endl;
         cout << "(computed) L/P ratio = " << lf1chi_abs_lambda/P0 << endl << endl;

         if (!denom_norm_bound)
           continue;

         cout << "Computing extra periods I_F({0,b/d}) for N(d) <= "<<denom_norm_bound<<endl;
         Quad a, b, c, d;
         long gcd_multiple = 0;
         for (Quadlooper dl(2, denom_norm_bound, 1); dl.ok(); ++dl)
           { d=(Quad)dl;
             Qideal D(d);
             if (N.is_coprime_to(D))
               {
                 vector<Quad> reslist = residues(d);
                 for( const auto& bi : reslist)
                   {
                     b = bi;
                     Qideal bN = b*N;
                     if (D.is_coprime_to(bN, a, c))
                       // found a candidate q=b/d: a+c=1 with d|a and b|c and c/b in N
                       {
                         c /= -b;
                         a /= d; // now a*d-b*c=1 with c in N
                         assert (a*d-b*c==Quad::one);
                         RatQuad q(b,d);
#ifdef RECOMPUTE_RATIOS
                         long matdot = I2long(abs((nf.h1->chain(q, 1))[i+1]) / nf.nflist[i].cuspidalfactor);
                         gcd_multiple = gcd(gcd_multiple, matdot);
#endif
                         REAL period__b_d = summer2.compute_period(a,b,c,d);
                         cout << " period " <<period__b_d
                              << " = I_F({0,"<< q <<"}) = P0 * " << period__b_d/P0
#ifdef RECOMPUTE_RATIOS
                              << " (expected multiple = " << matdot << ", gcd so far = "<<gcd_multiple<<")"
#endif
                              << endl;
                       } // b coprime to d test
                   } // loop over b
               } // d coprime to N test
           } // loop over d
         if (gcd_multiple)
           {
             if (gcd_multiple==1)
               cout << "These periods are all multiples of P0 by coprime integers, "
                    << "suggesting that P0 is correct and L/P0 = " << lf1chi_abs_lambda/P0 << endl;
             else
               cout << "These periods are all multiples of " <<gcd_multiple
                    <<"*P0 by coprime integers, suggesting that the correct base period is "<<gcd_multiple
                    <<"*P0 = " << gcd_multiple*P0
                    << ", and the correct L/P is " << lf1chi_abs_lambda/(gcd_multiple*P0) << endl;
           }
       }
   }
 cout << endl;
}
