// lf1.cc : classes lf1 and period_direct for computing L(F,1), L(F,chi,1) and periods of newforms

#include <values.h>
#include "eclib/interface.h"
#include "lf1.h"

// If Quad::t==0 then Quad (a,b) is x+y*sqrt(-d) so re,im parts are x, y*rootd
// If Quad::t==1 then Quad (a,b) is x+y*(1+sqrt(-d))/2 so re,im parts are x+y/2, y*rootd/2

// real part of Quad x
REAL real_part(const Quad& z)
{
  REAL x(real(z));
  if(Quad::t) // then x is coeff of 1/2
    x += REAL(imag(z))/2;
  return x;
}

// imag part of Quad x
REAL imag_part(const Quad& z, int norootd)
{
  REAL y(imag(z));
  if (norootd)  // return actual imaginary part divided by sqrt|D|
    return y/2;
  if(Quad::t) // then y is coeff of (1+sqrt(-d))/2
    y /= 2;
  return y*rootd();
}

// absolute value of Quad x
REAL realnorm(const Quad& z) {  return sqrt(z.norm());}

REAL lf1::K(REAL x)
{
  return (sfe==1? K1(x) : K0(x)/x);
}

int lf1::chi(const Quad& n)
{
  if (chi_is_trivial) return 1;
  return chitable[std::distance(lambdares.begin(), std::find(lambdares.begin(), lambdares.end(), n%lambda))];
}

void lf1::use(const Quad& n, int an)
{
  if (an ==0 || (n.norm() >= limitnorm))
    return;
  REAL rn = realnorm(n);
  REAL term = an * K(factor*rn) / rn;
  if (!chi_is_trivial && ar==0 && an)
    term *= chi(n);
  if (debug>1)
    {
      cout << "\nUsing term for n = " << n << ", |n| = "<<rn<<", a_n = " << an;
      if (!chi_is_trivial) cout << ", chi(n)="<<chi(n);
      cout << endl;
      REAL cn = REAL(an)/rn;
      cout << "Term    = " << term << " = "<<cn<<" * K("<<factor*rn<<") = "<<cn<<" * "<<K(factor*rn)<<endl;
      cout << "Sum was = " <<sum << endl;
    }
  sum += term;
  if (debug>1)
    cout << "Now sum = " <<sum << endl;
}

void lf1::add(const Quad& n, int pindex, int y, int z)
{
  int ip,istart=pindex;
  if ( y!=0 )
    {
      use(n,y);
      istart=0;
    }
  Quad p;
  long maxnorm=limitnorm/I2long(n.norm());
  for(p=quadprimes[ip=istart];
      (ip<=pindex) && (p.norm()<=maxnorm);
      p=quadprimes[++ip])
  {
   int x = y * aplist[ip];
   if ( (ip==pindex)  && ndiv(p,N.gen()))
     {
       x -=  I2long(p.norm())*z;
     }
   add(p*n,ip,x,y);
  }
}

//constructor -- does all the work after initializing:

lf1::lf1 (newform* f, int r, int db)
  :N(f->nf->N),
   debug(db),
   lambda(f->lambda),
   aplist(f->aplist),
   loverp(f->loverp),
   sfe(f->sfe),
   ar(r)
{
  Quad nu = N.gen();
  REAL modlambda = realnorm(lambda);
  factor = 4*REAL::Pi()/(rootdisc()*sqrt(realnorm(nu)));

  chi_is_trivial = lambda.is_unit();
  if (!chi_is_trivial && ar==0)
    {
      loverp = rational(f->lambdadot, f->cuspidalfactor);
      lambdares = residues(lambda);
      chitable = makechitable(lambda, lambdares);
      factor /= realnorm(lambda);
    }
  long maxnormp = I2long(quadprimes.back().norm());
  limitnorm = maxnormp;
  if(debug)
  {
    cout << "In lf1 constructor with ar = "<<ar<<endl;
    if (ar==0)
      cout << "chi_is_trivial = " << chi_is_trivial << endl;
    cout << "lambda = " << lambda << endl;
    cout << "loverp = " << loverp << endl;
    cout << "rootdisc = sqrt{|D_K|} = " << rootdisc() <<endl;
    if (chi_is_trivial || ar==1)
      {
        cout << "factor = (4*pi) / (rootdisc() * sqrt{|n|}) = " << factor <<endl;
      }
    else
      {
        cout << "factor = (4*pi) / (rootdisc() * sqrt{|n|} * |lambda|) = " << factor <<endl;
        cout << "Table of chi mod lambda (lambda = "<<lambda<<")\n";
        cout << "residue \t chi \n";
        for(unsigned int ires=0; ires<lambdares.size(); ires++)
          cout << lambdares[ires] << "\t" << chitable[ires] << endl;
      }
    cout<<"number of ap = "<<aplist.size()<<endl;
    cout<<"Summation using terms up to norm      "<<limitnorm<<endl;
    cout<<"          using a_p for norm(p) up to "<<maxnormp<<endl;
  }

//Initialize sum and use n=1 term:
  sum=0;
  use(Quad(1),1);

//add terms, one prime at a time:
  for(unsigned int i=0; i<aplist.size(); i++)
    {
      Quad p = quadprimes[i];
      int ap = aplist[i];
      if(debug>1)
        cout << "p= " << p << ",\ta_p = " << ap;
      add(p,i,ap,1);
      if(debug>1)
        cout << "\t\tSum = " << sum << endl;
      if (p.norm()>=maxnormp)
        break;
    }

//Scale to get values of L(f_chi,1) and the period, or L'(F,1):
  sum = 2*factor*sum;
  if (ar==0)
    {
      lf1chi = sum;
      if (loverp)
        {
          period=abs(lf1chi*modlambda);
          period/=num(loverp);
          period*=den(loverp);
        }
      else
        period = 0; // not determined
      if(debug)
        {
          cout<<"lf1chi = "<<lf1chi<<", L/P ratio = "<< loverp;
          if (loverp) cout<<", period = "<<period;
          cout<<endl;
        }
    }
  else
    {
      ldash1 = 2*sum;
      lf1chi = 0;
      period = 0; // not determined
    }
}

// psi(z) = exp(2*pi*i*Tr(z)) = cos(4*pi*x) + i*sin(4*pi*x)  where z=x+iy

// psi1(z) = psi(z)+psi(-z) = 2*cos(4*pi*x)

// psi_tilde(z) = sum of psi(eps*z) over units eps, so usually
// psi_tilde=psi1

REAL psi0(const REAL& x)
{
  return 2*cos(4*REAL::Pi()*x);
}

REAL psi1(const Quad& alpha, const Quad& a, const Quad& c)
{
  return psi0(imag_part(alpha*a*c.conj(), 1) / REAL(c.norm()));
}

REAL psi_tilde(const Quad& alpha, const Quad& a, const Quad& c)
{
  REAL ans = psi1(alpha, a, c);
  if (Quad::nunits == 2) // D not -3 or -4
    return ans;
  ans += psi1(fundunit*alpha, a, c);
  if (Quad::nunits == 4) // D=-4
    return ans;
  ans += psi1(fundunit*fundunit*alpha, a, c);
  return ans;
}

REAL period_direct::psi_factor(const Quad& alpha)
{
  // NB we divide by the number of units, correcting the formula (2.9) in Cremona-Whitley
  return (psi_tilde(alpha, Md, Mc) - psi_tilde(alpha, Ma, Mc)) / Quad::nunits;
}

period_direct::period_direct(newform*f, int db)
  :N(f->nf->N),
   debug(db),
   aplist(f->aplist),
   fa(f->a), fb(f->b),fc(f->c), fd(f->d),
   period_multiple(abs(f->matdot))
{
  maxnormp = I2long(quadprimes.back().norm());
  limitnorm = maxnormp;
  nu = N.gen();
}

// Compute the period along {.,g(.)} for g=[a,b;c,d] in Gamma_0(N)
REAL period_direct::compute_period(const Quad& a, const Quad& b, const Quad& c, const Quad& d)
{
  REAL t0 = realnorm(c).inv();
  Ma = a; Mb = b; Mc = c; Md = d;

  // The n'th term involves K_1(|n|*factor), and also the value of the
  // integral is factor*sum:

  factor = 4*REAL::Pi()*t0 / rootdisc();

  if (debug)
    {
      cout << "Matrix M = [a,b;c,d] = ["<<Ma<<", "<<Mb<<"; "<<Mc<<", "<<Md<<"]"<<endl;
      cout << "M in Gamma_0(N)? " << ((Ma*Md-Mb*Mc).is_one() && N.contains(Mc)) << endl;
      cout << "rootdisc = sqrt{|D_K|} = " << rootdisc() <<endl;
      cout << "factor = (4*pi*t0) / rootdisc = " << factor <<endl;
      cout << "number of ap =                          "<<aplist.size()<<endl;
      cout << "Integration using terms up to norm      "<<limitnorm<<endl;
      cout << "            using a_p for norm(p) up to "<<maxnormp<<endl;
    }

  //Initialize sum and use n=1 term:
  sum = REAL(0);
  use(Quad(1),1);

  //add terms, one prime at a time:
  for(unsigned int i=0; i<aplist.size(); i++)
    {
      Quad p = quadprimes[i];
      add(p,i,aplist[i],1);
      if (debug>1)
        cout << "\tSum = " << sum << endl;
      if (p.norm()>=maxnormp)
        break;
    }
  REAL period = abs(factor * sum);
  if (debug)
    {
      cout << "Sum      = " << abs(sum) << endl;
      cout << "Period = " << period << endl;
      cout << "Period*  = " << 2*period/rootdisc() << endl;
    }
  return period;
}

REAL period_direct::compute_base_period()
{
  REAL period = compute_period(fa,fb,fc,fd);
  REAL base_period = period/period_multiple;
  if (debug)
    {
      cout << "period along ["<<fa<<","<<fb<<";"<<fc<<","<<fd<<"] = " << period << endl;
      cout << "base period multiple = " << period_multiple << endl;
      cout << "base period          = " << base_period << endl;
    }
  return base_period;
}

void period_direct::use(const Quad& n, int an)
{
  if (an ==0 || (n.norm() >= limitnorm))
    return;
  REAL rn = realnorm(n);
  REAL cn = REAL(an)/rn;
  REAL term = cn * K1(factor*rn) * psi_factor(n);
  if (debug>1)
    {
      cout << "\nUsing term for n = " << n << ", |n| = "<<rn<<", a_n = " << an << endl;
      cout << "Term    = " << term << " = "<<cn<<" * K1("<<factor*rn<<") * psi_factor(n) = "<<cn<<" * "<<K1(factor*rn)<<" * "<< psi_factor(n) << endl;
      cout << "Sum was = " <<sum << endl;
    }
  sum += term;
  if (debug>1)
    cout << "Now sum = " <<sum << endl;
}

// NB the code here is identical to lf1::add()

void period_direct::add(const Quad& n, int pindex, int y, int z)
{
  if ( y!=0 )
    use(n,y);

  long maxnorm=limitnorm/I2long(n.norm());
  for(int ip = (y==0? pindex : 0); ip<=pindex; ip++)
  {
    Quad p = quadprimes[ip];
    long normp = I2long(p.norm());
    if (normp > maxnorm)
      break;
    int x = y * aplist[ip];
    if ( (ip==pindex)  && ndiv(p,N.gen()))
      {
        x -=  normp*z;
      }
    add(p*n,ip,x,y);
  }
}


