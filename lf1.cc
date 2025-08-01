// lf1.cc : class period_via_lfchi for integrating newforms

#include <values.h>
#include "eclib/interface.h"
#include "eclib/kbessel.h"
#include "lf1.h"

// absolute value of Quad x
double realnorm(const Quad& z) {  return sqrt(to_double(z.norm()));}

// convert to a complex number
bigcomplex to_bigcomplex(const Quad& a)
{
  static const bigfloat rootd = sqrt(to_bigfloat(Quad::d));
  bigfloat x = to_bigfloat(I2long(real(a))), y = to_bigfloat(I2long(imag(a)));
  if(Quad::t) // then y is coeff of (1+sqrt(-d))/2
    {
      y /= 2;
      x += y;
    }
  bigcomplex z(x, y*rootd);
  //cout << "a = " << a << " --> " << z << endl;
  return z;
}

int period_via_lf1chi::chi(const Quad& n)
{
  if (chi_is_trivial) return 1;
  return chitable[std::distance(lambdares.begin(), std::find(lambdares.begin(), lambdares.end(), n%lambda))];
}

void period_via_lf1chi::use(const Quad& n, int an)
{
  if (n.norm() >= limitnorm)
    return;
  double rn = realnorm(n);
  double cn = double(an)/rn;
  if (!chi_is_trivial) cn *= chi(n);
  if (debug>1)
    {
      cout << "\nUsing term n = " << n << ", a_n = " << an;
      if (!chi_is_trivial) cout << ", chi(n)="<<chi(n);
      cout << endl;
    }
  sum += cn * K1(factor*rn);
}

void period_via_lf1chi::add(const Quad& n, int pindex, int y, int z)
{
  int ip,istart=pindex;
  Quad p;
  long maxnorm=limitnorm/I2long(n.norm());
  if ( y!=0 )
    {
      use(n,y);
      istart=0;
    }
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

period_via_lf1chi::period_via_lf1chi (newform* f, int db)
  :N(f->nf->N),
   debug(db),
   lambda(f->lambda),
   aplist(f->aplist),
   loverp(f->loverp)
{
  Quad nu = N.gen();
  chi_is_trivial = lambda.is_unit();
  double rootdisc = sqrt((double)(I2long(Quad::absdisc)));
  double modlambda = realnorm(lambda);
  factor = 4*PI/(rootdisc*sqrt(realnorm(nu))*modlambda);
  lambdares = residues(lambda);
  chitable = makechitable(lambda, lambdares);
  long maxnormp = I2long(quadprimes.back().norm());
  limitnorm = maxnormp;
  if(debug)
  {
    cout << "lambda = " << lambda << endl;
    cout << "rootdisc = sqrt{|D_K|} = " << rootdisc <<endl;
    cout << "factor = (4*pi) / (rootdisc * sqrt{|n|} * |lambda|) = " <<factor <<endl;
    if (!chi_is_trivial)
      {
        cout << "Table of chi mod lambda (lambda = "<<lambda<<")\n";
        cout << "residue \t chi \n";
        for(unsigned int ires=0; ires<lambdares.size(); ires++)
          cout << lambdares[ires] << "\t" << chitable[ires] << endl;
      }
    cout<<"number of ap = "<<aplist.size()<<endl;
    cout<<"Integration using terms up to norm      "<<limitnorm<<endl;
    cout<<"            using a_p for norm(p) up to "<<maxnormp<<endl;
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

//Scale to get values of L(f_chi,1) and the period:
  lf1chivalue = 2*factor*sum;
  if (loverp)
    period=abs(lf1chivalue*modlambda/double(loverp));
  else
    period = 0;
  if(debug)
  {
    cout<<"lf1chivalue = "<<lf1chivalue<<", L/P ratio = "<< loverp;
    if (loverp) cout<<", period = "<<period;
    cout<<endl;
  }
}

// psi(z) = exp(2*pi*i*Tr(z)) = cos(4*pi*x) + i*sin(4*pi*x)  where z=x+iy

// psi1(z) = psi(z)+psi(-z) = 2*cos(4*pi*x)

// psi_tilde(z) = sum of psi(eps*z) over units eps, so usually
// psi_tilde=psi1

double psi1(bigcomplex z)
{
  return to_double(cos(4*Pi()*real(z)));
}

double psi_tilde(bigcomplex z)
{
  static const bigcomplex eps = to_bigcomplex(fundunit);
  switch (Quad::nunits) {
  case 4:
    return psi1(z) + psi1(eps*z);
  case 6:
    return psi1(z) + psi1(eps*z) + psi1(eps*eps*z);
  case 2:
  default:
    return psi1(z);
  }
}

double period_direct::psi_factor(const Quad& n)
{
  bigcomplex cn = to_bigcomplex(n);
  return psi_tilde(cn*z1)-psi_tilde(cn*z2);
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
  rootdisc = sqrt((double)(I2long(Quad::absdisc))); // = sqrt(|D|), real
  eta = bigcomplex(to_bigfloat(0),to_bigfloat(rootdisc));    // = sqrt(D), pure imagainry
}

// Compute the period along {.,g(.)} for g=[a,b;c,d] in Gamma_0(N)
double period_direct::compute_period(const Quad& a, const Quad& b, const Quad& c, const Quad& d)
{
  bigcomplex ca = to_bigcomplex(a);
  bigcomplex cc = to_bigcomplex(c);
  bigcomplex cd = to_bigcomplex(d);

  double t0 = 1/realnorm(c);
  z1 = - cd / (eta*cc);
  z2 = + ca / (eta*cc);
  // The path of integration is from (z1,t0) to (z2,t0) = M(z1,t0); we integrate
  // from (z1,t0) to oo and from (z2,t0) to oo and subtract

  // The n'th term involves K_1(|n|*factor), and also the value of the
  // integral is factor*sum:

  factor = 4*PI*t0 / rootdisc;

  if (debug)
    {
      cout << "Matrix M = [a,b;c,d] with a="<<a<<", b="<<b<<", c="<<c<<", d="<<d<<endl;
      cout << "M in Gamma_0(N)? " << ((a*d-b*c).is_one() && N.contains(c)) << endl;
      cout << "rootdisc = sqrt{|D_K|} = " << rootdisc <<endl;
      cout << "(z1,t0) = ("<<z1<<","<<t0<<")"<<endl;
      cout << "(z2,t0) = ("<<z2<<","<<t0<<")"<<endl;
      cout << "factor = (4*pi*t0) / rootdisc = " << factor <<endl;
      cout << "number of ap =                          "<<aplist.size()<<endl;
      cout << "Integration using terms up to norm      "<<limitnorm<<endl;
      cout << "            using a_p for norm(p) up to "<<maxnormp<<endl;
    }

  //Initialize sum and use n=1 term:
  sum = 0;
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
  double period = abs(factor * sum);
  if (debug)
    {
      cout << "Sum      = " << abs(sum) << endl;
      cout << "Period = " << period << endl;
      cout << "Period*  = " << 2*period/rootdisc << endl;
    }
  return period;
}

double period_direct::compute_base_period()
{
  double period = compute_period(fa,fb,fc,fd);
  double base_period = period/period_multiple;
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
  double rn = realnorm(n);
  sum += (double(an) * K1(rn*factor) * psi_factor(n) / rn);
}

// NB the code here is identical to period_via_lf1chi::add()

void period_direct::add(const Quad& n, int pindex, int y, int z)
{
  if ( y!=0 )
    use(n,y);

  long normn = I2long(n.norm());
  for(int ip = (y==0? pindex : 0); ip<=pindex; ip++)
  {
    Quad p = quadprimes[ip];
    long normp = I2long(p.norm());
    if (normp*normn > limitnorm)
      break;
    int x = y * aplist[ip];
    if ( (ip==pindex)  && ndiv(p,N.gen()))
      {
        x -=  normp*z;
      }
    add(p*n,ip,x,y);
  }
}


