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
  return bigcomplex(x, y*rootd);
}

void period_via_lf1chi::use(const Quad& n, int an)
{
  double rn = realnorm(n);
  if(0) // (debug)
  { cout << "\nUsing term n = " << n << ", a_n = " << an;
    cout << ", chi(n)="<<chi(n)<<endl;
  }
  if (n.norm()<limitnorm) sum += chi(n) * double(an) * K1(factor*rn) /rn;
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
   nap(f->aplist.size()),
   aplist(f->aplist),
   ratio(f->loverp)
{
  double rootdisc = sqrt((double)(I2long(Quad::absdisc))), modlambda = realnorm(lambda);
  factor = 4*PI/(rootdisc*sqrt(realnorm(N.gen()))*modlambda);
  lambdares = residues(lambda);
  chitable = makechitable(lambda, lambdares);
  long maxnormp = I2long(quadprimes[nap-1].norm());
  limitnorm = maxnormp;
  if(debug)
  {
    cout << "rootdisc = sqrt{|D_K|} = " << rootdisc <<endl;
    cout << "factor = (4*pi) / (rootdisc * sqrt{|n|} * |lambda|) = " <<factor <<endl;
    cout << "Table of chi mod lambda (lambda = "<<lambda<<")\n";
    cout << "residue \t chi \n";
    for(unsigned int ires=0; ires<lambdares.size(); ires++)
      cout << lambdares[ires] << "\t" << chitable[ires] << endl;
    cout<<"nap = "<<nap<<endl;
   cout<<"Integration using terms up to norm      "<<limitnorm<<endl;
   cout<<"            using a_p for norm(p) up to "<<maxnormp<<endl;
  }

//Initialize sum and use n=1 term:
  sum=0;
  use(Quad(1),1);

//add terms, one prime at a time:
  for(int i=0; i<nap; i++)
    {
      Quad p = quadprimes[i];
      if(debug>1)
        cout << "p= " << p << ",\ta_p = " << aplist[i];
      add(p,i,aplist[i],1);
      if(debug>1)
        cout << "\t\tSum = " << sum << endl;
      if (p.norm()>=maxnormp)
        break;
    }

//Scale to get values of L(f_chi,1) and the period:
  lf1chivalue = 2*factor*sum;
  period=abs(lf1chivalue*modlambda/double(ratio));
  if(debug)
  {
    cout<<"lf1chivalue = "<<lf1chivalue<<", ratio = "<<ratio;
    cout<<", period = "<<period<<endl;
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
  return psi_tilde(cn*theta1)-psi_tilde(cn*theta2);
}

period_direct::period_direct(newform*f, int db)
  :N(f->nf->N),
   debug(db),
   aplist(f->aplist)
{
  nap = aplist.size();
  long maxnormp = I2long(quadprimes[nap-1].norm());
  limitnorm = maxnormp;

  Quad nu = N.gen();
  Quad a=f->a, b=f->b, c=f->c, d=f->d; // matrix entries: [a,b;c*nu,d] is in Gamma_0(nu)
  if (debug)
    {
      cout << "Matrix [a,b;c*n,d] with a="<<a<<", b="<<b<<", c="<<c<<", d="<<d<<" and matdot="<<f->matdot<<endl;
    }
  double rootdisc = sqrt((double)(I2long(Quad::absdisc)));

  factor = 4*PI / (realnorm(nu*c)*rootdisc);
  bigcomplex ncrD =  to_bigcomplex(nu*c)*to_bigfloat(rootdisc);
  theta1 = to_bigcomplex(-d) / ncrD;
  theta2 = to_bigcomplex( a) / ncrD;

  //Initialize sum and use n=1 term:
  sum = 0;
  use(Quad(1),1);

  //add terms, one prime at a time:
  for(int i=0; i<nap; i++)
    {
      Quad p = quadprimes[i];
      if(debug>1)
        cout << "p= " << p << ",\ta_p = " << aplist[i];
      add(p,i,aplist[i],1);
      if (debug>1)
        cout << "\t\tSum = " << sum << endl;
    }
  period = abs(sum * factor / f->matdot);
}

void period_direct::use(const Quad& n, int an)
{
  double rn = realnorm(n);
  sum += double(an) * K1(rn*factor) * psi_factor(n) / rn;
}

// NB the code here is identical to period_via_lf1chi::add()

void period_direct::add(const Quad& n, int pindex, int y, int z)
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


