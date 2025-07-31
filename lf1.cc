// lf1.cc : class period_via_lfchi for integrating newforms

#include <values.h>
#include "eclib/interface.h"
#include "eclib/kbessel.h"
#include "lf1.h"

const double twopi = 2.0*PI;
const double eps = 1.0e-20;   // ?? mindouble;

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
  Quad p, pn;
  long maxnorm=limitnorm/I2long(n.norm());
  if ( y!=0 ) {use (n,y); istart=0;}
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
   ratio(f->loverp)
{
  double rootdisc = sqrt((double)(I2long(Quad::absdisc))), modlambda = realnorm(lambda);
  factor = 4*PI/(rootdisc*sqrt(realnorm(N.gen()))*modlambda);
  lambdares = residues(lambda);
  chitable = makechitable(lambda, lambdares);
  if(debug)
  {
    cout << "rootdisc = sqrt{|D_K|} = " << rootdisc <<endl;
    cout << "factor = (4*pi) / (rootdisc * sqrt{|n|} * |lambda|) = " <<factor <<endl;
    cout << "Table of chi mod lambda (lambda = "<<lambda<<")\n";
    cout << "residue \t chi \n";
    for(unsigned int ires=0; ires<lambdares.size(); ires++)
      cout << lambdares[ires] << "\t" << chitable[ires] << endl;
  }

//Extract the aplist from f->aplist which has all the coefficients at
//all primes in standard order:
  aplist = f->aplist;

  if(debug)
    {
      cout<<"nap = "<<nap<<endl;
    }
  long maxnormp = I2long(quadprimes[nap-1].norm());
  limitnorm = maxnormp;   //  for debugging only
  if(debug)
  {
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
      if(debug)
        cout << "p= " << p << ",\ta_p = " << aplist[i];
      add(p,i,aplist[i],1);
      if(debug)
        cout << "\t\tSum = " << sum << endl;
      if (p.norm()>=limitnorm)
        break;
    }

//Scale to get values of L(f_chi,1) and the period:
  lf1chivalue = 2*factor*sum;
  period=lf1chivalue*modlambda/double(ratio);
  if(debug)
  {
    cout<<"lf1chivalue = "<<lf1chivalue<<", ratio = "<<ratio;
    cout<<", period = "<<period<<endl;
  }
}



/*
Complex epi(double x)
{ double theta = 2*PI*x;
  return Complex(cos(theta),sin(theta));
}

void periods_direct::use(int n, int an)
{ double dn = (double)n, dan = (double)an;
  sum += (dan/dn) * exp(dn*efactor) * (epi(dn*theta1)-epi(dn*theta2));
}

void periods_direct::add(int n, int pindex, int y, int z)
{
  int p,pn,ip,istart=pindex;
  int triv = (y==0);
  if ( ! triv ) { use (n,y); istart=0; }

  for(pn=n*(p=quadprimes(ip=istart));
      (ip<=pindex) && (pn<=limit);
      pn=n*(p=quadprimes(++ip)))
    { int x = y * aplist[ip];
      if ( (ip==pindex)  && (N%p)) { x -=  p*z; }
      add(pn,ip,x,y);
    }
}


periods_direct::periods_direct(h1newform*f)
:N(level::modulus),nap(f->aplist.length)
{
//Sort out the aplist, since f_>aplist has the W_q eigs first:
  aplist=longlist(nap);
  int ip = level::npdivs;
  longvar p(quadprimes);
  for( ; p.ok(); p++)
    { int i = p.index;
      int j = level::plist.locate(p);  // = -1 if p is good
      if(j+1) // then p is j'th bad prime
        {
          if(div(p*p,N))
            {
              aplist[i]=0;
//cout << "p = " << p << "\ta_p = " << aplist[i] << endl;
            }
          else
            {
              aplist[i] = - (f->aplist)[j];
//cout << "p = " << p << "\ta_p = " << aplist[i] << endl;
            }
        }
        else  // p is good
        {
          aplist[i] = (f->aplist)[ip++];
//cout << "p = " << p << "\ta_p = " << aplist[i] << endl;
        }

    }
//  cout << "After rearranging, aplist = " << aplist << endl;
  long maxp = prime(nap);
  limit = maxp;
  double a=f->a,b=f->b,c=f->c,d=f->d;
  if (c<0) { a=-a;b=-b;c=-c;d=-d;}
  sum = 0.0;
  double cmodrecip =  1.0 / (c*N);
  theta1 = -d * cmodrecip;
  theta2 =  a * cmodrecip;
  efactor = -2 * PI * cmodrecip;
//n=1 term:
  use(1,1);
//add terms, one prime at a time:
  for (p.init(quadprimes); p.ok()&&(p<=limit); p++)
    { int ip=p.index;
//      cout << "p= " << p << ",\ta_p = " << aplist[ip];
      add(p,ip,aplist[ip],1);
//      cout << "\t\tSum now " << sum << endl;
    }
  double u = real(sum)/(f->dotplus);
  double v = imag(sum)/(f->dotminus);
  //and then the periods:
  periods = new Complex[2];
  switch (f->type) {
  case 2: periods[0]=Complex(u,0); periods[1]=Complex(0,v); break;
        case 1: periods[0]=Complex(2*u,0); periods[1]=Complex(u,v); break;
        }
}

*/

