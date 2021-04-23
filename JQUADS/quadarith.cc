// FILE quadarith.cc

#undef testbezout

//#ifndef testbezout
//#define testbezout    // define this to turn on self-verification of bezout
//#endif

#undef use_obs_bezout // define this to include obsolete code in compilation

#include <iostream>
//#include <builtin.h>

#include "arith.h"
#include "marith.h"
#include "quadarith.h"

Quad (*quadgcd)(const Quad& aa, const Quad& bb);
Quad (*quadbezout)(const Quad& aa, const Quad& bb, Quad& xx, Quad& yy);

/*
Quad::Quad(const complex& z)
{
  double x=z.real(), y=z.imag();
  if (Field::d>1) y/=Field::rootd;
  if (Field::t==1) {x-=y; y*=2.0;}
  r=(long)(x+0.5); i=(long)(y+0.5);    //Rounded
}

Quad::operator complex() const
{
  double x=r, y=i;
  if (Field::t) {y/=2.0; x+=y;}
  if (Field::d>1) y*=Field::rootd;
  return complex(x,y);
}
*/

Quad Quad::operator/ (const Quad& b) const
{
  if (b.i==0) {return (*this)/(b.r); }
  else
    {
      // (prodr,prodi)=(*this)*conj(b)
      Quad c=b.conj();
      double ar=(double)r, ai=(double)i, br=(double)c.r, bi=(double)c.i;
      double prodr = ar * br - ai * bi * Field::n;
      double prodi = ar * bi + ai * br;
      if (Field::t) prodi += ai * bi;
      
      double nb = b.dnorm();
      if (Field::t)
	{ long ansi = ldroundover(prodi,nb);
	  long ansr = ldroundover(2*prodr + prodi - nb*ansi, 2*nb);
	  return Quad(ansr,ansi);
	}
      else
	{ return Quad(ldroundover(prodr,nb), ldroundover(prodi,nb)); }
    }
}

double Quad::dnorm() const
// rtns the norm as a double
{
  double ar=(double)r, ai=(double)i;
  double ans = ar*ar + ai*ai*Field::n;
  if (Field::t)
    {
      ans +=  ar*ai;
    }
  return ans;
}

int Quad::norm(long&n) const
// rtns 1 if the norm is too big to fit in a long, 0 if okay with n=norm
{
  double temp=dnorm();
  n = (long)temp;
  return (temp>MAXLONG);
}

int Quad::small_integer_multiple(long&s, Quad&aux) const
{
  long g=gcd(real(),imag());
  Quad a=Quad(real()/g,imag()/g);
  aux=a.conj();
  double ds=g*a.dnorm();
  int errflag=(ds>MAXLONG);
  if (errflag) s=MAXLONG; else s=(long)ds;    // MAXLONG so min will choose ok
  return errflag;
}

int div(const long& a, const Quad& b)
{
    if (b==0) return 1;     // anything divides zero ...
    if (a==0) return 0;     // ... and zero divides nothing else
    return ((real(b)%a)==0) && ((imag(b)%a)==0);
}

int ndiv(const long& a, const Quad& b)
{ return !div(a,b); }

int Quad::div(const Quad& b) const
{
  if (b==0) return 1;     // anything divides zero.
  if (imag()==0) return ::div(real(),b);
  
  Quad temp=b/(*this);       // crude but works !
  return (temp*(*this)==b);
}
// neater might be to inline the following idea using doubles
//    long na = a.norm();
//    Quad c = b*(a.conj());
//    return ((real(c)%na)==0) && ((imag(c)%na)==0);

int div(const Quad& a, const Quad& b)
{ return a.div(b); }

int ndiv(const Quad& a, const Quad& b)
{ return !(a.div(b)); }

long val(const Quad& factor, const Quad& number)
{ 
  if ( (number==0) || (factor.norm()<=1) )
    {
      cerr<<"Warning: 9999 returned in Quad version of val"<<endl;
      return 9999;
    }
  long e = 0; Quad n=number, f=factor, nf;
  while (nf=n/f, f*nf==n) {e++; n=nf;}
  return e;
}

Quadlist residues(const Quad& a)
{
    long norma = a.norm(), m = gcd(a.real(),a.imag());
    long rednorma = (norma/m)/m;
    Quadlist ans(norma); long count=0;
    for(long j=0; j<m*rednorma; j++)
	for(long k=0; k<m; k++) {
	    Quad res = Quad(j,k)%a; //cout<<"res "<<count<<" = "<<res;
	    ans[count++]= res ;
	}
    return ans;
}

ostream& operator<<(ostream& s, const Quad& a)
{
  long  i = a.i, r = a.r;
  if (i==0) s<<r;
  else
    {
      if (r==0)
	{
	  if (i==-1) s << "-";
	  else if (i!=1) s << i;
	}
      else
	{
	  s<<r; 
	  if(i>0) s<<"+"; else s<<"-";
	  if (abs(i)>1) s<<abs(i);
	}
      s<<Field::name;
    }
  return s;
}

char* to_string(const Quad& a)
{
  ostringstream ans;
  ans << a;  
  return (char*)ans.str().c_str();
}

Quad quadgcd1(const Quad& aa, const Quad& bb)  
//Only works for Euclidean fields!
{
  Quad a=aa,b=bb,c,q;
  while (b!=0) {
    q = a/b; c = a-q*b; a = b; b = c;
    if (b.norm() >= a.norm())
      {cout<<"Error in quadgcd1 --norm not reduced!\n";
       exit(1);}
  }
  return a.pos_assoc();
}

Quad
quadbezout1(const Quad& alpha, const Quad& beta, Quad& coeff1, Quad& coeff2)
{Quad a=alpha,b=beta,c,x=0,oldx=1,newx,y=1,oldy=0,newy,q,g;
 while (b!=0)
   { q = a/b; 
     c    = a    - q*b; a    = b; b = c;
     newx = oldx - q*x; oldx = x; x = newx;
     newy = oldy - q*y; oldy = y; y = newy;
   }
 coeff1=oldx; coeff2=oldy; g=a;
 //Next two lines try to get coeff1,coeff2 as small as possible
 if(beta!=0)
   {
     coeff1 = coeff1%(beta/g);           //reduced
     coeff2 = (g-coeff1*alpha)/beta;     //should be exact
   }
 while (!g.is_pos())
   { g*=Field::fundunit; coeff1*=Field::fundunit; coeff2*=Field::fundunit;
   }
#ifdef testbezout
 //CHECK:
 if (div(g,alpha) && div(g,beta) && (g==coeff1*alpha+coeff2*beta)) {;}  //OK
 else 
   {cerr<<"Error in quadbezout1!"<<endl;
    cerr<<"alpha = "<<alpha<<endl;
    cerr<<"beta  = "<<beta<<endl;
    cerr<<"coeff1 = "<<coeff1<<endl;
    cerr<<"coeff2 = "<<coeff2<<endl;
    cerr<<"g   = "<<g<<endl;
  }
#endif
 return g;
}

Quad invmod(const Quad& a, const Quad& p)
{
    Quad x,y;
    Quad g=quadbezout(a,p,x,y);
    if (g==1) return x;
    else {
	cerr<<"invmod called with "<<a<<" and "<<p<<" -- not coprime!"<<endl;
	return 0;
    }
}

int coprime(const Quad& a, const Quad& b) 
{
    Quad g=quadgcd(a,b); 
    return g==1;
}

int invertible(const Quad& a, const Quad& b, Quad& inverse)
{ 
    Quad y; Quad g = quadbezout(a,b,inverse,y);
    return g==1;
}

//functions needed for non-euclidean fields to compute bezout/quadgcd

long vecbezout(long n, long* a, long* c) 
//returns g = content(a) = a.c
{
    long x=1,ai=0,ci=1,g=0;
//cout<<"In vecbezout with a="<<a<<", c="<<c<<endl;
    for(long i=0; i<n; i++) {
//cout<<"...calling bezout with g="<<g<<", a[i]="<<a[i]<<", x="<<x<<", ci="<<ci<<endl;
	g=bezout(g,a[i],x,ci);
//cout<<"...returns x="<<x<<", ci="<<ci<<endl;
	c[i]=ci;
	for(long j=0; j<i; j++) c[j]*=x;
    }
    return g;
}

long xmodvecbezout(long s, long n, long* a, long* c) 
//returns g = content(a) = a.c, with c reduced mod s
{
    long x=1,ai=0,ci=1,g=0;
    for(long i=0; i<n; i++)
      {
	g=bezout(g,a[i],x,ci);
	c[i] = xmod(ci,s);
	for(long j=0; j<i; j++) c[j]=xmodmul(c[j],x,s);
      }
    return g;
}
 
long vecgcd(long n, long* a)
//returns g = content(a)
{
    long g=0;
    for(long i=0; (i<n)&(g!=1); i++) g=gcd(g,a[i]);
    return g;
}
 
long dot(long n, long* a, long* c) 
//returns g = a.c
{
    long g=0;
    for(long i=0; i<n; i++) g+=a[i]*c[i];
    return g;
}

long xmoddot(long s, long n, long* a, long* c) 
//returns g = a.c, with g reduced mod s
{
    long g=0;
    for(long i=0; i<n; i++) g = xmod( g + xmodmul(a[i],c[i],s), s);
    return g;
}

double ddot(long n, long* a, long* c) 
//returns double g = a.c
{
    double g=0;
    for(long i=0; i<n; i++) g+=((double)a[i])*((double)c[i]);
    return g;
}

long special_solve_for_last_coeff(long d, long n, long*f, long*x, long s)
// returns ans s.t.  s*ans + (f dot x) equals d
{
  double tmp = d - ddot(n, f, x);
  if ( fmod(tmp, s) != 0 )
    {
      cerr << "Error: inexact division in special_solve_for_last_coeff!"<<endl;
    }
  tmp/=s;
  if (fabs(tmp)>MAXLONG)
    {
      cerr << "Error: overflow in special_solve_for_last_coeff!"<< endl;
    }
  long ans=(long)tmp;
  return ans;
}

// for debugging
void output_nlongs(ostream&op, long n, long*a, long new_line=1)
{
  op << a[0];
  for(long i=0; i<n; i++) op << " " << a[i];
  if (new_line==1) op << endl;
}

void specialfindzbasiscoeffsmod(long n, long* first, long* second,
				long* basis, long* x, long* y)
// Given: a (2xn) matrix a (of special form) with rows "first", "second"
// Returns: "basis" [e1,e2,f1] such that the cols of
//          (e1 f1)  are a Z-basis for the cols of a.
//          (e2  0)
//          "x","y": cols of coeffs (nx2), s.t. a.x=e, a.y=f
// Special form means that "first" is (... s 0) and "second" is (... 0 s)
// where s is a non-zero "help" integer.
//
// The values of x[n-2], x[n-1], y[n-2], y[n-1] may need to be treated
// with caution, since they were not calculated modulo s.
{
  if ((n<2)|| (second[n-2] != 0) || (first[n-1]!=0) ||
      (first[n-2]==0) || (first[n-2] != second[n-1]))
    { cerr << "Invalid call to specialfindzbasiscoeffsmod!" << endl;
      exit(1);
    }
  long s=first[n-2];
  long* u=new long[n];
  long* newfirst =new long[n];
  long* newsecond=new long[n];                    //temps

  long e2=xmodvecbezout(s,n,second,x);
//Now second.x=e2 ...
// ... except that x[n-2] is arbitrary (currently zero) and x[n-1] is wrong !

  long e1=xmoddot(s,n,first,x);                  //dot product
//Now  first.x=e1,  for suitable choice of x[n-2] (depending on correct
// value of x[n-1]

//That is, [e1,e2] is the x-combination of the data, with e2=gcd(second) ...
// ... except that the currently held values of x[n-2], x[n-1] are incorrect !

//newfirst = first-e1*(second/e2);
// (It was to stop this step overflowing that we chose a small e1.)
  for(long i=0; i<n; i++)
    {
      newsecond[i]=second[i]/e2;
      if (i==n-2) newfirst[i]=s;  // prevent reduction of s to zero !
      else newfirst[i]=xmod(first[i]-xmodmul(e1,newsecond[i],s),s);
    }
// Since newfirst[n-2]=s, f1 comes out OK even though rest of row is mod s

  long f1 = xmodvecbezout(s,n,newfirst,u);
//Now newfirst.u=f1, except that u[n-2] is wrong!

  basis[0] = e1;   basis[1] = e2;   basis[2] = f1;

  long t = xmoddot(s,n,u,newsecond);
// correct despite u[n-2] being wrong - because newsecond[n-2]=0 !

//  y = u - ((second/e2)*u)*x;
  for(long i=0; i<n-2; i++) y[i]=xmod(u[i]-xmodmul(t,x[i],s),s);

// So far, so good.
// Now try computing x[n-2], x[n-1], y[n-2], y[n-1] !!

// The trouble with the naive method is that it overflows too easily.
// x[n-2] = (e1 - dot(n-2, first,x))/s;
// x[n-1] = (e2 - dot(n-2,second,x))/s;

// y[n-2] = (f1 - dot(n-2, first,y))/s;
// y[n-1] = ( 0 - dot(n-2,second,y))/s;

  x[n-2] = special_solve_for_last_coeff(e1, n-2, first, x, s);
  x[n-1] = special_solve_for_last_coeff(e2, n-2,second, x, s);
  y[n-2] = special_solve_for_last_coeff(f1, n-2, first, y, s);
  y[n-1] = special_solve_for_last_coeff( 0, n-2,second, y, s);

if (0)
    { cerr <<"Reporting from specialfindzbasiscoeffsmod!"  <<endl;
      cerr <<"modulus s = "<<s<<endl;

      cerr <<"first = "; output_nlongs(cerr, n, first);

      cerr << "second= "; output_nlongs(cerr, n, second);

      cerr << "e1 = " << e1 << endl;
      cerr << "e2 = " << e2 << endl;
      cerr << "f1 = " << f1 << endl;

      cerr <<" x = "; output_nlongs(cerr, n, x);

      cerr <<" y = "; output_nlongs(cerr, n, y);

      cerr <<"Warning: x[n-2], x[n-1], y[n-2], y[n-1] are less robust!"<<endl;
    }

  delete[] newsecond;  delete[] newfirst;  delete[] u;
}



// The following procecdure is crucial in ideal arithmetic !!!
//
void findzbasis(long n, long* first, long* second, long* basis)
// much like old_specialfindzbasiscoeffsmod
// differences: n needn't be 4; j needn't be 2; coeffs are not needed
//
// Rewritten version - works modulo s.  Okay provided xmodmul can cope.
// Requires the existence of a "help" vector, and reports an error o/w.
//
{
  // look for a "help" vector modulo which to work
  //
  // look for j with second[j]=0 and first[j] (non-zero) minimal
  long j=-1;
  for (long i=0; i<n; i++)
    { if ((second[i]==0)&&(first[i]!=0))
	if (j==-1) j=i;
	else if (abs(first[i])<abs(first[j])) j=i;
    }

  if (j==-1)
    { cerr << "Error: findzbasis called without a help vector!"<<endl;
      exit(1);
    }
  long s = first[j];

  long* x=new long[n];   
  long* newfirst=new long[n];                    //temps

  long e2=xmodvecbezout(s,n,second,x);  // x only computed mod s
  // e2 = gcd(second) = second.(true x), where x= (true x)mod s
  // Note that x[j] = (true x)[j] = 0

  // (true e1) = first.(true x) might be huge, but the help vector allows
  // us to adjust (true x)[j] in order to reduce (true e1) modulo s.

  long e1=xmoddot(s,n,first,x);         //dot product, reduced mod s

  //  Now [e1,e2] is the (adjusted true x)-combination of the data.
  //  x= (true x)mod s  --  except that x[j] is incorrect !
  //  We cannot determine the true x[j] - but it is not needed.

  //newfirst = first-e1*(second/e2);
  for(long i=0; i<n; i++)
    {
      if (i==j) newfirst[i]=s;  // prevent reduction of s to zero !
      else newfirst[i]=xmod(first[i]-xmodmul(e1,second[i]/e2,s),s);
    }
  // It was to stop this step overflowing that we chose a small e1.
  // Since newfirst[j]=s, f1 comes out OK even though rest of row is mod s

  long f1 = vecgcd(n,newfirst);
  basis[0] = e1;   basis[1] = e2;   basis[2] = f1;
  delete[] newfirst;
  delete[] x;
}

long nearest_long_to_Quad_quotient ( const Quad& alpha, const Quad& beta)
// returns round(true_real_part_of(alpha/beta))  
{
  // (prodr,prodi)=(alpha)*conj(beta)
  Quad c=beta.conj();
  double ar=(double)alpha.real(), ai=(double)alpha.imag();
  double br=(double)c.real(), bi=(double)c.imag();
  double prodr = ar * br - ai * bi * Field::n;
  double prodi = ar * bi + ai * br;
  if (Field::t) prodi += ai * bi;
	
  double nb = beta.dnorm();
  if (Field::t)
    {
      return ldroundover(2*prodr + prodi, 2*nb);
    }
  else
    {
      return ldroundover(prodr,nb);
    }
}

// there follow 4 "flavours" of findminQuad returning different amounts of data

long* findminquadscoeffs(const Quad&al, const Quad&be, Quad& beta, Quad& alpha)
{ alpha=al;
  beta=be;
  long n;
  double normalpha,normbeta=beta.dnorm();
  Quad temp;
  long* c = new long[2];
  long* d = new long[2];  
  long v;
  c[0] = 1; c[1] = 0; d[0] = 0; d[1] = 1;
  while (
	 n = nearest_long_to_Quad_quotient(alpha,beta),
         //= round(true_real_part_of(alpha/beta))
         alpha -= n*beta,
//       d     -= n*c,
         d[0]     -= n*c[0],
         d[1]     -= n*c[1],
         normalpha = alpha.dnorm(),
         (normbeta > normalpha)
         )
    {
      temp = alpha; alpha = -beta; beta = temp;
      normbeta = normalpha;
      v    = d[0];     d[0]     = -c[0];    c[0]    = v;
      v    = d[1];     d[1]     = -c[1];    c[1]    = v;
    }
// beta is now the non-zero Quad of least norm in the lattice [al,be]
// alpha is another Quad such that [al,be]=[alpha,beta]
  delete[] d;
#ifdef testbezout
  if (beta != c[1]*al + c[0]*be)
    {
      cerr << "Error in findminquadscoeffs!" << endl;
      cerr << "[al,be] = ["<<al<<","<<be<<"]"<<endl;
      cerr << "c1,c0 = "<<c[1]<<","<<c[0]<<endl;
      cerr << "beta = "<<beta<< " not equal to "<<c[1]*al + c[0]*be<<endl;
    }
#endif
  return c;
}

long* findminquadcoeffs(const Quad& alpha, const Quad& beta, Quad& gen0)
{ Quad gen1;
  return findminquadscoeffs(alpha,beta,gen0,gen1);
}

void findminquads(const Quad&al, const Quad&be, Quad& beta, Quad& alpha)
// same as findminquadscoeffs but don't need coeffs
{ alpha=al;
  beta=be;
  long n;
  double normalpha,normbeta=beta.dnorm();
  Quad temp;
  while (
	 n = nearest_long_to_Quad_quotient(alpha,beta),
         //= round(true_real_part_of(alpha/beta))
         alpha -= n*beta,
         normalpha = alpha.dnorm(),
	 (normbeta > normalpha)
	 )
    {
      temp = alpha; alpha = -beta; beta = temp;
      normbeta = normalpha;
    }
// beta is now the non-zero Quad of least norm in the lattice [al,be]
// alpha is another Quad such that [al,be]=[alpha,beta]
}
 
void findminquad(const Quad& alpha, const Quad& beta, Quad& gen0)
// same as findminquadcoeffs but don't need coeffs
{ Quad gen1;
  findminquads(alpha,beta,gen0,gen1);
}

#ifdef use_obs_bezout
#include "obs_bez.cc"
#endif

// special follows!!!
Quad special_bezout(const Quad& alpha, const long& beta,
		    Quad& coeff1, Quad& coeff2)
// rtns non-zero g of minl norm in <alpha,beta>, st
//    g=(coeff1)alpha + (coeff2)beta
// in ptic, if <alpha,beta> is principal, then <alpha,beta>=<g>.
{
  Quad g;  
  if (div(beta,alpha)) { g=beta; coeff1=0; coeff2=1;}
  else if (alpha.div(beta)) { g=alpha; coeff1=1; coeff2=0;}
  else                                        // in ptic: alpha, beta non-zero
    {
      long* rv = new long[4];
      long* iv = new long[4];
      long* basis = new long[3];
      long* x = new long[4];
      long* y = new long[4];   
      long* z = new long[4];
      rv[0] = alpha.real();       iv[0] = alpha.imag();
      rv[1] = -Field::n*iv[0];    iv[1] = rv[0] + Field::t*iv[0];
      rv[2] = beta;               iv[2] = 0;
      rv[3] = 0;                  iv[3] = beta;

      // in the following call, x,y are only computed modulo s
      // x[0],x[1],y[0] and y[1] may be used subsequently
      // x[2],x[3],y[2] and y[3] should also be quite robust now
      specialfindzbasiscoeffsmod(4,rv,iv,basis,x,y);

      Quad al(basis[0],basis[1]);
      Quad be=basis[2];

      long s=beta;
      // in the following call, coeff need only be computed modulo s
      // (for now, work out coeff exactly)
      long* coeff = findminquadcoeffs(al,be,g);

      // make sure there will be no overflow in computation of z[i]
      for(long i=0; i<2; i++) coeff[i]=xmod(coeff[i],s);

      for(long i=0; i<2; i++)
	{ z[i] = xmod(xmodmul(coeff[0],y[i],s) + xmodmul(coeff[1],x[i],s),s); }
      coeff1 = Quad(z[0],z[1]);


// Now try to get coeff1,2 as small as possible.
// (Dividing by beta is okay, since the case beta==0 was dealt with above.)

// Since we have worked mod s, coeff1 is ALREADY reduced mod s.
// The line computing coeff2 tends to overflow.  So we use a trick.
//    coeff1 = coeff1%beta;            //reduced
//    coeff2 = (g-coeff1*alpha)/beta;  //should be exact

      z[2] = special_solve_for_last_coeff(g.real(), 2, rv, z, s);
      z[3] = special_solve_for_last_coeff(g.imag(), 2, iv, z, s);
      coeff2 = Quad(z[2],z[3]);

      delete[] coeff; delete[] z; delete[] y; delete[] x;
      delete[] basis; delete[] iv; delete[] rv;

      // for debugging the improvement
      // to trigger this event, try e.g. hecketest, level of norm 16, t_npnp
      if (alpha == Quad(-73045,146111))
	{
	  cerr << "For debug info set breakpoint at this statement!" << endl;
	}

    }

  while (!g.is_pos())
    {
      g*=Field::fundunit; coeff1*=Field::fundunit; coeff2*=Field::fundunit;
    } 

/*
#ifdef testbezout
//CHECK:
  if ( (Field::not_princ() || (g.div(alpha) && g.div(beta)) ) && (g==coeff1*alpha+coeff2*beta)) {;}  //OK
  else 
    {cerr<<"Error in special_bezout!"<<endl;
     cerr<<"alpha = "<<alpha<<endl;
     cerr<<"beta  = "<<beta<<endl;
     cerr<<"coeff1 = "<<coeff1<<endl;
     cerr<<"coeff2 = "<<coeff2<<endl;
     cerr<<"a*c1 + b*c2 = "<< coeff1*alpha+coeff2*beta << endl;
     cerr<<"g   = "<<g<<endl;
     }
#endif
*/

  return g;
} // end of special_bezout

// new quadbezout2 follows!!!
Quad quadbezout2(const Quad& q1, const Quad& q2, Quad& coeff1, Quad& coeff2)
// rtns non-zero g of minl norm in <q1,q2>, st  g=(coeff1)q1 + (coeff2)q2
// in ptic, if <q1,q2> is principal, then <q1,q2>=<g>.
{
  Quad g;  
  if (q2.div(q1)) { g=q2; coeff1=0; coeff2=1;}
  else if (q1.div(q2)) { g=q1; coeff1=1; coeff2=0;}
  else                                        // in ptic: q1,q2 non-zero
    {
      Quad aux1,aux2,alpha;
      long beta;
      long s1,s2;

      int swapflag;
      int errflag1=q1.small_integer_multiple(s1,aux1);
      int errflag2=q2.small_integer_multiple(s2,aux2);
      if (errflag1)
	if (errflag2)
	  { cerr <<"Error: quadbezout2 can't find small help integer!\n";}
	else
	  { swapflag=0; }
      else
	if (errflag2)
	  { swapflag=1; }
	else
	  { swapflag = ((aux1==1)&&((aux2!=1)||(s1<s2))) ||
	      ((aux1!=1)&&(aux2!=1)&&(s1<s2)); }
      if (swapflag)
	{ alpha=q2*aux1; beta=s1; }
      else
	{ alpha=q1*aux2; beta=s2; }
      Quad gaux = special_bezout(alpha,beta,coeff1,coeff2);
      if (swapflag)
	{ Quad temp=coeff1; coeff1=coeff2; coeff2=temp; }
      g = coeff1*q1 + coeff2*q2;
      Quad check=g;
      if (swapflag) check*=aux1; else check*=aux2;
      if (check!=gaux)
	{ cerr << "aux error in quadbezout2!"<<endl; }
    }

  while (!g.is_pos())
    {
      g*=Field::fundunit; coeff1*=Field::fundunit; coeff2*=Field::fundunit;
    } 
#ifdef testbezout
//CHECK:
  if ( (Field::not_princ() || (g.div(q1) && g.div(q2)) ) && (g==coeff1*q1+coeff2*q2)) {;}  //OK
  else 
    {cerr<<"Error in quadbezout2!"<<endl;
     cerr<<"alpha = "<<q1<<endl;
     cerr<<"beta  = "<<q2<<endl;
     cerr<<"coeff1 = "<<coeff1<<endl;
     cerr<<"coeff2 = "<<coeff2<<endl;
     cerr<<"a*c1 + b*c2 = "<< coeff1*q1+coeff2*q2 << endl;
     cerr<<"g   = "<<g<<endl;
     }
#endif
  return g;
}

// rewritten version - includes a "help" vector
//
Quad quadgcd2(const Quad& alpha, const Quad& beta)
// Much the same as quadbezout2 except don't need coeff1, coeff2
{
  Quad g;
  if (beta.div(alpha)) { g = beta; }
  else if (alpha.div(beta)) { g= alpha; }
  else
    {
      Quad aux1, aux2;
      long s=0,s1,s2;
      if (!alpha.small_integer_multiple(s1,aux1)) s=gcd(s,s1);
      if (!beta.small_integer_multiple(s2,aux2)) s=gcd(s,s2);
      if (s==0)
	{ cerr <<"Error: quadgcd2 can't find small enough help integer!\n";}
      long* rv = new long[5];
      long* iv = new long[5];
      long* basis = new long[3];
      rv[0] = alpha.real();       iv[0] = alpha.imag();   
      rv[1] = -Field::n*iv[0];    iv[1] = rv[0] + Field::t*iv[0];
      rv[2] = beta.real();        iv[2] = beta.imag();
      rv[3] = -Field::n*iv[2];    iv[3] = rv[2] + Field::t*iv[2];
      rv[4] = s;                  iv[4] = 0;
      findzbasis(5,rv,iv,basis);
      Quad al(basis[0],basis[1]);
      Quad be=basis[2];
      findminquad(al,be,g);
      delete[] basis; delete[] iv; delete[] rv;
    }
  g = g.pos_assoc();
#ifdef testbezout
  //CHECK:
  if (Field::not_princ() || (div(g,alpha) && div(g,beta)) ) {;}  //OK
  else 
    {cerr<<"Error in quadgcd2!"<<endl;
     cerr<<"alpha = "<<alpha<<endl;
     cerr<<"beta  = "<<beta<<endl;
     cerr<<"g   = "<<g<<endl;
   }
#endif
  return g;
}

 
//-------------------------------------------------------------------------
// N.B.  The point of the following ghastly functions is that the
// built-in gcc division truncates towards 0, while we need
// rounding, with a consistent behaviour for halves (they go up here). 
//-------------------------------------------------------------------------

long roundover_away_from_zero(long a, long b)  // halves round away from zero
{ //cout<<"roundover("<<a<<","<<b<<") = ";

  // wlog the divisor is positive
  if (b<0) { a=-a; b=-b; }

  // first obtain the remainder in range   -b/2 < rem <= b/2
  long c = mod (a,b) ;

  // with no adjustment below, halves round toward minus infinity

  // if a>0 adjust to range  -b/2 <= rem < b/2
  // this adjusts so positive halves round towards plus infinity
  // the result is rounding of halves away from zero
  if ((a>0)&&(c+c==b)) { c-=b; }

  // // alternatively
  // // if a<0 adjust to range -b/2 <= rem < b/2
  // // this adjusts so negative halves round towards plus infinity
  // // the result would be rounding of halves towards zero
  // // if ((a<0)&&(c+c==b)) { c-=b; }

  long ans= (a-c)/b;

  //cout<<ans<<endl;
  return ans;
}

long roundover(long a, long b)      // my emulation of John's version
{
  // wlog the divisor is positive
  // combine with negation because default rounding of halves is down
  if (b<0) { b=-b; } else { a=-a;}

  // first obtain the remainder in range   -b/2 < rem <= b/2
  long c = mod (a,b) ;

  long ans= (c-a)/b;     // second negation to cancel first
  return ans;
}

// John's old version :  halves rounded "up" i.e. towards plus infinity
// overflows more readily than mine because of the step c=a+a+b;
long old_roundover_by_JEC(long aa, long bb)
{ //cout<<"roundover("<<aa<<","<<bb<<") = ";
    long ans;
    if(aa%bb==0) ans = aa/bb;   //exact case always OK
    else { 
	long a=aa,b=bb;
	if(b<0){a=-a; b=-b;}
	long c = a+a+b, b2=b+b; 
	ans = c/b2;
	if((c<0)&&(c%b2))ans--;
    }
    //cout<<ans<<endl;
    return ans;
}

long ldroundover(double a, double b)
{
  if (b<0) { b=-b; } else { a=-a; }

  // need to obtain the remainder in range   -b/2 < c <= b/2
  // crude code, but should work for since 0.5 is exact in binary
  double c = fmod(a,b);
  double bd = b/2;
  if (c>bd) c-=b; else if (c<=-bd) c+=b;
//  cout << "DEBUG: a="<<a<<"  b="<<b<<"  mymod(a,b)="<<c<<endl;

  double ans = (c-a)/b;
  return (long)ans;
}

bigint floor(const bigint &a, const bigint&b )
// returns floor of quotient a/b
{
  if (is_zero(b))
    {
      cerr << "Error: floor called with b=0\n" ;
      long e,f=0; e/=f;
    }
  bigint ans, r;
  int dummy = ::divides(a,b,ans,r);
  if ( (sign(a) == -(sign(b))) && (! is_zero(r)) )
    {
      ans-=1;
// next line commented out since not required for return value
//    r+=b;
    }
// next line commented out since for diagnostics only
//  cout<<"Debug:  " << a << " = " << ans << " * " << b << " + " << r << endl;
  return ans;
}

// END OF FILE quadarith.cc
