// FILE QUADS.CC

#define testbezout    // define this to turn on self-verification of bezout

#include <iostream>

#include "intprocs.h"
#include "quads.h"
#include "primes.h"
#include "geometry.h"

//Declare static data members of class Quad:

long Quad::d;
INT Quad::disc;
INT Quad::absdisc;
vector<INT> Quad::prime_disc_factors;
vector<INT> Quad::all_disc_factors;
long Quad::t;
INT Quad::n;
char Quad::name;
long Quad::maxnorm;
int Quad::nunits;
int Quad::is_Euclidean;
int Quad::class_number;
int Quad::class_group_2_rank;
Quad Quad::w;
Quad Quad::zero;
Quad Quad::one;

//Primes
vector<Quad> quadprimes;  //Initialised by initquadprimes, see below
long nquadprimes;         //The number of them.
vector<Quad> quadunits, squareunits;
Quad fundunit;

vector<long> euclidean_fields = {1,2,3,7,11};

// fields for which geometry is defined (all 305 with |disc|<1000 *except* 246 (984) and 249 (996), which include all h<=3)

vector<long> valid_fields = {1, 2, 3, 5, 6, 7,
                             10, 11, 13, 14, 15, 17, 19, 21, 22, 23, 26, 29, 30, 31, 33, 34, 35, 37, 38, 39, 41, 42, 43, 46, 47, 51, 53, 55,
                             57, 58, 59, 61, 62, 65, 66, 67, 69, 70, 71, 73, 74, 77, 78, 79, 82, 83, 85, 86, 87, 89, 91, 93, 94, 95, 97,
                             101, 102, 103, 105, 106, 107, 109, 110, 111, 113, 114, 115, 118, 119, 122, 123, 127, 129, 130, 131, 133, 134,
                             137, 138, 139, 141, 142, 143, 145, 146, 149, 151, 154, 155, 157, 158, 159, 161, 163, 165, 166, 167, 170, 173,
                             174, 177, 178, 179, 181, 182, 183, 185, 186, 187, 190, 191, 193, 194, 195, 197, 199,
                             201, 202, 203, 205, 206, 209, 210, 211, 213, 214, 215, 217, 218, 219, 221, 222, 223, 226, 227, 229, 230, 231,
                             233, 235, 237, 238, 239, 241, 246, 247, 249, 251, 255, 259, 263, 267, 271, 283, 287, 291, 295,299,
                             303, 307, 311, 319, 323, 327, 331, 335, 339, 347, 355, 359, 367, 371, 379, 383, 391, 395, 399,
                             403, 407, 411, 415, 419, 427, 431, 435, 439, 443, 447, 451, 455, 463, 467, 471, 479, 483, 487, 491, 499,
                             503, 511, 515, 519, 523, 527, 535, 543, 547, 551, 555, 559, 563, 571, 579, 583, 587, 591, 595, 599,
                             607, 611, 615, 619, 623, 627, 631, 635, 643, 647, 651, 655, 659, 663, 667, 671, 679, 683, 687, 691, 695, 699,
                             703, 707, 715, 719, 723, 727, 731, 739, 743, 751, 755, 759, 763, 767, 771, 779, 787, 791, 795, 799,
                             803, 807, 811, 815, 823, 827, 831, 835, 839, 843, 851, 859, 863, 871, 879, 883, 887, 895, 899,
                             903, 907, 911, 915, 919, 923, 935, 939, 943, 947, 951, 955, 959, 967, 971, 979, 983, 987, 991, 995};

vector<long> class_number_one_fields   = {1,2,3,7,11,19,43,67,163};                                     // 9
vector<long> class_number_two_fields   = {5,6,10,13,15,22,35,37,51,58,91,115,123,187,235,267,403,427};  // 18
vector<long> class_number_three_fields = {23,31,59,83,107,139,211,283,307,331,379,499,547,643,883,907}; // 16
vector<long> class_number_four_fields  = {14,17,21,30,33,34,39,42,46,55,57,70,73,78,82,85,93,97,102,130,
                                         133,142,155,177,190,193,195,203,219,253,259,291,323,355,435,483,
                                         555,595,627,667,715,723,763,795,955,1003,1027,1227,1243,1387,1411,
                                         1435,1507, 1555};                                                 // 54
vector<long> class_number_five_fields = {47,79,103,127,131,179,227,347,443,523,571,619,683,691,739,787,947,
                                         1051,1123,1723,1747,1867,2203,2347, 2683};                        // 25

vector<long> valid_field_discs(long max_disc)
{
  vector<long> discs;
  for(auto di = valid_fields.begin(); di!=valid_fields.end(); ++di)
    {
      long d = *di; // positive, square-free
      long D = (d%4==3? d : 4*d);
      if ((max_disc==0) || (D<=max_disc))
        discs.push_back(D);
    }
  ::sort(discs.begin(), discs.end());
  return discs;
}

int check_field(long d, vector<long> fields)
{
  return (std::find(fields.begin(), fields.end(), d) != fields.end());
}

// declaration of "extern" functions declared in quads.h:
Quad (*mult)(const Quad& a, const Quad& b);
Quad (*qdivi)(const Quad& a, INT c);
int (*pos)(const Quad& a);
Quad (*quadgcd)(const Quad& aa, const Quad& bb);
Quad (*quadbezout)(const Quad& aa, const Quad& bb, Quad& xx, Quad& yy);
Quad (*quadconj)(const Quad& a);

Quad mult0(const Quad& a, const Quad& b)
{
  return Quad(a.r*b.r - Quad::n*a.i*b.i, a.r*b.i + a.i*b.r, a.nm*b.nm);
}

Quad mult1(const Quad& a, const Quad& b)
{
  return Quad(a.r*b.r-Quad::n*a.i*b.i, a.r*b.i+a.i*b.r+a.i*b.i, a.nm*b.nm);
}

Quad qdivi0(const Quad& a, INT c) // c>0,    // used when t=0
{
  Quad ans;
  if (c>0)
    {
      ans.r = rounded_division(a.r,c);
      ans.i = rounded_division(a.i,c);
    }
  else
    {
      ans.r = rounded_division(-a.r,-c);
      ans.i = rounded_division(-a.i,-c);
    }
  ans.setnorm();
  return ans;
}

Quad qdivi1(const Quad& a, INT c) // used when t=1
{
  Quad ans;
  if (c>0)
    {
      ans.i = rounded_division(a.i,c);
      ans.r = rounded_division(2*a.r+a.i-c*ans.i,2*c);
    }
  else
    {
      ans.i = rounded_division(-a.i,-c);
      ans.r = rounded_division(-2*a.r-a.i+c*ans.i,-2*c);
    }
  ans.setnorm();
  return ans;
}

// static function (one for the class, not per instance)
void Quad::field(long dd, long max)
{
  // Clear these in case this is not the first field run

  quadunits.clear();
  squareunits.clear();
  quadprimes.clear();
  prime_disc_factors.clear();
  all_disc_factors.clear();
  class_group.clear();
  class_group_2_torsion.clear();
  class_group_2_cotorsion.clear();
  class_group_2_torsion_gens.clear();
  class_group_2_cotorsion_gens.clear();
  Qideal_lists::clear();

  d = squarefree_part(dd);
  if (d!=dd)
    cout << "Replacing d = " << dd << " with " << d << endl;
  is_Euclidean = check_field(d, euclidean_fields);
  class_number = 0;
  if (check_field(d, class_number_one_fields))
    class_number=1;
  else
    if (check_field(d, class_number_two_fields))
      class_number=2;
    else
      if (check_field(d, class_number_three_fields))
        class_number=3;
      else
        if (check_field(d, class_number_four_fields))
          class_number=4;
        else
          if (check_field(d, class_number_five_fields))
            class_number=5;

  // else class number is set in fill_class_group()

  INT odd(d&1?d:d/2); // odd part of d

  switch (d%4)
    {
    case 1:
      {
        t=0; absdisc=4*d; disc=-absdisc; n=d;
        prime_disc_factors.push_back(INT(-4));
        quadconj=&quadconj0;
        mult=&mult0; qdivi=&qdivi0;
        break;
      }
    case 2:
      {
        t=0; absdisc=4*d; disc=-absdisc; n=d;
        prime_disc_factors.push_back(INT(odd%4==1 ? -8 : 8));
        quadconj=&quadconj0;
        mult=&mult0; qdivi=&qdivi0;
        break;
      }
    case 3:
    default:
      {
        t=1; absdisc=d; disc=-d;   n=(d+1)/4;
        quadconj=&quadconj1;
        mult=&mult1; qdivi=&qdivi1;
      }
    }
  vector<INT> pp = pdivs(odd);
  for (auto pi=pp.begin(); pi!=pp.end(); ++pi)
    {
      INT p = *pi;
      prime_disc_factors.push_back((p%4==1?p:-p));
    }

  INT i0(0), i1(1);
  w = Quad(i0, i1, n);
  zero = Quad(i0,i0, i0);
  one = Quad(i1,i0, i1);

  switch (d) {
  case 1:  pos=&pos13; name='i'; nunits=4; fundunit=w; break;
  case 2:  pos=&pos2;  name='t'; nunits=2; fundunit=-one; break;
  case 3:  pos=&pos13; name='w'; nunits=6; fundunit=w; break;
  default: pos=&pos2;  name='a'; nunits=2; fundunit=-one;
  }

  quadgcd=&quadgcd_psea;
  quadbezout=&quadbezout_psea;

  int i;

  quadunits.push_back(one);
  quadunits.push_back(fundunit);
  for(i=2; i<nunits; i++)
    quadunits.push_back(fundunit*quadunits[i-1]);
  for(i=0; 2*i<nunits; i++)
    squareunits.push_back(quadunits[2*i]);
  maxnorm=max;
  if(class_number==1)
    initquadprimes();
  Quadprimes::init(max);
  if (check_field(d)) // test whether d is in the list of valid_fields
                      // for which we have geometry set up
    setup_geometry();

  fill_class_group();

  int n2r = class_group_2_rank; // which was set in fill_class_group()
  int nchi = (1<<n2r);
  for(int chi_index=0; chi_index<nchi; chi_index++)
    {
      INT D(1);
      for (int i=0; i<n2r; i++)
        if (bit(chi_index,i)==1)
          D *= prime_disc_factors[i];
      if (D>1) // negate positive ones except D=1 itself
        D = Quad::disc/D;
      all_disc_factors.push_back(D);
    }
}

// static function (one for the class, not per instance)
void Quad::displayfield(ostream& s, int info2)
{s<<"Q(sqrt("<<-d<<"))\tdiscriminant = "<<disc;
 s<<"\tmin poly("<<name<<") = "<<name<<"^2"; if(t) s<<"-"<<name; 
 s<<"+"<<n<<".\n";
 switch (d) {
 case 1: case 2: case 3: case 7: case 11:       // Euclidean
   cout << "Euclidean" << endl;
   break;
 case 19: case 43: case 67: case 163:           // Non-Euclidean
   cout << "Non-Euclidean, class number 1" << endl;
   break;
 default:
   cout << "Class number " << class_number << endl;
   cout << "Ideal class group representatives: " << class_group << endl;
   if (info2)
     {
       cout << "2-rank of class group = " << class_group_2_rank << endl;
       if (class_group_2_rank>0)
         {
           cout << " discriminant = ";
           for (auto di=prime_disc_factors.begin(); di!=prime_disc_factors.end(); ++di)
             {
               if (di!=prime_disc_factors.begin()) cout<<"*";
               cout<<"("<<(*di)<<")";
             }
           cout<<endl;
           cout << " 2-torsion   generators " << class_group_2_torsion_gens
                << ", 2-torsion   elements "<< class_group_2_torsion<<endl;
           cout << " 2-cotorsion generators " << class_group_2_cotorsion_gens
                << ", 2-cotorsion elements "<< class_group_2_cotorsion<<endl;
         }
     }
 }
 if (class_number==1)
   s<<nquadprimes<<" primes initialised, max norm = " << maxnorm << endl;
}

int Quad::chi(INT p)
{
  return (p==2? (d%4==3? (d%8==3? -1: +1): 0):  legendre(disc,p));
}

Quad makepos(const Quad& a)
{Quad ans=a;
  while(!pos(ans)) {ans*=fundunit; }
 return ans;
}

istream& operator>>(istream& s, Quad& x)
{
  s >> x.r >> x.i;
  x.setnorm();
  return s;
}

void Quad::setnorm()
{
  nm = r*r + n*i*i;
  if (t) {nm += r*i;};
  assert (nm>=0);
}

Quad Quad::operator/ (const Quad& b) const
{
  if (!::is_zero(b.i))
    return qdivi(mult(*this,quadconj(b)), b.nm);
  else
    return qdivi(*this,b.r);
}

void Quad::operator/=(const Quad& b)
{
  if (!::is_zero(b.i))
    *this=qdivi(mult(*this,quadconj(b)), b.nm);
  else
    *this=qdivi(*this,b.r);
}

// binary index i of [I] in C/C^2 (0 <= i < 2^{2-rank}).
//
// "binary" means that the bits of I give the coordinates w.r.t. the
// basis class_group_2_cotorsion_gens
int Quad::ideal_class_mod_squares(const Qideal& I)
{
  return find_ideal_class_mod_squares(I, class_group_2_cotorsion);
}

// image (+1/-1) of I under i'th quadratic character (0 <= i < 2^{2-rank})
int Quad::unramified_character(int i, const Qideal& I)
{
  return (i & ideal_class_mod_squares(I) ? -1 : +1); // bitwise &
}

// Quad::Quad(const bigcomplex& z)
// {bigfloat x=real(z), y=imag(z);
//  if(d>1) y/=sqrt(to_bigfloat(d));
//  if(d>2) {x-=y; y*=2.0;}
//  Iasb(r,x); Iasb(i,y);    //Rounded
//  //longify(x, r); longify(y, i);    //Rounded
//  setnorm();
// }

// Quad::operator bigcomplex() const
// {bigfloat x = to_bigfloat(r), y = to_bigfloat(i);
//  if(d>2) {y/=2.0; x+=y;}
//  if(d>1) y*=sqrt(to_bigfloat(d));
//  return bigcomplex(x,y);
// }

int div(const Quad& a, const Quad& b)
{
 if (a.nm==0) return (b.nm==0);
 if (b.nm==0) return 1;
 if (b.nm%a.nm!=0) return 0;
 Quad c = b*quadconj(a);
 return (((c.r)%a.nm)==0) && (((c.i)%a.nm)==0);
}

int ndiv(const Quad& a, const Quad& b)
{
 return !div(a,b);
}

int val(const Quad& factor, const Quad& number)
{
  if ((number.nm==0) || (factor.nm<=1))
    {
      cout << "Error in val(): factor = "<<factor<< " should not be a unit"<<endl;
      exit(1);
    }
  int e = 0; Quad n=number, f=factor, nf;
  while (nf=n/f, f*nf==n) {e++; n=nf;}
  return e;
}

vector<Quad> residues(const Quad& a)
{
  INT norma = a.norm(), m = gcd(a.re(), a.im());
  INT rednorma = (norma/m)/m;
  vector<Quad> ans;
  for(int j=0; j<m*rednorma; j++)
    {
      INT J(j);
      for(int k=0; k<m; k++)
        ans.push_back(Quad(J,INT(k))%a);
    }
  return ans;
}

Quad::operator string() const
{
  ostringstream s;
  if (i==0)
    s<<r;
  else
    {
      if (r==0)
        {
          if (i==1) ;
          else if (i==-1) s << "-";
          else s << i;
        }
      else
        {
          s<<r;
          if(i>0) s<<"+"; else s<<"-";
          if (abs(i)>1) s<<abs(i);
        }
      s<<Quad::name;
    }
  return s.str();
}

ostream& operator<<(ostream& s, const Quad& a)
{
  s << (string)a;
  return s;
}


//Functions for computing quad-primes, initializing the vector<Quad>
//quadprimes.  NB all primes are "pos" i.e. normalized w.r.t. units

void factorp0(long p, INT& a, INT& b, INT d)
// finds a,b s.t. a^2+d*b^2=0 (mod p)
{ int found=0;
  for (int ib=1; !found; ib++)
  {
    b = ib;
    INT a2 = p - d*b*b;
    a = a2.isqrt();
    found = (a*a == a2);
  }
}

void factorp1(long p, INT& a, INT& b, INT d)
// finds a,b s.t. a^2+a*b+((d+1)/4)*b^2=0 (mod p)
{ int found=0; long fourp = 4*p;
  for (int ib=1; !found; ib++)
  {
    b = ib;
    INT a2 = fourp -d*b*b;
    a = a2.isqrt();
    found = (a*a == a2);
  }
  a=(a-b)/2;
}

vector<Quad> Quad::primes_above(long p, int& sig)
{
  INT d(Quad::d), P(p);
  int t=Quad::t;
  INT a,b;  Quad pi, piconj;
  vector<Quad> list;
  sig = Quad::chi(P);
  //  cout<<"disc = "<<Quad::disc<<", p="<<p<<", chi(p)="<<sig<<endl;
  INT i0(0), i1(1), i2(2), i3(3);
  switch (sig) {
  case  0: // ramified
    pi =  (d==1 ? Quad(i1,i1, i2) :
           d==2 ? Quad(i0,i1, i2) :
           d==3 ? Quad(i1,i1, i3) :
           Quad(-i1,i2, d));
    list.push_back(pi);
    break;
  case -1: // inert
    pi = Quad(P,i0, P*P);
    list.push_back(pi);
    break;
  case 1: // split
    if(t==0) factorp0(p,a,b,d); else factorp1(p,a,b,d);
    pi = makepos(Quad(a,b, P));
    piconj = makepos(quadconj(pi));
    // We choose the ordering so the HNFs are [p,c,1], [p,c',1] with c<c'
    long c = posmod((a%p)*invmod(b,p),p);
    if (2*c<p)
      {
        list.push_back(pi);
        list.push_back(piconj);
      }
    else
      {
        list.push_back(piconj);
        list.push_back(pi);
      }
  }
  //cout << "primes_above("<<p<<") = " << list << endl;
  return list;
}

void Quad::initquadprimes()
{
  int sig;
  vector<Quad> list, list1, list2;
  vector<Quad>::iterator pi, alpha, beta;
  for (primevar pr; pr.ok()&&pr<maxnorm; pr++)
    { long p=pr;
      list = Quad::primes_above(p, sig);
      for (pi = list.begin(); pi!=list.end(); )
        switch (sig) {
        case  0:
          list1.push_back(*pi++);
          break;
        case -1:
          if(p*p<=maxnorm)
            list2.push_back(*pi);
          ++pi;
          break;
        case 1:
          {
            list1.push_back(*pi++);
            list1.push_back(*pi++);
          }
        }
    }

  // Now list1 contains the degree 1 primes and list2 the degree 2
  // primes, each ordered by norm; we merge these lists so that they
  // are still sorted by norm:

  alpha=list1.begin(); beta=list2.begin();
  while ((alpha!=list1.end()) && (beta!=list2.end()))
    if (quadnorm(*alpha)<quadnorm(*beta))
      quadprimes.push_back(*alpha++);
    else
      quadprimes.push_back(*beta++);

  while (alpha!=list1.end()) quadprimes.push_back(*alpha++);
  while ( beta!=list2.end()) quadprimes.push_back(*beta++);

  nquadprimes = quadprimes.size();
}

// Quad::fill_class_group() is implemented in qidloop.cc

Quad primdiv(const Quad& a)
{
  INT na=quadnorm(a);
  if (na<2) return Quad::zero;   // must return something!
  vector<Quad>::const_iterator pr;
  for (pr=quadprimes.begin(); pr!=quadprimes.end(); ++pr)
    {
      Quad p=*pr;
      if (div(p,a)) return p;
      INT np=quadnorm(p);
      if (np*np>na) return makepos(a);
    }
  cout<<"No prime divisor found for "<<a<<" so assuming prime!\n";
  return makepos(a);
}

vector<Quad> pdivs(const Quad& aa)
{ Quad a=aa; INT norma=quadnorm(a);
  vector<Quad> plist; // will hold prime factors
  if (norma<2) return plist;
  vector<Quad>::const_iterator pr;
  for (pr=quadprimes.begin(); pr!=quadprimes.end(); ++pr)
     { Quad p = *pr;
       if (div(p,a)) 
	 {
	   plist.push_back(p);
	   while (div(p,a)) a/=p; 
	   norma=quadnorm(a);	 
	   if (norma==1) return plist;
	 }
       else 
	 {
	   INT normp=quadnorm(p);
	   if (normp*normp>norma) 
	     {
	       plist.push_back(makepos(a));
	       return plist;
	     }
	 }
     }
  //In case of p-factors outside range, assume the cofactor is prime:
  if (quadnorm(a)>1)
    plist.push_back(makepos(a));  
  return plist;
}

vector<Quad> posdivs(const Quad& a)   // all "positive" divisors (up to units)
{
  vector<Quad> plist=pdivs(a); Quad p; 
  int e, nu = 1; int nd=nu;
  vector<int> elist;
  vector<Quad>::iterator pr;
  for (pr=plist.begin(); pr!=plist.end(); ++pr)
    {
      e=val(*pr,a); 
      elist.push_back(e);
      nd*=(1+e);
    }
  vector<Quad> dlist(nd);
  dlist[0]=Quad::one;
  nd=nu;
  vector<int>::iterator ei;
  for (pr=plist.begin(), ei=elist.begin(); pr!=plist.end(); ++pr, ++ei)
   {
     p = *pr;
     e = *ei;
     for (int j=0; j<e; j++)
       for (int k=0; k<nd; k++)
         dlist[nd*(j+1)+k] = makepos(p*dlist[nd*j+k]);
     nd*=(e+1);
   } 
 return dlist;
}

vector<Quad> alldivs(const Quad& a)       // all divisors
{
  vector<Quad> plist=pdivs(a); Quad p; 
  int e, nu = Quad::nunits; int nd=nu;
  vector<int> elist;
  vector<Quad>::iterator pr;
  for (pr=plist.begin(); pr!=plist.end(); ++pr)
   {
     e=val(*pr,a); 
     elist.push_back(e);
     nd*=(1+e);
   }
  vector<Quad> dlist(nd);
  for(int i=0; i<nu; i++) dlist[i]=quadunits[i];
  nd=nu;
  vector<int>::iterator ei;
  for (pr=plist.begin(), ei=elist.begin(); pr!=plist.end(); ++pr, ++ei)
   {
     p = *pr; 
     e = *ei;
     for (int j=0; j<e; j++)
       for (int k=0; k<nd; k++)
         dlist[nd*(j+1)+k] = p*dlist[nd*j+k];
     nd*=(e+1);
   } 
 return dlist;
}

vector<Quad> sqdivs(const Quad& a) // all divisors whose square divides a, up to +/-
{
  vector<Quad> plist=pdivs(a); Quad p;
  int e, nu = Quad::nunits/2; int nd=nu;
  vector<int> elist;
  vector<Quad>::iterator pr;
  for (pr=plist.begin(); pr!=plist.end(); ++pr)
   {
     e=val(*pr,a)/2; 
     elist.push_back(e);
     nd*=(1+e);
   }
  vector<Quad> dlist(nd);
  for(int i=0; i<nu; i++) dlist[i]=quadunits[i];
  nd=nu;
  vector<int>::iterator ei;
  for (pr=plist.begin(), ei=elist.begin(); pr!=plist.end(); ++pr, ++ei)
   {
     p = *pr; 
     e = *ei;
     for (int j=0; j<e; j++)
       for (int k=0; k<nd; k++)
         dlist[nd*(j+1)+k] = p*dlist[nd*j+k];
     nd*=(e+1);
   } 
 return dlist;
}

vector<Quad> sqfreedivs(const Quad& a)       // all square-free divisors
{
  vector<Quad> plist=pdivs(a); Quad p;
  int e, nu = 2; int nd=nu;
  vector<int> elist;
  vector<Quad>::iterator pr;
  for (pr=plist.begin(); pr!=plist.end(); ++pr)
   {
     e=1; 
     elist.push_back(e);
     nd*=(1+e);
   }
  vector<Quad> dlist(nd);
  for(int i=0; i<nu; i++) dlist[i]=quadunits[i];
  nd=nu;
  vector<int>::iterator ei;
  for (pr=plist.begin(), ei=elist.begin(); pr!=plist.end(); ++pr, ++ei)
   {
     p = *pr; 
     e = *ei;
     for (int j=0; j<e; j++)
       for (int k=0; k<nd; k++)
         dlist[nd*(j+1)+k] = p*dlist[nd*j+k];
     nd*=(e+1);
   } 
 return dlist;
}

Quad invmod(const Quad& a, const Quad& p)
{Quad x,y;
 Quad g=quadbezout(a,p,x,y);
 if (g==Quad::one) return x;
 else {cerr<<"invmod called with "<<a<<" and "<<p<<" -- not coprime!"<<endl;
   return Quad::zero;
   exit(1);
      }
}

int coprime(const Quad& a, const Quad& b) 
{
  Quad g = quadgcd(a,b); 
  return g==Quad::one;
}

int invertible(const Quad& a, const Quad& b, Quad& inverse)
{ Quad y; Quad g = quadbezout(a,b,inverse,y);
  return g==Quad::one;
}

//#define test_reduce

// returns n = round(true_real_part_of(alpha/beta)), so alpha-n*beta
// is reduced mod Z<beta>
INT nearest_INT_to_Quad_quotient ( const Quad& alpha, const Quad& beta)
{
  return rounded_division((alpha*quadconj(beta)).re(), beta.norm());
}

// reduction of gamma modulo Z<alpha,beta>
Quad reduce_mod_zbasis(const Quad& gamma, const Quad& alpha, const Quad& beta)
{
#ifdef test_reduce
  cout << "reduction of "<<gamma<<" mod <"<<alpha<<","<<beta<<">"<< flush;
#endif
  INT d = (quadconj(alpha)*beta).im();
  assert (d>0);
  INT x = rounded_division((-gamma*quadconj(beta)).im(), d);
  INT y = rounded_division(( gamma*quadconj(alpha)).im(), d);
  Quad ans = gamma - (x*alpha + y*beta);
#ifdef test_reduce
  cout << " is "<< ans << " (d="<<d<<", x="<<x<<", y="<<y<<")"<<endl;
  cout << " x*alpha =  "<< x*alpha << endl;
  cout << " y*beta =  "<< y*beta << endl;
  cout << " x*alpha+y*beta =  "<< x*alpha+y*beta << endl;
#endif
  INT gn = quadnorm(ans);
  vector<Quad> tests = {ans+alpha, ans-alpha, ans+beta, ans-beta,
                        ans-alpha-beta, ans-alpha+beta, ans+alpha-beta, ans+alpha+beta};
  for (vector<Quad>::const_iterator t=tests.begin(); t!=tests.end(); ++t)
    {
      if (gn>quadnorm(*t))
        {
#ifdef test_reduce
          cout<<"reduction "<<ans<<" has larger norm than shift "<<*t<<", switching"<<endl;
#endif
          ans = *t;
          gn = quadnorm(*t);
        }
    }
  return ans;
}

#undef test_reduce

// Replace alpha, beta by an SL2Z-equivalent pair with the same Z-span.
// The new alpha has the smallest norm.
void sl2z_reduce(Quad& alpha, Quad& beta)
{
#ifdef test_reduce
  cout<<"SL2Z-reducing ["<<alpha<<","<<beta<<"]..."<<endl;
#endif
  int s=1;
  Quad t;
  while (s)
    {
      s = 0; // will be set to 1 if anything changes
      INT n = nearest_INT_to_Quad_quotient(beta,alpha);
      if(n!=0)
        {
          s=1;
          beta -= n*alpha;
#ifdef test_reduce
          cout<<" -- shift by " << n <<": alpha="<<alpha<<", beta="<<beta<< endl;
#endif
        }
      if (beta.nm < alpha.nm)
        {
          s=1;
          t = -alpha; alpha = beta; beta = t;
#ifdef test_reduce
          cout<<" -- invert: alpha="<<alpha<<", beta="<<beta<< endl;
#endif
        }
    }
  // alpha is now the non-zero Quad of least norm in the lattice [alpha,beta]
  // and beta/alpha is in the closed fundamental region

  // We want to orient so that beta/alpha has positive imaginary part
  // e.g. so that the basis [1,w] is reduced
  // NB im(x+y*w)>0 iff y>0 for both t=0 and t=1 cases

  if ((quadconj(alpha)*beta).im() < 0)
    beta=-beta;

  // If we wanted to make sure that beta/alpha is unambiguously defined when on the boundary:
  // if (2*(quadconj(alpha)*beta).re() == -quadnorm(alpha))
  //   beta+=alpha;

#ifdef test_reduce
  cout<<"After reduction, alpha="<<alpha<<", beta="<<beta
      <<" with norms "<<quadnorm(alpha)<<", "<<quadnorm(beta)<<endl;
  assert (quadnorm(alpha)<=quadnorm(beta));
  //assert (nearest_INT_to_Quad_quotient(alpha,beta)==0);
  assert ((quadconj(alpha)*beta).im() > 0);
#endif
}

// Replace alpha, beta by an SL2Z-equivalent pair with the same Z-span.
// The new alpha has the smallest norm.
// U holds the unimodular transform (not implemented for FLINT INTs as the unimod class uses NTL ZZ)
// void sl2z_reduce(Quad& alpha, Quad& beta, unimod&U)
// {
//   Quad alpha0=alpha, beta0=beta;
//   INT U11, U12, U21, U22;
// #ifdef test_reduce
//   cout<<"SL2Z-reducing ["<<alpha<<","<<beta<<"]..."<<endl;
// #endif
//   U.reset(); // to the identity
//   int s=1;
//   Quad t;
//   while (s)
//     {
//       s = 0; // will be set to 1 if anything changes
//       INT n = nearest_INT_to_Quad_quotient(beta,alpha);
//       if(n!=0)
//         {
//           s=1;
//           beta -= n*alpha;
//           U.y_shift(INT(n));
// #ifdef test_reduce
//           cout<<" -- shift by " << n <<": alpha="<<alpha<<", beta="<<beta<< endl;
// #endif
//           U11 = I2long(U(1,1)); U12 = I2long(U(1,2));
//           U21 = I2long(U(2,1)); U22 = I2long(U(2,2));
//           assert (U11*alpha+U12*beta == alpha0);
//           assert (U21*alpha+U22*beta == beta0);
//         }
//       if (beta.nm < alpha.nm)
//         {
//           s=1;
//           t = -alpha; alpha = beta; beta = t;
//           U.invert();
// #ifdef test_reduce
//           cout<<" -- invert: alpha="<<alpha<<", beta="<<beta<< endl;
//           U11 = I2long(U(1,1)); U12 = I2long(U(1,2));
//           U21 = I2long(U(2,1)); U22 = I2long(U(2,2));
//           assert (U11*alpha+U12*beta == alpha0);
//           assert (U21*alpha+U22*beta == beta0);
// #endif
//         }
//     }
//   // alpha is now the non-zero Quad of least norm in the lattice [alpha,beta]
//   // and beta/alpha is in the closed fundamental region

//   // We want to orient so that beta/alpha has positive imaginary part
//   // e.g. so that the basis [1,w] is reduced
//   // NB im(x+y*w)>0 iff y>0 for both t=0 and t=1 cases

//   if ((quadconj(alpha)*beta).im() < 0)
//     beta=-beta;

//   // If we wanted to make sure that beta/alpha is unambiguously defined when on the boundary:
//   // if (2*(quadconj(alpha)*beta).re() == -quadnorm(alpha))
//   //   beta+=alpha;

// #ifdef test_reduce
//   cout<<"After reduction by U="<<U<<", alpha="<<alpha<<", beta="<<beta
//       <<" with norms "<<quadnorm(alpha)<<", "<<quadnorm(beta)<<endl;
//   U11 = I2long(U(1,1)); U12 = I2long(U(1,2));
//   U21 = I2long(U(2,1)); U22 = I2long(U(2,2));
//   assert (U11*alpha+U12*beta == alpha0);
//   assert (U21*alpha+U22*beta == beta0);
//   assert (quadnorm(alpha)<=quadnorm(beta));
//   assert (nearest_INT_to_Quad_quotient(alpha,beta)==0);
//   assert ((quadconj(alpha)*beta).im() > 0);
// #endif
// }


// HNF of ideal (alpha) as a triple [a c d] where [a,c+d*w] is a Z-basis with
//
// a,d>0; c>=0
// N=a*d = Norm(alpha)
// d|a and d|b
// 0 <=c < a

vector<INT> HNF(const Quad& alpha)
{
  INT N = alpha.norm(), xa = alpha.re(), ya = alpha.im(), u, v;
  INT g = bezout(xa,ya,u,v);  // g=u*xa+v*ya=gcd(xa,ya)
  INT x = xa/g, y = ya/g;
  // Now the HNF is g*[a, b+w] for some b mod a=N/g
  INT a = N/(g*g);
  INT b = posmod(((v-(Quad::t)*u)*x - (Quad::n)*u*y), a);
  vector<INT> ans;
  ans.push_back(a*g);
  ans.push_back(b*g);
  ans.push_back(g);
  return ans;
}

// Ideal label: formed from the Norm and HNF of the ideal (alpha)
// (subject to change!)

string old_ideal_label(const Quad& alpha)  // returns label of ideal (alpha)
{
  vector<INT>H = HNF(alpha);
  stringstream s;
  s << "[" << alpha.norm() << "," << H[1] << "," << H[2] << "]";
  return s.str();
}

string ideal_label(const Quad& alpha)  // returns label of ideal (alpha)
{
  vector<INT>H = HNF(alpha);
  stringstream s;
  s << alpha.norm() << "." << H[1] << "." << H[2];
  return s.str();
}

string field_label() // returns field label, e.g. '2.0.4.1'
{
  stringstream s;
  s << "2.0." << Quad::absdisc << ".1";
  return s.str();
}

int are_associate(const Quad& a, const Quad& b)
{
  if(a.is_zero()) return (b.is_zero());
  if(a.norm() != b.norm()) return 0;
  vector<Quad>::const_iterator eps;
  for(eps=quadunits.begin(); eps!=quadunits.end(); ++eps)
    if(a*(*eps)==b) return 1;
  return 0;
}

int is_ideal_Galois_stable(const Quad& a)
{
  Quad b(quadconj(a));
  return are_associate(a, b);
}

string ideal_code(const Quad& N) // string code for a (principal) ideal N
{
  vector<INT>H = HNF(N);
  stringstream s;
  s << quadnorm(N) << "." << H[1] << "." << H[2];
  return s.str();
}

vector<int> makechitable(const Quad& lambda, const vector<Quad>& reslist)
{
  vector<int> chi;
  if(reslist.size()==1)
    chi.push_back(1);
  else
    {
      vector<Quad>::const_iterator r=reslist.begin();
      while(r!=reslist.end())
	chi.push_back(squaremod(*r++,lambda,reslist));
    }
  return chi;
}

// brute force test whether a is a square of some element of reslist, mod m

int squaremod(const Quad& a, const Quad& m, const vector<Quad>& reslist)
{
  if (div(m,a)) return 0;
  vector<Quad>::const_iterator r=reslist.begin();
  while(r!=reslist.end())
    {
      Quad res=*r++;
      if(div(m,res*res-a)) return +1;
    }
  return -1;
}

// bigfloat gauss(const Quad& m, const vector<Quad>& reslist)
// {
// //cout<<"Computing g(chi) for lambda = " << m << endl;
//   bigfloat ans1(0); //double ans2=0;
//   bigfloat rootdisc;
//   rootdisc = sqrt(to_bigfloat(Quad::absdisc));
//   bigcomplex lrd = bigcomplex(m)*bigcomplex(to_bigfloat(0), rootdisc);
//   vector<Quad>::const_iterator r=reslist.begin();
//   while(r!=reslist.end())
//   {
//     Quad res = *r++;
//       bigfloat term1 = squaremod(res,m,reslist)*psif(bigcomplex(Quad(res))/lrd);
//       ans1+=term1;      //    ans2+=term2;
//   }
//   return ans1;
// }

// convert a binary vector into a discriminant dividing Quad::disc
INT discchar(vector<int> c)
{
  INT D(1);
  for(int i=0; i<(int)c.size(); i++)
    if (c[i])
      D *= Quad::prime_disc_factors[i];
  return D;
}

// convert a discriminant dividing Quad::disc into a binary vector

vector<int> chardisc(INT D)
{
  vector<int> ans;
  for(auto di=Quad::prime_disc_factors.begin(); di!=Quad::prime_disc_factors.end(); ++di)
    ans.push_back(div_disc(*di,D));
  // cout<<"prime_disc_factors: "<< Quad::prime_disc_factors << endl;
  // cout<<"chardisc("<<D<<") = "<<ans<<endl;
  return ans;
}

// all_disc_factors modulo D mod squares, i.e. factoring out D.  D
// should be in all_disc_factors.  Returns a list of half the length
// unless D=1.
vector<INT> disc_factors_mod_D(const INT& D)
{
  if ((D==ONE)||(D==ZERO))
    return Quad::all_disc_factors;
  vector<int> Dv = chardisc(D);
  // cout<<"D="<<D<<" Dv="<<Dv<<endl;
  int i = std::find(Dv.begin(), Dv.end(), 0) - Dv.begin();
  int j = std::find(Dv.begin(), Dv.end(), 1) - Dv.begin();
  vector<INT> ans;
  for( auto Di=Quad::all_disc_factors.begin(); Di!=Quad::all_disc_factors.end(); ++Di)
    {
      // cout<<"Di="<<(*Di)<<":"<<chardisc(*Di)<<endl;
      if((chardisc(*Di)[i]==0) && (chardisc(*Di)[j]==0))
        ans.push_back(*Di);
    }
  // cout<<"All discs "<<Quad::all_disc_factors<<" mod "<<D<<": "<<ans<<endl;
  return ans;
}

#if 0
// there follow 4 "flavours" of findminQuad returning different amounts of data

vector<INT> findminquadcoeffs(const Quad&al, const Quad&be, Quad& alpha, Quad& beta)
{ alpha=al;
  beta=be;
  INT i0(0), i1(1);
  vector<INT> c = {i1, i0}; // alpha = c[0]*al + c[1]*be  always
  vector<INT> d = {i0, i1}; // beta  = d[0]*al + d[1]*be  always

  while (1)
    {
      INT n = nearest_INT_to_Quad_quotient(alpha,beta);
      alpha -= n*beta;
      c[0]  -= n*d[0];
      c[1]  -= n*d[1];
      Quad temp = alpha; alpha = -beta; beta = temp;
      INT t = c[0]; c[0] = -d[0]; d[0] = t;
            t = c[1]; c[1] = -d[1]; d[1] = t;
      if (quadnorm(beta) >= quadnorm(alpha))
        {
          // alpha is now the non-zero Quad of least norm in the lattice [al,be]
          // beta is another Quad such that [al,be]=[alpha,beta]
#ifdef testbezout
          if (alpha != c[0]*al + c[1]*be)
            {
              cerr << "Error in findminquadscoeffs!" << endl;
              cerr << "[al,be] = ["<<al<<","<<be<<"]"<<endl;
              cerr << "c0,c1 = "<<c[0]<<","<<c[1]<<endl;
              cerr << "alpha = "<<alpha<< " not equal to "<<c[0]*al + c[1]*be<<endl;
            }
          exit(1);
#endif
          return c;
        }
    }
}

vector<INT> findminquadcoeffs(const Quad& alpha, const Quad& beta, Quad& gen0)
{ Quad gen1;
  return findminquadcoeffs(alpha,beta,gen0,gen1);
}

//#define DEBUG_FINDMINQUAD

void findminquad(const Quad&al, const Quad&be, Quad& alpha, Quad& beta)
// same as findminquadscoeffs but don't need coeffs
{ alpha=al;
  beta=be;
#ifdef  DEBUG_FINDMINQUAD
  unimod U; Quad alpha1=alpha, beta1=beta;
  sl2z_reduce(alpha1, beta1, U);
  cout<<"In findminquad with alpha="<<alpha<<" (norm "<<quadnorm(alpha) <<"), beta="
     <<beta<<" (norm "<<quadnorm(beta)<<")"<<endl;
#endif
  while (1)
    {
      alpha -= nearest_INT_to_Quad_quotient(alpha,beta)*beta;
      Quad temp = alpha; alpha = -beta; beta = temp;
#ifdef  DEBUG_FINDMINQUAD
      cout<<" ... alpha="<<alpha<<" (norm "<<quadnorm(alpha)<<"), beta="
          <<beta<<" (norm "<<quadnorm(beta)<<")"<<endl;
#endif
      if (quadnorm(beta) >= quadnorm(alpha))
        {
          // alpha is now the non-zero Quad of least norm in the lattice [al,be]
          // beta is another Quad such that [al,be]=[alpha,beta]
#ifdef  DEBUG_FINDMINQUAD
          cout<<" ... on return, alpha = "<<alpha<<" has minimal norm "<<quadnorm(alpha)<<endl;
#endif
          return;
        }
    }
}
 
void findminquad(const Quad& alpha, const Quad& beta, Quad& gen0)
// same as findminquadcoeffs but don't need coeffs
{ Quad gen1;
  findminquad(alpha,beta,gen0,gen1);
}

#endif
