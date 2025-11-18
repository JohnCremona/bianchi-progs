// FILE QUADS.CC

#define testbezout    // define this to turn on self-verification of bezout

#include <iostream>

#include "intprocs.h"
#include "quads.h"
#include "primes.h"
#include "geometry.h"
#include "swan.h"

//Declare static data members of class Quad (except SD which is declared in geometry.cc)

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
vector<Quad> Quad::shifts_by_one;
int Quad::geometry_initialised;

//Primes
vector<Quad> quadprimes;  //Initialised by initquadprimes, see below
long nquadprimes;         //The number of them.
vector<Quad> quadunits, squareunits;
Quad fundunit;

vector<long> euclidean_fields = {1,2,3,7,11};

// fields for which geometry is defined (all 305 with |disc|<1000 *except* 246 (984) and 249 (996), which include all h<=3)

vector<long> valid_fields = {1, 2, 3, 5, 6, 7,
                             10, 11, 13, 14, 15, 17, 19, 21, 22, 23, 26, 29, 30, 31, 33, 34, 35, 37, 38, 39,
                             41, 42, 43, 46, 47, 51, 53, 55, 57, 58, 59, 61, 62, 65, 66, 67, 69, 70, 71, 73,
                             74, 77, 78, 79, 82, 83, 85, 86, 87, 89, 91, 93, 94, 95, 97,
                             101, 102, 103, 105, 106, 107, 109, 110, 111, 113, 114, 115, 118, 119, 122, 123,
                             127, 129, 130, 131, 133, 134, 137, 138, 139, 141, 142, 143, 145, 146, 149, 151,
                             154, 155, 157, 158, 159, 161, 163, 165, 166, 167, 170, 173, 174, 177, 178, 179,
                             181, 182, 183, 185, 186, 187, 190, 191, 193, 194, 195, 197, 199,
                             201, 202, 203, 205, 206, 209, 210, 211, 213, 214, 215, 217, 218, 219, 221, 222,
                             223, 226, 227, 229, 230, 231, 233, 235, 237, 238, 239, 241, 246, 247, 249, 251,
                             253, 254, 255, 257, 258, 259, 262, 263, 265, 266, 269, 267, 271, 273, 274, 277,
                             278, 281, 282, 283, 285, 286, 287, 290, 291, 293, 295, 298, 299,
                             301, 302, 303, 305, 307, 309, 310, 311, 314, 313, 317, 318, 319, 322, 321, 323,
                             326, 327, 329, 330, 331, 334, 335, 337, 339, 341, 345, 346, 347, 349, 353, 354,
                             355, 357, 358, 359, 362, 365, 366, 367, 370, 371, 373, 374, 377, 379, 381, 382,
                             383, 385, 386, 389, 390, 391, 393, 394, 395, 397, 398, 399,
                             401, 402, 403, 406, 407, 409, 410, 411, 413, 415, 417, 418, 419, 421, 422, 426,
                             427, 429, 430, 431, 433, 434, 435, 437, 438, 439, 442, 443, 445, 446, 447, 449,
                             451, 453, 454, 455, 457, 458, 461, 462, 463, 465, 466, 467, 469, 470, 471, 473,
                             474, 478, 479, 481, 482, 483, 485, 487, 489, 491, 493, 494, 497, 498, 499,
                             501, 502, 503, 505, 506, 509, 510, 511, 514, 515, 517, 518, 519, 521, 523, 526,
                             527, 530, 533, 534, 535, 537, 538, 541, 542, 543, 545, 546, 547, 551, 553, 554,
                             555, 557, 559, 561, 562, 563, 565, 566, 569, 570, 571, 573, 574, 577, 579, 581,
                             582, 583, 586, 587, 589, 590, 591, 593, 595, 597, 598, 599,
                             601, 602, 606, 607, 609, 610, 611, 613, 614, 615, 617, 618, 619, 622, 623, 627,
                             631, 635, 643, 647, 651, 655, 659, 663, 667, 671, 679, 683, 687, 691, 695, 699,
                             703, 707, 715, 719, 723, 727, 731, 739, 743, 751, 755, 759, 763, 767, 771, 779,
                             787, 791, 795, 799,
                             803, 807, 811, 815, 823, 827, 831, 835, 839, 843, 851, 859, 863, 871, 879, 883,
                             887, 895, 899,
                             903, 907, 911, 915, 919, 923, 935, 939, 943, 947, 951, 955, 959, 967, 971, 979,
                             983, 987, 991, 995,
                             1003, 1007, 1011, 1015, 1019, 1023, 1027, 1031, 1039, 1043, 1047, 1051, 1055, 1059,
                             1063, 1067, 1079, 1087, 1091, 1095, 1099,
                             1103, 1111, 1115, 1119, 1123, 1131, 1135, 1139, 1147, 1151, 1155, 1159, 1163, 1167,
                             1171, 1187, 1191, 1195, 1199,
                             1203, 1207, 1211, 1219, 1223, 1227, 1231, 1235, 1239, 1243, 1247, 1255, 1259, 1263,
                             1267, 1271, 1279, 1283, 1291, 1295, 1299,
                             1303, 1307, 1311, 1315, 1319, 1327, 1335, 1339, 1343, 1347, 1351, 1355, 1363, 1367,
                             1371, 1379, 1383, 1387, 1391, 1399,
                             1403, 1407, 1411, 1415, 1419, 1423, 1427, 1435, 1439, 1443, 1447, 1451, 1455, 1459,
                             1463, 1471, 1479, 1483, 1487, 1491, 1495, 1499,
                             1507, 1511, 1515, 1523, 1527, 1531, 1535, 1543, 1547, 1551, 1555, 1559, 1563, 1567,
                             1571, 1579, 1583, 1591, 1595, 1599,
                             1603, 1607, 1615, 1619, 1623, 1627, 1631, 1635, 1639, 1643, 1651, 1655, 1659, 1663,
                             1667, 1671, 1679, 1687, 1691, 1695, 1699,
                             1703, 1707, 1711, 1723, 1727, 1731, 1735, 1739, 1743, 1751, 1747, 1759, 1763, 1767,
                             1771, 1779, 1783, 1787, 1795, 1799,
                             1803, 1807, 1819, 1811, 1823, 1831, 1835, 1839, 1843, 1847, 1851, 1855, 1867, 1871,
                             1879, 1883, 1887, 1891, 1895,
                             1903, 1907, 1915, 1919, 1923, 1927, 1931, 1939, 1943, 1947, 1951, 1955, 1959, 1963,
                             1967, 1979, 1983, 1987, 1991, 1995, 1999,
                             2003, 2011, 2015, 2019, 2027, 2031, 2035, 2039, 2047, 2051, 2055, 2059, 2063, 2067,
                             2071, 2083, 2087, 2091, 2095, 2099,
                             2103, 2111, 2119, 2123, 2127, 2131, 2135, 2139, 2143, 2147, 2155, 2159, 2163, 2167,
                             2171, 2179, 2183, 2191, 2195, 2199,
                             2203, 2207, 2211, 2215, 2219, 2227, 2231, 2235, 2239, 2243, 2247, 2251, 2255, 2263,
                             2267, 2271, 2279, 2283, 2287, 2291,
                             2307, 2311, 2315, 2319, 2323, 2327, 2335, 2339, 2343, 2347, 2351, 2355, 2359, 2363,
                             2371, 2379, 2383, 2387, 2391, 2395, 2399,
                             2407, 2411, 2415, 2419, 2423, 2427, 2431, 2435, 2443, 2447, 2451, 2455, 2459, 2463,
                             2467, 2471, 2479, 2483, 2487, 2491, 2495,
                             2683};

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
  for(auto d : valid_fields)
    {
      long D = (d%4==3? d : 4*d);
      if ((max_disc==0) || (D<=max_disc))
        discs.push_back(D);
    }
  ::sort(discs.begin(), discs.end());
  return discs;
}

// declaration of "extern" functions declared in quads.h:
Quad (*mult)(const Quad& a, const Quad& b);
Quad (*qdivi)(const Quad& a, const INT& c);
int (*pos)(const Quad& a);
Quad (*quadgcd)(const Quad& aa, const Quad& bb);
Quad (*quadbezout)(const Quad& aa, const Quad& bb, Quad& xx, Quad& yy);
Quad (*quadconj)(const Quad& a);

Quad mult0(const Quad& a, const Quad& b)
{
  return Quad(fmms(a.r,b.r, Quad::n*a.i,b.i), fmma(a.r,b.i, a.i,b.r), a.nm*b.nm);
}

Quad mult1(const Quad& a, const Quad& b)
{
  INT x = a.i*b.i;
  return Quad(fmms(a.r,b.r, Quad::n,x), fmma(a.r,b.i, a.i,b.r) + x, a.nm*b.nm);
}

// compute a*b+c*d
Quad mma(const Quad& a, const Quad& b, const Quad& c, const Quad& d)
{
  INT x = fmma(a.i,b.i, c.i,d.i); // a.i*b.i + c.i*d.i;
  INT r = fmma(a.r,b.r, c.r,d.r) - Quad::n*x; // a.r*b.r + c.r*d.r - Quad::n*x;
  INT i = fmma(a.r,b.i, a.i,b.r) + fmma(c.r,d.i, c.i,d.r); // a.r*b.i + a.i*b.r + c.r*d.i + c.i*d.r;
  if (Quad::t)
    i += x;
  return Quad(r,i);
}

// compute a*b-c*d
Quad mms(const Quad& a, const Quad& b, const Quad& c, const Quad& d)
{
  INT x = fmms(a.i,b.i, c.i,d.i); // a.i*b.i - c.i*d.i;
  INT r = fmms(a.r,b.r, c.r,d.r) - Quad::n*x; // a.r*b.r - c.r*d.r - Quad::n*x;
  INT i = fmma(a.r,b.i, a.i,b.r) - fmma(c.r,d.i, c.i,d.r); // a.r*b.i + a.i*b.r - (c.r*d.i + c.i*d.r);
  if (Quad::t)
    i += x;
  return Quad(r,i);
}

void Quad::addprod(const Quad& a, const Quad& b) // this +=a*b
{
  if (a.nm==0 || b.nm==0) return;
  INT c = a.i*b.i;
  r += fmms(a.r,b.r, n,c);
  i += fmma(a.r,b.i, a.i,b.r);
  if (t)
    i += c;
  setnorm();
}

void Quad::addprod(long a, const Quad& b) // this +=a*b
{
  if (a==0 || b.nm==0) return;
  r += a*b.r;
  i += a*b.i;
  setnorm();
}

void Quad::subprod(const Quad& a, const Quad& b) // this -=a*b
{
  if (a.nm==0 || b.nm==0) return;
  INT c = a.i*b.i;
  r -= fmms(a.r,b.r, n,c);
  i -= fmma(a.r,b.i, a.i,b.r);
  if (t)
    i -= c;
  setnorm();
}

void Quad::subprod(long a, const Quad& b) // this -=a*b
{
  if (a==0 || b.nm==0) return;
  r -= a*b.r;
  i -= a*b.i;
  setnorm();
}

Quad qdivi0(const Quad& a, const INT& c) // c>0,    // used when t=0
{
  return (c>0?
          Quad(rounded_division(a.r,c), rounded_division(a.i,c))
          :
          Quad(rounded_division(-a.r,-c), rounded_division(-a.i,-c)));
}

Quad qdivi1(const Quad& a, const INT& c) // used when t=1
{
  INT ansi = (c>0? rounded_division(a.i,c) : rounded_division(-a.i,-c));
  return (c>0?
          Quad(rounded_division(2*a.r+a.i-c*ansi,2*c), ansi)
          :
          Quad(rounded_division(-2*a.r-a.i+c*ansi,-2*c), ansi));
}

// static function (one for the class, not per instance)
void Quad::field(long dd, long max)
{
  // Clear these in case this is not the first field run

  geometry_initialised = 0;
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
  SD.clear();

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

  // append odd primes divisors of d with sign to be 1 mod 4:
  vector<INT> pp = pdivs(odd);
  std::for_each(pp.begin(), pp.end(), [] (INT& p) { if (p%4==3) p=-p;});
  prime_disc_factors.insert(prime_disc_factors.end(), pp.begin(), pp.end());

  INT i0(0), i1(1);
  w = Quad(i0, i1, n);
  zero = Quad(i0,i0, i0);
  one = Quad(i1,i0, i1);
  shifts_by_one = {-one-w, -one, -one+w, -w, zero, w, one-w, one, one+w};

  switch (d) {
  case 1:  pos=&pos13; name='i'; nunits=4; fundunit=w; break;
  case 2:  pos=&pos2;  name='t'; nunits=2; fundunit=-one; break;
  case 3:  pos=&pos13; name='w'; nunits=6; fundunit=w; break;
  default: pos=&pos2;  name='a'; nunits=2; fundunit=-one;
  }

  quadgcd=&quadgcd_default;
  quadbezout=&quadbezout_default;

  quadunits.push_back(one);
  quadunits.push_back(fundunit);
  for(int i=2; i<nunits; i++)
    quadunits.push_back(fundunit*quadunits[i-1]);
  for(int i=0; 2*i<nunits; i++)
    squareunits.push_back(quadunits[2*i]);
  maxnorm=max;
  if(class_number==1)
    initquadprimes();
  Quadprimes::init(max);
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

int Quad::chi(const INT& p)
{
  return (p==2? (d%4==3? (d%8==3? -1: +1): 0):  legendre(disc,p));
}

Quad makepos(const Quad& a)
{
  Quad ans=a;
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
  nm = fmma(r,r, n*i,i);
  if (t) {nm += r*i;};
  // INT N = r*r+t*r*i+n*i*i;
  // if (N!=nm) cerr<<"In setnorm(): Quad "<<(*this)<<" with r = "<<r<<", i = "<<i
  //                <<" and norm "<<N<<" has nm = "<<nm<<endl;
  // assert(N==nm);
  assert (nm>=0);
}

Quad Quad::operator/ (const Quad& b) const
{
  if (!::is_zero(b.i))
    return qdivi(mult_conj(*this,b), b.nm);
  else
    return qdivi(*this,b.r);
}

void Quad::operator/=(const Quad& b)
{
  if (!::is_zero(b.i))
    *this=qdivi(mult_conj(*this,b), b.nm);
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
  return ((i & ideal_class_mod_squares(I)) ? -1 : +1); // bitwise &
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
 Quad c = mult_conj(b,a);
 return (((c.r)%a.nm)==0) && (((c.i)%a.nm)==0);
}

// as above but return quotient b/a when a|b
int div(const Quad& a, const Quad& b, Quad& quo)
{
  if (a.nm==0)
    {
      quo=0;
      return (b.nm==0);
    }
  if (b.nm==0)
    {
      quo=0;
      return 1;
    }
  if (b.nm%a.nm!=0)
    return 0;
  Quad c = mult_conj(b,a);
  INT qr, qi, rr, ri;
  if ( divrem(c.r, a.nm, qr, rr) && divrem(c.i, a.nm, qi, ri) )
    {
      quo = Quad(qr,qi);
      return 1;
    }
  return 0;
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

vector<Quad> invertible_residues(const Quad& a)
{
  vector<Quad> res = residues(a);
  vector<Quad> ires(res.size()); // will be shrunk later
  auto it = std::copy_if(res.begin(), res.end(), ires.begin(),
                         [a] (const Quad& r) {return coprime(a,r);});
  ires.resize(std::distance(ires.begin(), it));
  return ires;
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

void factorp0(long p, INT& a, INT& b, const INT& d)
// finds a,b s.t. a^2+d*b^2=0 (mod p)
{ int found=0;
  for (int ib=1; !found; ib++)
  {
    b = ib;
    found = isqrt(p - d*b*b, a);
  }
}

void factorp1(long p, INT& a, INT& b, const INT& d)
// finds a,b s.t. a^2+a*b+((d+1)/4)*b^2=0 (mod p)
{ int found=0; long fourp = 4*p;
  for (int ib=1; !found; ib++)
  {
    b = ib;
    found = isqrt(fourp -d*b*b,a);
  }
  a=(a-b)/2;
}

vector<Quad> Quad::primes_above(long p, int& sig)
{
  INT d(Quad::d), P(p), a, b;
  Quad pi, piconj;
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
    if(Quad::t) factorp1(p,a,b,d); else factorp0(p,a,b,d);
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
  for (primevar pr; pr.ok()&&pr<maxnorm; pr++)
    { long p=pr;
      list = Quad::primes_above(p, sig);
      if (sig==-1)
        {
          if (p*p<=maxnorm)
            list2.push_back(list[0]);
        }
      else
        list1.insert(list1.end(), list.begin(), list.end());
    }

  // Now list1 contains the degree 1 primes and list2 the degree 2
  // primes, each ordered by norm; we merge these lists so that they
  // are still sorted by norm:

  auto alpha=list1.begin(), beta=list2.begin();
  while ((alpha!=list1.end()) && (beta!=list2.end()))
    if ((alpha->nm)<(beta->nm))
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
  INT na = a.norm();
  if (na<2) return Quad::zero;   // must return something!
  for ( const auto& p : quadprimes)
    {
      if (div(p,a)) return p;
      INT np = p.norm();
      if (np*np>na) return makepos(a);
    }
  cout<<"No prime divisor found for "<<a<<" so assuming prime!\n";
  return makepos(a);
}

vector<Quad> pdivs(const Quad& aa)
{ Quad a=aa; INT norma = a.norm();
  vector<Quad> plist; // will hold prime factors
  if (norma<2) return plist;
  for ( const auto& p : quadprimes)
     {
       if (div(p,a))
	 {
	   plist.push_back(p);
	   while (div(p,a, a));
	   norma = a.norm();
	   if (norma==1) return plist;
	 }
       else
	 {
	   INT normp = p.norm();
	   if (normp*normp>norma)
	     {
	       plist.push_back(makepos(a));
	       return plist;
	     }
	 }
     }
  //In case of p-factors outside range, assume the cofactor is prime:
  if (a.norm()>1)
    plist.push_back(makepos(a));
  return plist;
}

vector<Quad> posdivs(const Quad& a)   // all "positive" divisors (up to units)
{
  vector<Quad> plist=pdivs(a);
  int nu = 1; int nd=nu;
  vector<int> elist;
  for ( const auto& p : plist)
    {
      int e=val(p,a);
      elist.push_back(e);
      nd*=(1+e);
    }
  vector<Quad> dlist(nd);
  dlist[0]=Quad::one;
  nd=nu;
  auto pr=plist.begin();
  auto ei=elist.begin();
  for ( ; pr!=plist.end(); ++pr, ++ei)
   {
     Quad p = *pr;
     int e = *ei;
     for (int j=0; j<e; j++)
       for (int k=0; k<nd; k++)
         dlist[nd*(j+1)+k] = makepos(p*dlist[nd*j+k]);
     nd*=(e+1);
   }
 return dlist;
}

vector<Quad> alldivs(const Quad& a)       // all divisors
{
  vector<Quad> plist=pdivs(a);
  int nu = Quad::nunits; int nd=nu;
  vector<int> elist;
  for ( const auto& p :plist)
   {
     int e=val(p,a);
     elist.push_back(e);
     nd*=(1+e);
   }
  vector<Quad> dlist(nd);
  for(int i=0; i<nu; i++) dlist[i]=quadunits[i];
  nd=nu;
  auto pr=plist.begin();
  auto ei=elist.begin();
  for ( ; pr!=plist.end(); ++pr, ++ei)
   {
     Quad p = *pr;
     int e = *ei;
     for (int j=0; j<e; j++)
       for (int k=0; k<nd; k++)
         dlist[nd*(j+1)+k] = p*dlist[nd*j+k];
     nd*=(e+1);
   }
 return dlist;
}

vector<Quad> sqdivs(const Quad& a) // all divisors whose square divides a, up to +/-
{
  vector<Quad> plist=pdivs(a);
  int nu = Quad::nunits/2; int nd=nu;
  vector<int> elist;
  for ( const auto& p : plist)
   {
     int e=val(p,a)/2;
     elist.push_back(e);
     nd*=(1+e);
   }
  vector<Quad> dlist(nd);
  for(int i=0; i<nu; i++) dlist[i]=quadunits[i];
  nd=nu;
  auto pr=plist.begin();
  auto ei=elist.begin();
  for ( ; pr!=plist.end(); ++pr, ++ei)
   {
     Quad p = *pr;
     int e = *ei;
     for (int j=0; j<e; j++)
       for (int k=0; k<nd; k++)
         dlist[nd*(j+1)+k] = p*dlist[nd*j+k];
     nd*=(e+1);
   }
 return dlist;
}

vector<Quad> sqfreedivs(const Quad& a)       // all square-free divisors
{
  vector<Quad> plist=pdivs(a);
  int nu = 2; int nd=pow(2,plist.size()+1);
  vector<int> elist(plist.size(), 1);
  vector<Quad> dlist(nd);
  for(int i=0; i<nu; i++) dlist[i]=quadunits[i];
  nd=nu;
  auto pr=plist.begin();
  auto ei=elist.begin();
  for ( ; pr!=plist.end(); ++pr, ++ei)
   {
     Quad p = *pr;
     int e = *ei;
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
 else {
   cerr<<"invmod called with "<<a<<" and "<<p<<" -- not coprime!"<<endl;
   return Quad::zero;
 }
}

int invertible(const Quad& a, const Quad& b, Quad& inverse)
{ Quad y; Quad g = quadbezout(a,b,inverse,y);
  return g.is_one();
}

//#define test_reduce

// returns n = round(true_real_part_of(alpha/beta)), so alpha-n*beta
// is reduced mod Z<beta>
INT nearest_INT_to_Quad_quotient ( const Quad& alpha, const Quad& beta)
{
  return rounded_division(racb(alpha,beta), beta.norm());
}

// reduction of gamma modulo Z<alpha,beta>
Quad reduce_mod_zbasis(const Quad& gamma, const Quad& alpha, const Quad& beta)
{
#ifdef test_reduce
  cout << "reduction of "<<gamma<<" mod <"<<alpha<<","<<beta<<">"<< flush;
#endif
  INT d = iacb(beta, alpha); // = (beta*quadconj(alpha)).im();
  assert (d>0);
  INT x = rounded_division(-iacb(gamma, beta), d);
  INT y = rounded_division( iacb(gamma, alpha), d);
  Quad ans = gamma - (x*alpha + y*beta);
#ifdef test_reduce
  cout << " is "<< ans << " (d="<<d<<", x="<<x<<", y="<<y<<")"<<endl;
  cout << " x*alpha =  "<< x*alpha << endl;
  cout << " y*beta =  "<< y*beta << endl;
  cout << " x*alpha+y*beta =  "<< x*alpha+y*beta << endl;
#endif
  INT gn = ans.norm();
  vector<Quad> tests = {ans+alpha, ans-alpha, ans+beta, ans-beta,
                        ans-alpha-beta, ans-alpha+beta, ans+alpha-beta, ans+alpha+beta};
  for ( const auto& t : tests)
    {
      if (gn > t.norm())
        {
#ifdef test_reduce
          cout<<"reduction "<<ans<<" has larger norm than shift "<<t<<", switching"<<endl;
#endif
          ans = t;
          gn = t.norm();
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
          Quad temp = -alpha; alpha = beta; beta = temp;
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

  if (iacb(beta, alpha) < 0)
    beta=-beta;

  // If we wanted to make sure that beta/alpha is unambiguously defined when on the boundary:
  // if (2*racm(beta, alpha) == -alpha.norm())
  //   beta+=alpha;

#ifdef test_reduce
  cout<<"After reduction, alpha="<<alpha<<", beta="<<beta
      <<" with norms "<<alpha.norm()<<", "<<beta.norm()<<endl;
  assert (alpha.norm()<=beta.norm());
  //assert (nearest_INT_to_Quad_quotient(alpha,beta)==0);
  assert (iacb(beta, alpha) > 0);
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
//           Quad temp = -alpha; alpha = beta; beta = temp;
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

//   if (iacb(beta,alpha) < 0)
//     beta=-beta;

//   // If we wanted to make sure that beta/alpha is unambiguously defined when on the boundary:
//   // if (2*tacm(beta, alpha) == -alpha.norm())
//   //   beta+=alpha;

// #ifdef test_reduce
//   cout<<"After reduction by U="<<U<<", alpha="<<alpha<<", beta="<<beta
//       <<" with norms "<<alpha.norm()<<", "<<beta.norm()<<endl;
//   U11 = I2long(U(1,1)); U12 = I2long(U(1,2));
//   U21 = I2long(U(2,1)); U22 = I2long(U(2,2));
//   assert (U11*alpha+U12*beta == alpha0);
//   assert (U21*alpha+U22*beta == beta0);
//   assert (alpha.norm()<=beta.norm());
//   assert (nearest_INT_to_Quad_quotient(alpha,beta)==0);
//   assert (iacb(beta, alpha) > 0);
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
  return {a*g, b*g, g};
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
  if (a.is_zero() || b.is_zero()) return 0;
  if ((a==b) || (a==-b)) return 1;
  if (Quad::nunits !=2) return 0;
  if(a.norm() != b.norm()) return 0;
  return std::any_of(quadunits.begin(), quadunits.end(),
                     [a,b] (const Quad& eps) {return a*eps==b;});
}

int is_ideal_Galois_stable(const Quad& a)
{
  return are_associate(a, a.conj());
}

string ideal_code(const Quad& N) // string code for a (principal) ideal N
{
  vector<INT>H = HNF(N);
  stringstream s;
  s << N.norm() << "." << H[1] << "." << H[2];
  return s.str();
}

vector<int> makechitable(const Quad& lambda, const vector<Quad>& reslist)
{
  if (lambda.is_unit())
    return {1};
  vector<int> chi(reslist.size());
  std::transform(reslist.begin(), reslist.end(), chi.begin(),
                 [lambda, reslist] (const Quad& r) {return squaremod(r,lambda,reslist);});
  return chi;
}

// brute force test whether a is a square of some element of reslist, mod m

int squaremod(const Quad& a, const Quad& m, const vector<Quad>& reslist)
{
  if (div(m,a))
    return 0;
  if (std::any_of(reslist.begin(), reslist.end(), [a,m] (const Quad& r) {return div(m,r*r-a);}))
    return +1;
  else
    return -1;
}

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

vector<int> chardisc(const INT& D)
{
  vector<int> ans(Quad::prime_disc_factors.size());
  std::transform(Quad::prime_disc_factors.begin(), Quad::prime_disc_factors.end(),
                 ans.begin(),
                 [D] (const INT& d) {return div_disc(d,D);});
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
  vector<INT> ans(Quad::all_disc_factors.size());
  auto it = std::copy_if(Quad::all_disc_factors.begin(), Quad::all_disc_factors.end(), ans.begin(),
               [i,j] (const INT& Di) {return ((chardisc(Di)[i]==0) && (chardisc(Di)[j]==0));});
  ans.resize(std::distance(ans.begin(), it));
  // cout<<"All discs "<<Quad::all_disc_factors<<" mod "<<D<<": "<<ans<<endl;
  return ans;
}

void Quad::setup_geometry(string subdir, int debug)
{
  if (!geometry_initialised)
    {
      SD.read_or_create(subdir, debug);
      geometry_initialised = 1;
    }
}

// return name of newforms directory for this field; if
// create_if_necessary, creates the directory if it does not yet exist
string newforms_directory(int create_if_necessary)
{
  stringstream s;
  s << getenv_with_default("NF_DIR", "newforms") << "/" << field_label();
  string dirname = s.str();
  if (create_if_necessary)
    {
      int res = std::system(("mkdir -p "+dirname).c_str());
      if (res)
        {
          cerr << "mkdir -p "<<dirname<<" failed, writing newforms files in current directory"<<endl;
          dirname = ".";
        }
    }
  return dirname;
}

// return name of Newspace directory for this field; if
// create_if_necessary, creates the directory if it does not yet exist
string newspaces_directory(int create_if_necessary)
{
  stringstream s;
  s << getenv_with_default("NS_DIR", "newspaces") << "/" << field_label();
  string dirname = s.str();
  if (create_if_necessary)
    {
      int res = std::system(("mkdir -p "+dirname).c_str());
      if (res)
        {
          cerr << "mkdir -p "<<dirname<<" failed, writing newforms files in current directory"<<endl;
          dirname = ".";
        }
    }
  return dirname;
}

string eigfile(const Quad& N, long p)    //returns filename for eigs at level N, characteristic p
{
  stringstream s;
  s << newforms_directory() << "/" << ideal_code(N);
  if (p)
    s << "_mod_"<<p;
  return s.str();
}

string eigfile(Qideal& N, long p)    //returns filename for eigs at level N, characteristic p
{
  stringstream s;
  s << newforms_directory() << "/";
  if (Quad::class_number==1) // for backwards compatibility of data file names
    s << ideal_code(N.gen());
  else
    s << ideal_label(N);
  if (p)
    s << "_mod_"<<p;
  return s.str();
}

