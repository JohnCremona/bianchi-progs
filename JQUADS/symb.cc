// FILE symb.cc   ---  Implementation for cusps, M-symbols, modular symbols
// new version based on matrix representation of cusps

#include "symb.h"

#undef DEBUG_CUSP_EQ

n_cusp* geometry::basic_cusps;
long geometry::num_etypes;
n_modsym* geometry::basic_edges;
long geometry::aitch2_nonprinc_det_value;
long geometry::max_length_of_face_rel;
long geometry::max_relmat_shape_ratio;

long pseudo_euclid::bubblenum;
n_mat* pseudo_euclid::mats;
n_mat* pseudo_euclid::imats;

abstract_mat::abstract_mat(const Quad& ia, const Quad& ib, const Quad& ic, const Quad& id) : aa(ia), bb(ib), cc(ic), dd(id)
{
  Quad tmp = det();
  if (tmp.pos_assoc()==1) iclass=0;
  else
    if (tmp==geometry::aitch2_nonprinc_det_val()) iclass=1;
    else
      { cerr<<"Error:det(abstract_mat) = "<<tmp<<" needs explicit typeflag!\n";
	long d=0,e=0;
	d/=e;             // for debugging - generates SIGFPE exception
	exit(1);
      }
}

long class_group_product(long c1, long c2, long& cf)
{
  long ans;

  if ((c1==-1)||(c2==-1)) { ans = -1; cf = 1; }
  else

    // for now, never called when class group is trivial
    // if (Field::is_princ())  { ans=0; cf=1; }
    // else

  // for now, wma Field::classnum()==2

  if (c1!=c2) { ans=1; cf=1; }           // princ * non-princ = non-princ
  else
    {
      ans=0;                             // product is principal
      if (c1==0)
	cf=1; 
      else                               // non-princ * non-princ = princ ... 
	cf = geometry::aitch2_nonprinc_det_val();  // ... with a common factor
    }
  return ans;
}

void abstract_mat::operator/=(const long&cf)
{
  if (cf!=1) { aa/=cf; bb/=cf; cc/=cf; dd/=cf; }
}

abstract_mat abstract_mat::conj() const
{
  abstract_mat ans(aa.conj(), bb.conj(), cc.conj(), dd.conj(), iclass);
  return ans;
}

abstract_mat abstract_mat::operator*(const abstract_mat&m) const    // mx multn
{
  abstract_mat ans;
  ans.aa = aa*m.aa + bb*m.cc;
  ans.bb = aa*m.bb + bb*m.dd;            //   ( a  b )(m.a m.b)
  ans.cc = cc*m.aa + dd*m.cc;            //   ( c  d )(m.c m.d)
  ans.dd = cc*m.bb + dd*m.dd;

  if ((cl()==-1)||(m.cl()==-1))
    { ans.iclass=-1; }
  else
    if (Field::is_princ())
      { ans.iclass=0; }
    else
      {
	long cf;
	ans.iclass = class_group_product(iclass, m.iclass, cf);
	ans/=cf;
      }
  return ans;
}

n_cusp n_mat::operator()(const n_cusp& q) const
{
  n_cusp ans=(*this)*q;              // uses multn of underlying abstract_mats

// Comment out next bit - although our cusps won't have the nice standard
// representation after being hit by non-principal matrices coming from Hecke
// operators, we don't want to call quadbezout unless we really have to.
/*
  if (ans.cl()==-1)
    {
      ans = ans.infmat();
    }
*/

  // finally, could try to simplify right-hand column

  return ans;
}

n_cusp n_cusp::infmat() const
{
  abstract_mat ans;
  abstract_mat temp;
  Quad x,y,g;
  int found=0;
  long i=0, imax=Field::classnum();
  do
    {
      temp = ((n_mat)(abstract_mat)(geometry::bc(i))).inv_mod_scalars()
	* (abstract_mat)(*this);
      g = quadbezout(temp.a(), temp.c(), x, y);
      if (g.div(temp.a()) && g.div(temp.c()))
	found=1;
      else i++;
    }
  while ((!found) && (i<imax));
  if (found)
    {
      ans = geometry::bc(i) * abstract_mat(temp.a()/g, -y, temp.c()/g, x, 0);
    }
  else
    {
      cerr << "Error: infmat can't find orbit." << endl;
      exit(1);
    }
  
  if ( (n_cusp)ans != *this )
    {
      cerr << "Error: infmat of infty not equal to given cusp!"<<endl;
      exit(1);
    }
  return (n_cusp)ans;
}

n_modsym n_modsym::conj() const
{
  n_modsym ans(a.conj(), b.conj());
  return ans;
}

n_mat n_mat::specialmat(const Qideal&p)
// INPUT: an ideal p whose square is principal
// OUTPUT: a matrix giving an isomorphism   R (+) R  --->  p (+) p
{
  Quad x,y;
  Quad g1 = p.cee() * p.ay();
  Quad g2 = p.cee() * Quad(p.bee(),1);      // i.e.  p = <g1,g2>
  Quad h1 = g1*g1;
  Quad h2 = g2*g2;                          // FACT: p^2=<h1,h2>
  Quad beta = quadbezout(h1, h2, x, y);
  if (ndiv(beta,h1) || ndiv(beta,h2))
    {
      cerr << "Error in specialmat." << endl;
      exit(1);
    }
  n_mat ans(g1, -g2*y, g2, g1*x, -1);
  return ans;
}

void pseudo_euclid::init()
{
  // initialise bubblenum, the no. of spheres needed modulo translations
  //
  if (Field::is_euc()) bubblenum=1;
  else
  switch (Field::d)
    {
    case 5:
      bubblenum = 2;
      break;
    default:
      cerr<< "Error: pseudo-Euclid not implemented for field "<<Field::d<<endl;
      exit(1);
    }

  // allocate memory

  mats = new n_mat[bubblenum];
  imats = new n_mat[bubblenum];

  // initialise the list of inversion matrices

  if (Field::is_euc())
    {
      mats[0] = n_mat(0,-1, 1, 0, 0);
    }
  else
  switch (Field::d)
    {
    case 5:
      mats[0] = n_mat(0,-1, 1, 0, 0);
      mats[1] = n_mat(Quad(1,1), Quad(1,-1), 2, Quad(-1,-1), 1);
      break;
    default:
      cerr<< "Error: pseudo-Euclid not implemented for field "<<Field::d<<endl;
      exit(1);
    }

  // compute their inverses (modulo scalars)
  for (long i=0; i<bubblenum; i++)
    {
      imats[i] = mats[i].inv_mod_scalars();
    }
}

n_cusp pseudo_euclid::co_centre(long j)
// return the co-centre of the j'th inversion hemisphere, i.e. image of
//   infinity under the j'th inversion matrix
{
  if ((j<0)||(j>=bubblenum))
    { cerr << "Error: index out of range in pseudo_euclid::co_centre"<<endl;
      exit(1);
    }
  return (abstract_mat)(mats[j]);
}

n_cusp pseudo_euclid::centre(long j)
// return the centre of the j'th inversion hemisphere, i.e. the point mapped
//   to infinity by the j'th inversion matrix, i.e. the image of infinity
//   under the inverse of the j'th inversion matrix
{
  if ((j<0)||(j>=bubblenum))
    { cerr << "Error: index out of range in pseudo_euclid::centre"<<endl;
      exit(1);
    }
  return (abstract_mat)(imats[j]);
}

/*
OLD (1996) VERSION
void pseudo_euclid::translation_step()
// find q,j such that z-q is best covered by hemisphere j
// replace z by z-q and m by m*T(q)
{
  if (z.c() == 0)
    {
      cerr << "Error: translation_step called for cusp at infinity!"<<endl;
      exit(1);
    }

  Quad q;
  if (Field::is_euc())
    {
      q = z.a() / z.c();                // rounded division
      j=0;
    }
  else
  switch (Field::d)
    {
    case 5:
      {
	long den;
	Quad conj;
	if (z.c().small_integer_multiple(den,conj))
	  {
	    cerr << "Error: overflow in translation_step."<<endl;
	    long e,f=0; e/=f;
	  }
	Quad num=z.a()*conj;
	long x=num.real();
	long y=num.imag();

	long qx=roundover(x,den);
	long newx=x-qx*den;

	long qy=(long)floor( ((double)(2*den -abs(newx) +5*y)) / (5*den) );
	long newy=y-qy*den;

	// check that this has worked
	if ((2*abs(newx)>den)||
	    (abs(newx)-5*newy <=-3*den)||
	    (abs(newx)-5*newy  > 2*den))
	  {
	    cerr << "Error: translation failed to find inversion region"<<endl;
	    long e,f=0; e/=f;
	  }

	if (newx<0)
	  if (newx-5*newy +2*den < 0) { newx+=den; qx-=1; }

	if (newx+5*newy <= 2*den ) j=0; else j=1;      // which sphere covers?
	
	q=Quad(qx,qy);
      }
      break;

    default:
      cerr<<"Error: translation_step not implemented for field "<<Field::d<<endl;
      exit(1);
    }

  z = n_mat::T(-q)(z);
  m = m * n_mat::T(q);
}
*/

// BIGINT VERSION (1997)
void pseudo_euclid::translation_step()
// find q,j such that z-q is best covered by hemisphere j
// replace z by z-q and m by m*T(q)
{
  if (z.c() == 0)
    {
      cerr << "Error: translation_step called for cusp at infinity!"<<endl;
      exit(1);
    }

  Quad q;
  if (Field::is_euc())
    {
      q = z.a() / z.c();                // rounded division
      j=0;
    }
  else
  switch (Field::d)
    {
    case 5:
      {
	bigint a1=BIGINT(z.a().real()), a2=BIGINT(z.a().imag());
	bigint c1=BIGINT(z.c().real()), c2=BIGINT(z.c().imag());
	bigint den = c1*c1 + 5*c2*c2;
	bigint x = a1*c1 + 5*a2*c2;
	bigint y = a2*c1 - a1*c2;

	bigint qx;
	nearest(qx, x, den);      // set qx to nearest integer to zx/den
	bigint newx = x - qx*den;

	bigint qy = floor( (2*den - abs(newx) + 5*y), 5*den );
	bigint newy = y - qy*den;

	// check that this has worked
	if ((2*abs(newx) > den) ||
	    (abs(newx)-5*newy <= -3*den)||
	    (abs(newx)-5*newy  > 2*den))
	  {
	    cerr << "Error: translation failed to find inversion region"<<endl;
	    long e,f=0; e/=f;
	  }

	if (newx<0)
	  if (newx-5*newy +2*den < 0) { newx+=den; qx-=1; }

	if (newx+5*newy <= 2*den ) j=0; else j=1;      // which sphere covers?

	long ansx, ansy;
	if ( !is_long(qx) || !is_long(qy))
	  {
	    cerr << "Error: bigint too large for long in pseudo_euclid::translation_step\n";
	    long e,f=0; e/=f;
	  }
	longasI(ansx,qx);  longasI(ansy,qy);
	q=Quad(ansx,ansy);
      }
      break;

    default:
      cerr<<"Error: translation_step not implemented for field "<<Field::d<<endl;
      exit(1);
    }

  z = n_mat::T(-q)(z);
  m = m * n_mat::T(q);
}

/*
OBSOLETE (1997) VERSION USING DOUBLES --- hopeless truncation errors!
long round_double_to_long(const double &x)
// round to nearest integer (halves away from zero)
// assumes built-in conversion truncates towards zero, which is
// implementation dependent, not guaranteed by the language 
{
  long ans;
  if (x>=0)
    { ans = (long)(x+0.5); }
  else
    { ans = (long)(x-0.5); }
  return ans;
}

void pseudo_euclid::translation_step()
// find q,j such that z-q is best covered by hemisphere j
// replace z by z-q and m by m*T(q)
{
  if (z.c() == 0)
    {
      cerr << "Error: translation_step called for cusp at infinity!"<<endl;
      exit(1);
    }

  Quad q;
  if (Field::is_euc())
    {
      q = z.a() / z.c();                // rounded division
      j=0;
    }
  else
  switch (Field::d)
    {
    case 5:
      {
	double a1=z.a().real(), a2=z.a().imag();
	double c1=z.c().real(), c2=z.c().imag();
	double den = c1*c1 + 5*c2*c2;
	double zx = (a1*c1 + 5*a2*c2)/den;
	double zy = (a2*c1 - a1*c2)/den;

	long qx= round_double_to_long(zx);
	double newzx = zx-qx;

	double temp = abs(newzx) - 5*zy;
	long qy = (long)floor( (2-temp)/5 );
	double newzy = zy - qy;

	// check that this has worked
	if ((2*abs(newzx)>1)||
	    (abs(newzx)-5*newzy <=-3)||
	    (abs(newzx)-5*newzy  > 2))
	  {
	    cerr << "Error: translation failed to find inversion region"<<endl;
	    long e,f=0; e/=f;
	  }

	if (newzx<0)
	  if (newzx-5*newzy +2 < 0) { newzx+=1; qx-=1; }

	if (newzx+5*newzy <= 2 ) j=0; else j=1;      // which sphere covers?
	
	q=Quad(qx,qy);
      }
      break;

    default:
      cerr<<"Error: translation_step not implemented for field "<<Field::d<<endl;
      exit(1);
    }

  z = n_mat::T(-q)(z);
  m = m * n_mat::T(q);
}
*/

void pseudo_euclid::reexpress()
//
// So far, each step yields an expression
//      M_i { beta(j), infty }
// which needs to be re-expressed in the form
//     (+/-)  (matrix in SL2(O)) acting on (standard edge representative)
//      eps      M_i * A              A^{-1} { beta(j), infty }
// the latter term being the t'th standard edge for some t.
//
// One way would be to pre-initialise an array measuring h x bubblenum,
// with the (c,j) entry being the re-expression triple (A, eps, t) for
// a matrix of class c acting on the symbol { beta(j), infty }.
//
// For now, here is an ad-hoc method using if-statements.
{
  if (Field::is_euc())
    {
      eps=+1;
      t=0;
      ma=m;
    }
  else
  switch (Field::d)
    {
    case 5:
      {
	if (m.cl()==0)
	  {
	    eps=+1;
	    t=j;
	    ma=m;
	  }
	else
	  {
	    if (j==0)
	      {
		eps=+1;
		t=2;
		ma=m* n_mat( Quad(-1,1), 2, 2, Quad(-1,-1), 1);        // *B^3
	      }
	    else
	      {
		eps=-1;
		t=1;
		ma=m*n_mat( Quad(1,1), Quad(1,-1), 2, Quad(-1,-1), 1); // *TB^3
	      }
	  }
      }
      break;
    default:
      {
	cerr << "Error: pseudo_euclid cannot re-express for field "<<Field::d<<endl;
	exit(1);
      }
    }
}

int pseudo_euclid::not_infty()
// if not there, takes one step towards infty
{
  int ans = (z.c() != 0);
  if (ans)
    {
      translation_step();
      z = mats[j](z);
      m = m * imats[j];
      reexpress();
    }
  return ans;
}

void geometry::init()
{
  init_params();
  init_cusp_info();
  init_edge_info();
  pseudo_euclid::init();
}

void geometry::init_params()
{
  if (Field::is_princ())
    {
      aitch2_nonprinc_det_value=0;   // should never be used!! 
    }

  if (Field::is_euc())                 // Euclidean fields only
    {
      max_relmat_shape_ratio=2;
    }

  switch (Field::d)
    {
    case 1: case 3:
      max_length_of_face_rel=3;
      break;
    case 2: case 7:
      max_length_of_face_rel=4;
      break;
    case 11:
      max_length_of_face_rel=6;
      break;
    case 19:
      max_length_of_face_rel=4;
      max_relmat_shape_ratio=3;
      break;
    case 5:
      aitch2_nonprinc_det_value=2;
      max_length_of_face_rel=4;
      max_relmat_shape_ratio=2;
      break;
    default:
      cerr<< "Error: Geometry not implemented for field "<<Field::d<<" !!!\n";
      exit(1);
    }
}

void geometry::init_cusp_info()
{
  basic_cusps = new n_cusp[ Field::classnum() ];

  if (Field::is_princ())
    {
      basic_cusps[0] = n_cusp::infty();      // (would be default anyway!)
    }
  else
  switch (Field::d)
    {
    case 5:
      basic_cusps[0] = n_cusp::infty();      // (would be default anyway!)
      basic_cusps[1] = n_cusp( Quad(1,1), Quad(1,-1), 2, Quad(-1,-1), 1);
      break;
    default:
      cerr<< "Error: cusp info not implemented for field "<<Field::d<<" !!!\n";
      exit(1);
    }
}

void geometry::alloc_edge_info(long n)
{
  num_etypes = n;
  basic_edges = new n_modsym[n];
}

void geometry::init_edge_info()
{
  if (Field::is_euc())                 // Euclidean fields only
    {
      alloc_edge_info(1);

      basic_edges[0] = n_modsym(n_cusp::zero(),                      // 0
				n_cusp::infty() );                   // infty

//      basic_edges[0] = n_modsym(n_cusp(0,-1,1,0,0),               // 0
//				n_cusp(1,0,0,1,0) );                // infty
    }
  else
  switch (Field::d)
    {
    case 19:
      alloc_edge_info(2);
      basic_edges[0] = n_modsym(n_cusp( 0,-1, 1, 0, 0),               // 0
				n_cusp( 1, 0, 0, 1, 0) );             // infty
      basic_edges[1] = n_modsym(n_cusp(Quad(0,1),2,2,Quad(1,-1), 0),  // w/2
				n_cusp( 1, 0, 0, 1, 0) );             // infty
      break;
    case 5:
      alloc_edge_info(3);
      basic_edges[0] = n_modsym(n_cusp(0,-1,1,0,0),                   // 0
				n_cusp(1,0,0,1,0) );                  // infty
      basic_edges[1] = n_modsym(n_cusp(Quad(1,1),2,2,Quad(1,-1) ,1),  //(1+w)/2
				n_cusp(1,0,0,1, 0) );                 // infty
      basic_edges[2] = n_modsym(n_cusp(2,Quad(-1,-1),Quad(1,-1),-2,1),//(1+w)/3
				n_cusp(Quad(1,1),2,2,Quad(1,-1), 1)); //(1+w)/2
      break;
    default:
      cerr<< "Error: edge info not implemented for field "<<Field::d<<" !!!\n";
      exit(1);
    }
}

void geometry::display(ostream& s, int verbose)
{
  s << "Geometry for field d=" <<  Field::d << endl;
  s << "aitch2_nonprinc_det_value: "<<aitch2_nonprinc_det_val();
  if ((Field::is_princ())&&(aitch2_nonprinc_det_val()!=0))
    s << " --- unexpected !";
  s << endl;
  s << "Number of types of edges: "<< nt() << endl;
  s << "The basic edges:"<<endl;
  for (long i=0; i< nt(); i++) s << i << "\t" << be(i) << endl;
  s << "Length of longest face rel: " << max_length_of_rel() << endl << endl;
  s << "Max num of rels as multiple of ngens: " << relmat_shape() << endl;
}

//Members of former class symblist:
/*
symblist::symblist(long n) : maxnum(n),num(0) {list=new cdsymb[n];}

void symblist::add(const cdsymb& s, long start)
{
  if (index(s,start)==-1) 
    {
      if (num<maxnum) list[num++]=s; 
      else
	{
	  cerr << "Error in symblist::add: ";
	  cerr << "attempt to add too many symbols to list!\n";
	}
    }
}

long symblist::index(const cdsymb& s, long start) const
{
 long i,ans;
 for (i=start,ans=-1; ((i<num)&&(ans==-1)); i++)
   if (equal(list[i],s)) ans=i;
 return ans;
}

cdsymb symblist::item(long n) const
{
 if ((n>num)||(n<0)) 
 {cerr<<"Error in symblist::item: index out of range!\n";
  return cdsymb(0,1);
 }
 else return list[n];
}
*/

//Member functions for class symbdata:
symbdata::symbdata(const Qideal &n) : moddata(n), quotient_ring(n,phi)
{
  nsymb = psi;
  nsymb1=2*normod-phi;
  nsymb2=nsymb-nsymb1;
// "nsymb"xx moved here from moddata:: by JSB

//cout<<normod<<" residues, "<<phi<<" invertible and "<<normod-phi<<" non.\n";
//cout<<nsymb<<" symbols, of which "<<nsymb2<<" are special"<<endl;

specials = new cdsymb[nsymb2];
specialsnum=0;

// cout << "In constructor symbdata::symbdata.\n";
// cout << "nsymb2 = " << nsymb2 << "\n";
  if (nsymb2>0)
    {
      Quad dummy1, dummy2;
      long nonl=noninvlist.size();
      for (long ic=0; (ic<nonl); ic++)
      //nicer:  for (long ic=0; (ic<nonl)&&(specialsnum<nsymb2); ic++)

	{
	  Quad c = resnum(noninvlist[ic]);
	  Qideal qidc = Qideal(c);
	  //cout<<"Looking for specials with c = " << c << endl;
	  for (long id=0; id<nonl; id++)
	  //nicer:  for (long id=0; (id<nonl)&&(specialsnum<nsymb2); id++)

	    {
	      Quad d = resnum(noninvlist[id]);		  
	      Qideal fa = qidc + d;
	      if (comax(fa, n, dummy1, dummy2))
		{
		  cdsymb s = cdsymb(c,d);
		  //cout<<"s = "<<s<<": ";
		  //long oldnum=specialsnum;

		  // Is it there already? ...
		  long i,ind;
		  for (i=0,ind=-1; ((i<specialsnum)&&(ind==-1)); i++)
		    if (equal(specials[i],s)) ind=i;
		  if (ind==-1)         // ... if not, add it!
		    {
		      if (specialsnum<nsymb2) specials[specialsnum++]=s; 
		      else
			{
			  cerr << "Error in constructor symbdata --- attempt";
			  cerr << " to add too many specials to list!\n";
			}
		    }

		  //if(oldnum<specialsnum)
		  //  {cout<<"new one\n";}
		  //else {cout<<"old one\n";}
		}
	    }     // end of d loop
	}      // end of c loop
      if (specialsnum<nsymb2)
	{ cerr << "Problem: makesymbols found only " << specialsnum << " symbols ";
	  cerr << "out of " << nsymb2 << "\n";
	}
    }
}

n_modsym symbdata::modsymoftys(const long&jt) const
{
  n_mat mm=symbtomat(symbnum(softys(jt)));
  return mm(geometry::be(toftys(jt)));
}

int symbdata::equal(const cdsymb&s, const cdsymb&t) const
{ return modulus.divides( s.c * t.d - t.c * s.d );}

n_mat symbdata::symbtomat(const cdsymb&s) const
{
  n_mat ans;
  Quad a,b, c=s.cee(), d=s.dee();
  Quad g = quadbezout(c, d, b, a);               // now g=bc+ad so that ...
  if ( g.div(c) && g.div(d) )
    { ans = n_mat ( a, -b, c/g, d/g, 0 ); }      // 1 = ad/g + bc/g = det(ans)
  else
    {
      Qideal help=Qideal(d);
      for (long i=0; i<nfactn.num_primes(); i++)
	{
	  Qideal p=nfactn.prime(i);
	  while (p.divides(help)) { help/=p;}
	}
      Quad cc12, cc34;
      if (!comax(help, modulus, cc12, cc34))
	{ cerr << "Error in symbtomat -- help and n not comax!"<< endl; }
      c=c*cc12 + cc34;
      g = quadbezout(c, d, b, a);
      if (g!=1)
	{ cerr << "Error in symbtomat: g="<<g<<", should be 1."<< endl;}
      ans = n_mat ( a, -b, c, d, 0 );
    }

#ifdef testbezout
  cdsymb test (ans.c(), ans.d());
  if (!equal(test,s))
    {
      cerr << "Error in symbtomat - returning mx which is not a lift!"<<endl;
    }
#endif

  return ans;
}

n_mat symbdata::special_adjustmat(const n_mat&m) const
// assumes that the first column of m is a valid M-symbol at this level
{
  n_mat ans;
  cdsymb s = symbnum(numsymb2(-m.c(),m.a()));
  n_mat u = symbtomat(s);
  ans = u*m;

  if (!modulus.divides(ans.c()))
    {
      cerr << "Error in special_adjustmat." << endl;
      exit(1);
    }

  // during debugging only - use the old-fashioned adjustmat
  if (ans != adjustmat(m))
    {
      cerr << "Warning: adjustmat and special_adjustmat disagree!" << endl;
    }

  return ans;
}

n_mat symbdata::adjustmat(const n_mat&m) const
{
  n_mat ans;
  int notfound=1;
  for (long i=0; (i<nsymb)&&(notfound); i++)
    {
      cdsymb s = symbnum(i);
      if (modulus.divides( s.cee()*m.a() + s.dee()*m.c() ))
	{
	  notfound=0;
	  n_mat u = symbtomat(s);
	  ans = u*m;
	}
    }
  if (notfound || (!modulus.divides(ans.c())))
    {
      cerr << "Error in adjustmat." << endl;
      exit(1);
    }
  return ans;
}
 
long symbdata::numsymb2(const Quad& cc, const Quad& dd) const
{
  Quad c=cc;
  Quad d=dd;

// for debugging
  if ((c==Quad(10080857,-2172602))&&(d==Quad(988276,-5998394)))
    {
      cerr << "Hup - break here!" <<endl;
    }

 long kd = code(d);
  if (kd>0)                // d invertible, with inverse res[kd]
    return numres(c*resnum(kd));   // (c:d) = (c*res[kd]:1)
  else
  {
    // make d smaller!!
  //  d = resnum(noninvlist[-kd]);

    long kc = code(c);
    if (kc>0)              // (c:d) = (1:res[kc]*d) if c invertible
      return   normod-code(resnum(kc)*d);
    else
      {
	// make c smaller!!
   //	c = resnum(noninvlist[-kc]);
	
//cout<<"\nkc="<<kc<<" kd="<<kd;
	cdsymb s(c,d);
	// find its index in list of specials
	long i,ind;
	for (i=0,ind=-1; ((i<specialsnum)&&(ind==-1)); i++)
	  if (equal(specials[i],s)) ind=i;
	if (ind==-1)
	  { cerr<<"Error in numsymb2 --- special not found!\n";
	    return 0;
	  }
	else return nsymb1+ind;
      }
  }
}

cdsymb symbdata::symbnum(long i) const
{ if (i<normod) return cdsymb(resnum(i),1);
  else
    {
      if (i<nsymb1) return cdsymb(1,resnum(noninvlist[i-normod]));
      else      // return specials.item(i-nsymb1);
	{
	  long n=i-nsymb1;
	  if ((n>specialsnum)||(n<0)) 
	    { cerr<<"Error in symbdata::symbnum --- index out of range!\n";
	      return cdsymb(0,1);
	    }
	  else return specials[n];
	}
    }
}

void symbdata::display(ostream& s, int verbose) const
{ moddata::display(s);
  quotient_ring::display(s, verbose);
//s<<normod<<" residues, "<<phi<<" invertible and "<<normod-phi<<" non.\n";
  s << "Number of symbols = " << nsymb << endl;
  if (nsymb2>0)
    {
      s << "Number of special symbols = " << nsymb2 << endl;
      s << "The specials:"<< endl;
      for(long i=0; i<specialsnum; i++) {s <<i<<":\t"<<specials[i]<<"\n";}
    }
  else { s << "No special symbols."<< endl;}
  if (verbose)
    {
      s << "The full list of symbols:" << endl;
      for (long i=0; i<nsymb; i++) { s << i<<":\t"<<symbnum(i)<<"\n";}
      s << endl;
    }
}

int symbdata::check(int verbose) const
{ int moddataok = set_of_residues::check(verbose);
  int ok=1;
  long i,j;
  cdsymb s;
  for (i=0; i<nsymb; i++)
    {
      //  cout<<i<<": "<<flush;
      s = symbnum(i);  
      //  cout<<s<<": "<<flush;
      j = numsymb(s); 
      //  cout<<j<<endl;
      ok&=(i==j);
      if (i!=j) cout << i << "-->" << s << "-->" << j << endl;
    }
  if(verbose)
    {
      if (ok) cout << "symbols check OK!"<<endl;
      else cout << "symbols check found errors!"<<endl;
    }
  return ok&&moddataok;
}

// Operator to compute the normaliser character
// JSB: Don't expect it to depend on the non-principal pivot p supplied
matop symbdata::Chi_p(const Qideal&pp) const
{
  Qideal p = pp;            // UGLY: make non-const copy for isprincipal()
  if (p.isprincipal())
    { cerr << "Error: Chi_p called when p principal"<<endl;
      exit(1);
    }
  if (p.divides(modulus))
    { cerr << "Error: Chi_p called when p divides level"<<endl;
      exit(1);
    }
  matop ans(1);
  ans[0] = special_adjustmat(n_mat::specialmat(p));
  return ans;
}

// Hecke operator for good principal prime
matop symbdata::T_p(const Qideal&pp) const
{
  Qideal p = pp;            // UGLY: make non-const copy for isprincipal()
  long normp=p.norm();
  long n=1+normp;
  matop ans(n);             // allocates memory for n matrices

  if (!p.isprincipal())
    { cerr << "Error: T_p called when p not principal"<<endl;
      exit(1);
    }

  Quad beta=p.gee0();     // principal generator of p
  if (p.divides(modulus))
    { cerr << "Error: T_p called when p divides level"<<endl;
      exit(1);
    }

  set_of_residues r(p);

  ans[--n] = n_mat(beta, 0, 0, 1, -1);

  for (long i=0; i<normp; i++)
    ans[--n] = n_mat(1, r.resnum(i), 0, beta, -1);

  if (n != 0)
    { cerr << "Error: miscount in T_p."<< endl;
      exit(1);
    }

  return ans;
}

// Added 20-07-1998
// Hecke operator for bad principal prime!
matop symbdata::U_p(const Qideal&pp) const
{
  Qideal p = pp;            // UGLY: make non-const copy for isprincipal()
  long normp=p.norm();
  long n=normp;
  matop ans(n);             // allocates memory for n matrices

  if (!p.isprincipal())
    { cerr << "Error: U_p called when p not principal"<<endl;
      exit(1);
    }

  Quad beta=p.gee0();     // principal generator of p
  if (!p.divides(modulus))
    { cerr << "Error: U_p for bad p called when p is not bad!"<<endl;
      exit(1);
    }

  set_of_residues r(p);

  for (long i=0; i<normp; i++)
    ans[--n] = n_mat(1, r.resnum(i), 0, beta, -1);

  if (n != 0)
    { cerr << "Error: miscount in U_p."<< endl;
      exit(1);
    }

  return ans;
}

// Hecke operator for square of non-principal prime
matop symbdata::T_npnp(const Qideal&p) const
{
  long normp=p.norm();
  long normpp=normp*normp;
  long n=1+normp+normpp;
  matop ans(n);             // allocates memory for n matrices

  Qideal psq=fill(p*p);

  if (!psq.isprincipal())
    { cerr << "Error: T_npnp called when p^2 not principal"<<endl;
      exit(1);
    }

  Quad beta=psq.gee0();     // principal generator of p^2
  Quad x,v;

  if (!comax(psq,modulus,x,v))
    { cerr << "Error: T_npnp called when p divides level"<<endl;
      exit(1);
    }

  set_of_residues r(psq);

  for (long i=0; i<normpp; i++)
    ans[--n] = n_mat(1, r.resnum(i), 0, beta, -1);

  for (long i=0; i<normpp; i++)
    {
      if (p.divides(r.resnum(i)))
	ans[--n] = n_mat(beta, 0, v*r.resnum(i), 1, -1);
    }

  if (n != 1)
    { cerr << "Error: miscount in T_npnp."<< endl;
      exit(1);
    }

  // finally, the specialmat

  ans[0] = special_adjustmat(n_mat::specialmat(p));
  return ans;
}

// Hecke operator for product of distinct non-principal primes
matop symbdata::T_npnq(const Qideal&p, const Qideal&q) const
{
  long normp=p.norm();
  long normq=q.norm();
  long n=(1+normp)*(1+normq);
  matop ans(n);             // allocates memory for n matrices

  Quad x,y;
  Qideal pq=fill(p.princprod(q, x, y));
  // recall:  x,y \in q satisfy  (ac)x + (bc+cw)y = principal g'tor of pq

  if (!pq.isprincipal())
    { cerr << "Error: T_npnq called when pq not principal"<<endl;
      exit(1);
    }

  Quad beta=pq.gee0();     // principal generator of pq
  Quad u,v;

  if (!comax(pq,modulus,u,v))
    { cerr << "Error: T_npnq called when p or q divides level"<<endl;
      exit(1);
    }

  // the specialmats
  ans[0] = special_adjustmat( n_mat(p.ay()*p.cee(), -y, p.cee()*Quad(p.bee(), 1), x, -1) );
  ans[1] = special_adjustmat( n_mat(-y, p.ay()*p.cee(), x, p.cee()*Quad(p.bee(), 1), -1) );
  ans[2] = n_mat(beta, 0, 0, 1, -1);

  set_of_residues r(pq);

  for (long i=0; i<normp*normq; i++)
    ans[--n] = n_mat(1, r.resnum(i), 0, beta, -1);

  for (long i=1; i<normp*normq; i++)      // i=0 would give the residue 0
    {
      if (p.divides(r.resnum(i)))
	ans[--n] = n_mat(beta, 0, v*r.resnum(i), 1, -1);
    }

  for (long i=1; i<normp*normq; i++)      // i=0 would give the residue 0
    {
      if (q.divides(r.resnum(i)))
	ans[--n] = n_mat(beta, 0, v*r.resnum(i), 1, -1);
    }

  if (n != 3)
    { cerr << "Error: miscount in T_npnq."<< endl;
      exit(1);
    }
  return ans;
}

/*
matop::matop(const Quad& p, const Quad& n)
{
 if (p==n)
   {
     length=1;  mats=new mat22[1]; mats[0]=mat22(0,-1,n,0);
   }
 else
 if (div(p,n))   // W involution, 1 term
   {
      length=1;
      Quad u,v,a,b;
      for (u=1, v=n; div(p,v); v/=p, u*=p) ;
      quadbezout(u,v,a,b);
      mats = new mat22[1];
      mats[0]=mat22(u*a,-b,n,u);
   }
else                 // Hecke operator, p+1 terms
  {
    Quadlist resmodp = residues(p);
    int normp = resmodp.length;
    length=normp+1;
    mats = new mat22[length];
    for (Quadvar j(resmodp); j.ok(); j++) 
      mats[j.index] = mat22(1,(Quad)j,0,p);
    mats[normp] = mat22(p,0,0,1);
  }
}
*/

// WAS IN SEPARATE FILE cusp.cc

int cusplist::eq(const n_cusp&m1, const n_cusp&m2) const
{
  int ans=0;
#ifdef DEBUG_CUSP_EQ
  cout<<"Testing equivalence of cusps " << m1 << " and " << m2;
  cout<<" (N="<< modulus <<")"<<endl;
  cout << "Corr. infmats: " << (abstract_mat)m1 << "   "<<(abstract_mat)m2 << endl;
#endif
  long cl1=m1.cl();
  if (cl1==m2.cl())
    {

      Quad x1 = m1.c() * m2.d();
      Quad x2 = m1.d() * m2.c();
      Quad x3 = m1.c() * m2.c();

      if (cl1!=0)
	{
	  long cf = geometry::aitch2_nonprinc_det_val();
	  
	  // debug code:
	  if (ndiv(cf,x1)||ndiv(cf,x2)||ndiv(cf,x3))
	    {
	      cerr << "Error: wrong det in cusplist::eq !" <<endl;
	      cerr << "x1 = " << x1 << "   x2 = " << x2 << "   x3 = "<<x3;
	      cerr << "   --- but cf = " << cf << endl;
	    }
	  
	  x1/=cf; x2/=cf; x3/=cf;
	}
#ifdef DEBUG_CUSP_EQ
      cout << "x1= "<<x1 << "   x2= "<<x2 << "   x3= "<<x3 <<endl;
#endif

      Qideal test(x3);
      test += modulus;
      
#ifdef DEBUG_CUSP_EQ
      cout << "<x3>+modulus = " << test << endl;
#endif
      
      long imax=Field::nunits;
      Quad coeff=Field::fundunit;
      if (GL_vs_SL_flag==0)
	{
	  imax/=2;
	  coeff*=coeff;
	}
      for (long i=0; i<imax; i++)
	{
	  ans = ans || test.divides(x2-x1);
	  x2*=coeff;
	}
    }      
#ifdef DEBUG_CUSP_EQ
      cout<<"Returning equiv="<< ans<<endl;   
      cout << endl;
#endif
  return ans;
}

long cusplist::index(const n_cusp& c)
{  
  // adds c to list if not there already, and return index
  long ans=-1;
  for (long i=0; (i<num) && (ans<0); i++) if (eq(c,list[i]))  ans=i;
  if (ans==-1)
    if (num<maxnum)
      {
#ifdef DEBUG_CUSP_EQ
      cout<< "Adding cusp "<< c << " to cusplist"<<endl;   
      cout << endl;
#endif
	list[num]=c;
	ans=num;
	num++;
      }
    else { cerr << "Error: Too many cusps in cusplist!" << endl; }
  return ans;
}

// END OF FILE symb.cc
