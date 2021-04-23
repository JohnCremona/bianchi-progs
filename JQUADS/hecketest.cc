// FILE hecketest2.cc  ---  test for Hecke operators (version 2)

#include "homspace.h"
#include "qidloop.h"

int verbose;
int plusflag=1;

class heck_test_obj : homspace {
  int verbose;
  Primelist plist_bad, plist_goodprin, plist_goodnonprin;// initialised by ctor
  mat chi0;                                           // initialised by ctor
  void make_plists();
  void proc_Tp_solo();
  void proc_Tp_loop();
  void proc_Uprinc_loop();
  void proc_Tpp_solo();
  void proc_Tpp_loop();
  void proc_Tpp_loop_princ();
  void proc_Tpq_solo();
  void proc_Tpq_loop();
  void proc_Tppchi_loop();
  void proc_Tpqchi_loop();
  void proc_chi_loop_np();
  void proc_cx_conj();
  void menu();
public:
  heck_test_obj(const Qideal&level, int vflag=0);
};

void init()
{
  long d;
  cout << "Enter field: " << flush;  cin >> d;
  Field::init(d);
  geometry::init();
  Quadprimes::init(1000);
  cout << "Verbose? (0/1) "; cin >> verbose;
  cout << endl;
  if (verbose)
    { Field::display(cout);
      geometry::display(cout,verbose);
    }
}

void get_and_test_level()
{
 Qideal alpha,zero=0;
 while(cout<<"Enter level (1 0 0 to finish): ", cin>>alpha, alpha!=zero)
   { heck_test_obj h(alpha, verbose); }
}

int main ()
{
  cout << endl << "TEST PROGRAM FOR HECKE MATRICES (Version 2)" << endl;
  init();
  get_and_test_level();
  geometry::dealloc();
}

void heck_test_obj::make_plists()
{
  long* flagarray=new long[Quadprimes::list.getlength()];
  long n_bad=0, n_goodprin=0, n_goodnonprin=0, n_all=0;
  long currflag;
  for (Primevar pr(Quadprimes::list); pr.ok(); pr++)
    {
      Quadprime p=pr;
      if (p.divides(modulus))
	{ n_bad++; currflag=0; }
      else
	{
	  if (p.isprincipal())
	    { n_goodprin++; currflag=1; }
	  else
	    { n_goodnonprin++; currflag=2; }
	}
      flagarray[n_all++]=currflag;
    }
  if (n_bad+n_goodprin+n_goodnonprin != n_all)
    { cerr << "Error: miscount in make_plists; aborting."<<endl; exit(1);}
  plist_bad = Primelist(n_bad);
  plist_goodprin = Primelist(n_goodprin);
  plist_goodnonprin = Primelist(n_goodnonprin);
  n_bad=0;
  n_goodprin=0;
  n_goodnonprin=0;
  n_all=0;
  for (long i=0; i<Quadprimes::list.getlength(); i++)
    switch(flagarray[i])
      {
      case 0:
	plist_bad[n_bad++]=Quadprimes::list[i];
	break;
      case 1:
	plist_goodprin[n_goodprin++]=Quadprimes::list[i];
	break;
      case 2:
	plist_goodnonprin[n_goodnonprin++]=Quadprimes::list[i];
      }
  if (verbose)
    { cout << "The bad primes: " << plist_bad << endl << endl;
      cout << "The good principal primes: " << plist_goodprin << endl<< endl;
      cout << "The good non-principal primes: "<<plist_goodnonprin<<endl<<endl;
    }
}
	
heck_test_obj::heck_test_obj(const Qideal&level, int vflag) : homspace(level, plusflag, 0)
{
  make_plists();
  verbose=vflag;

  long d = h1dim();
  long den = h1denom();
  cout << "Dimension = " << d << endl;
//  if(den!=1)
  cout << "denominator = " << den << endl;

  Quadprime p0 = plist_goodnonprin(0);

// the normaliser character
  chi0 = chi_p(p0);
  if (verbose)
    {
      cout << "The normaliser character computed via prime " << p0 << endl;
      cout << chi0;
    }
// compute complex conjugation where this makes sense
  if (level==level.conj())
    {
      cout << "Self-conjugate level: ";
      proc_cx_conj();
    }

  menu();
}

void heck_test_obj::menu()
{
  long ch;
  cout << endl << "Options for analysing level "<< modulus << endl;
  cout << " 0) quit" << endl;
  cout << " 1) T(p) for given p" << endl;
  cout << " 2) T(p) for loop of p" << endl;
  cout << " 3) T(p^2) for given p" << endl;
  cout << " 4) T(p^2) for loop of p" << endl;
  cout << " 5) T(pq) for given p,q" << endl;
  cout << " 6) T(pq) for loop of p,q" << endl;
  cout << " 7) principal good primes of norm < 50" << endl;
  cout << " 8) hecke_npnp for *all* (sic) good primes of norm < 50" << endl;
  cout << " 9) chi_p for loop of non-principal p" << endl;
  cout << "10) T(p^2) for loop of principal (sic!) p" << endl;
  cout << "11) T(p^2) with chi-adjust, for loop of p" << endl;
  cout << "12) T(pq) with chi-adjust, for loop of p,q" << endl;
  cout << "13) complex conjugation (may fail if level not self-conj)" << endl;
  cout << "14) experimental U(p) for bad princ p" << endl;
  do {
    do { cout << endl << "Choice of proc: ";  cin >> ch; }
    while ((ch<0)||(ch>14));
    cout << endl;
    switch(ch)
      {
      case 1:
	proc_Tp_solo();
	break;
      case 2:
	proc_Tp_loop();
	break;
      case 3:
	proc_Tpp_solo();
	break;
      case 4:
	proc_Tpp_loop();
	break;
      case 5:
	proc_Tpq_solo();
	break;
      case 6:
	proc_Tpq_loop();
	break;
      case 7:
	cerr << "Not implemented." << endl;
	break;
      case 8:
	cerr << "Not implemented." << endl;
	break;
      case 9:
	proc_chi_loop_np();
	break;
      case 10:
	proc_Tpp_loop_princ();
	break;
      case 11:
	proc_Tppchi_loop();
	break;
      case 12:
	proc_Tpqchi_loop();
	break;
      case 13:
	proc_cx_conj();
	break;
      case 14:
	proc_Uprinc_loop();
	break;
      }
  }
  while (ch!=0);
}

void heck_test_obj::proc_Tp_solo()
{
  Quadprime p;
  cout << "T(p) for principal p --- enter p: ";
  cin >> p;
  cout << "T(p) " << p << endl;
  mat m = hecke_p(p);
  cout << m << endl;
}

void heck_test_obj::proc_Tp_loop()
{ // principal primes
  long imin, imax, max_imax=plist_goodprin.getlength();
  cout << "Computing (unadjusted) T(p) ..." << endl;
  cout << "Number of good princ primes available: " << max_imax << endl;
  cout << "Enter index (min 1) for start of loop: ";
  cin >> imin;
  cout << "Enter index (max " << max_imax << ") for end of loop: ";
  cin >> imax;
  cout << endl;
  for (long i=imin-1; i<imax; i++)
    {
      cout << "T(p) " << plist_goodprin(i) << endl;
      mat temp=hecke_p(plist_goodprin(i));
      cout << temp << endl;
    }
}

// Added 20-07-1998
void heck_test_obj::proc_Uprinc_loop()
{ // bad principal primes
  long imax=plist_bad.getlength();
  cout << "Computing U(p) for bad princ p ..." << endl;
 
  for (long i=0; i<imax; i++)
    {
      cout << "U(p) " << plist_bad(i) << endl;
      if (plist_bad(i).isprincipal())
	{
	  mat temp=hecke_u_princ(plist_bad(i));
	  cout << temp << endl;
	}
      else
	{
	  cout << "... prime not principal" << endl;
	}
    }
}

void heck_test_obj::proc_Tpp_solo()
{
  cerr<< "Not implemented!" << endl;
}

void heck_test_obj::proc_Tpp_loop()
{ // squares of non-principal primes
  long imin, imax, max_imax=plist_goodnonprin.getlength();
  cout << "Computing unadjusted T(p^2) ..." << endl;
  cout << "Number of good non-princ primes available: " << max_imax << endl;
  cout << "Enter index (min 1) for start of loop: ";
  cin >> imin;
  cout << "Enter index (max " << max_imax << ") for end of loop: ";
  cin >> imax;
  cout << endl;
  for (long i=imin-1; i<imax; i++)
    {
      cout << "T(p^2)nochi " << plist_goodnonprin(i) << endl;
      mat temp=hecke_npnp(plist_goodnonprin(i));
      cout << temp << endl;
    }
}

void heck_test_obj::proc_Tpp_loop_princ()
{ // squares of principal (sic!) primes
  long imin, imax, max_imax=plist_goodprin.getlength();
  cout << "Computing T(p^2) for principal (sic!) primes: " << endl;
  cout << "Number of good princ primes available: " << max_imax << endl;
  cout << "Enter index (min 1) for start of loop: ";
  cin >> imin;
  cout << "Enter index (max " << max_imax << ") for end of loop: ";
  cin >> imax;
  cout << endl;
  for (long i=imin-1; i<imax; i++)
    {
      cout << "T(p^2)nochi " << plist_goodprin(i) << endl;
      mat temp=hecke_npnp(plist_goodprin(i));
      cout << temp << endl;
    }
}

void heck_test_obj::proc_Tppchi_loop()
{ // squares of non-principal primes
  long imin, imax, max_imax=plist_goodnonprin.getlength();
  Quadprime r = plist_goodnonprin(0);

  cout << "Computing T(p^2) with adjustment by chi_r: " << endl;
  cout << "Number of good non-princ primes available: " << max_imax << endl;
  cout << "Enter index (min 1) for start of loop: ";
  cin >> imin;
  cout << "Enter index (max " << max_imax << ") for end of loop: ";
  cin >> imax;
  cout << endl;

  for (long i=imin-1; i<imax; i++)
    {
      cout << "T(p^2) " << plist_goodnonprin(i) << endl;
      mat temp=hecke_npnpchi(plist_goodnonprin(i),r);
      cout << temp << endl;
    }
}

void heck_test_obj::proc_Tpq_solo()
{ // product of two distinct non-principal primes
  long i, j, max_i=plist_goodnonprin.getlength();
  cout << "Number of good non-princ primes available: " << max_i << endl;

  cout << "Enter index (min 1) for p: ";
  cin >> i; i-=1;
  Quadprime p = plist_goodnonprin(i);
  cout << "Taking p equal to " << p << endl;

  cout << "Enter index for q: ";
  cin >> j; j-=1;
  Quadprime q = plist_goodnonprin(j); 
  cout << "Taking q equal to " << q << endl;

  cout << "T(pq) " << p << " " << q << endl;
  mat temp = hecke_npnq(p,q);
  cout << temp << endl;
}

void heck_test_obj::proc_Tpq_loop()
{ // products of two distinct non-principal primes
  long imin, imax, jmin, jmax, max_imax=plist_goodnonprin.getlength();
  cout << "Computing T(pq) ..." << endl;
  cout << "Number of good non-princ primes available: " << max_imax << endl;
  cout << "Enter index (min 1) for start of p-loop: ";
  cin >> imin;
  cout << "Starting p-loop at " << plist_goodnonprin(imin-1) << endl;
  cout << "Enter index (max " << max_imax << ") for end of p-loop: ";
  cin >> imax;
  cout << "Ending p-loop at " << plist_goodnonprin(imax-1) << endl;
  cout << "Enter index for start of inner loop (q-loop): ";
  cin >> jmin;
  cout << "Enter index (max " << max_imax << ") for end of inner loop (q-loop): ";
  cin >> jmax;
  cout << endl;
  for(long i=imin-1; i<imax; i++)
    {
      Quadprime p = plist_goodnonprin(i);
      cout<< "Computing T(pq) with p = " << p << endl;
      
      //    j = ((i==0)? 1 : 0);
      long j = ((i+1> jmin-1)? i+1 : jmin-1);
      while (j<jmax)
	{
	  cout << "T(pq)nochi " << p << " " << plist_goodnonprin(j) << endl;
	  mat temp = hecke_npnq(p,plist_goodnonprin(j));
	  cout << temp << endl;
	  do j++; while (j==i);   // avoid i==j if j-loop started below i
	}
    }
}

void heck_test_obj::proc_Tpqchi_loop()
{ // products of two distinct non-principal primes
  long imin, imax, jmin, jmax, max_imax=plist_goodnonprin.getlength();
  cout << "Computing T(pq) with adjustment by chi_r ..." << endl;
  cout << "Number of good non-princ primes available: " << max_imax << endl;
  cout << "Enter index (min 1) for start of p-loop: ";
  cin >> imin;
  cout << "Starting p-loop at " << plist_goodnonprin(imin-1) << endl;
  cout << "Enter index (max " << max_imax << ") for end of p-loop: ";
  cin >> imax;
  cout << "Ending p-loop at " << plist_goodnonprin(imax-1) << endl;
  cout << "Enter index for start of inner loop (q-loop): ";
  cin >> jmin;
  cout << "Enter index (max " << max_imax << ") for end of inner loop (q-loop): ";
  cin >> jmax;
  Quadprime r = plist_goodnonprin(0);

  for(long i=imin-1; i<imax; i++)
    {
      Quadprime p = plist_goodnonprin(i);
      cout<< "Computing T(pq) with p = " << p << endl;
      
      //    j = ((i==0)? 1 : 0);
      long j = ((i+1> jmin-1)? i+1 : jmin-1);
      while (j<jmax)
	{
	  cout << "T(pq) " << p << " " << plist_goodnonprin(j) << endl;
	  mat temp = hecke_npnqchi(p,plist_goodnonprin(j),r);
	  cout << temp << endl;
	  do j++; while (j==i);   // avoid i==j if j-loop started below i
	}
    }
}

void heck_test_obj::proc_chi_loop_np()
{ // compute chi(p) for loop of p
  long imin, imax, max_imax=plist_goodnonprin.getlength();
  cout << "Computing chi(p) ... " << endl;
  cout << "Number of good non-princ primes available: " << max_imax << endl;
  cout << "Enter index (min 1) for start of loop: ";
  cin >> imin;
  cout << "Enter index (max " << max_imax << ") for end of loop: ";
  cin >> imax;
  cout << endl;
  for (long i=imin-1; i<imax; i++)
    {
      cout << "Chi() " << plist_goodnonprin(i) << endl;
      mat temp=chi_p(plist_goodnonprin(i));
      cout << temp << endl;
      if (temp!=chi0)
	{ cerr << "Warning: chi differs from chi0 !" << endl;}
    }
}

void heck_test_obj::proc_cx_conj()
{
  cout << "Computing effect of complex conjugation ... "<<endl;
  mat temp=calc_conj();
  cout << temp << endl;
}

/*
void n_50(const Qideal&n);            // good princ primes up to norm 50
void pp_50(const Qideal&n);           // hecke_npnp for primes up to norm 50
*/


void pp_50(const Qideal&n)
{
  homspace h(n,plusflag,0);  //level=n, plusflag, verbose=0

  long d = h.h1dim();
  long den = h.h1denom();
  cout << "Dimension = " << d << endl;
  if (d>0)
    {
      cout << "denominator = " << den << endl;

      int goon=1;
      for (Primevar pr(Quadprimes::list);
	   (goon)&&pr.ok();
	   pr++)
	{
	  Quadprime p=pr;
	  if (p.norm() >= 50) goon=0;
	  if ( goon && (!p.divides(n)) && p.isprincipal() )
	    {
	      cout << p << endl;
	      cout << h.hecke_npnp(p) << endl;
	    }
	}
    }
}

/*
// OBSOLETE
// find smallest good non-principal prime
Quadprime smallest_good_non_princ(const Qideal&n)
{
  int found=0;
  Quadprime p;
  Primevar pr(Quadprimes::list);
  while ((found==0)&&pr.ok())
    {
      p=pr;
      if ( (!p.divides(n)) && (!p.isprincipal()) )
	{
	  found=1;
	}
      pr++;
    }
  if (found==0)
    {
      cerr << "Error: Unable to find good non-principal prime." << endl;
      exit(1);
    }
  return p;
}
*/

/*

  Quadlist badprimes = level::plist;
  int nq = badprimes.length; int firstq=0;  // =0 for all W's
  if (d>0)
    {
      mat id = (den*den)*idmat(d);
      mat* wqlist = new mat[nq];
      for (Quadvar qvar=badprimes; qvar.ok(); qvar++)
	{Quad q=qvar.value(); int i=qvar.index; if(i<firstq) continue;
	 cout << "Computing W("<<q<<")...  " << flush;
	 wqlist[i] = h.heckeop(q,verbose);
	 cout << "done. " << flush;
	 if (wqlist[i]*wqlist[i]==id) cout << "Involution!" << "\n";
	 else                         cout << "NOT an involution...." << "\n";
       }
      int np=5,ip=0; 
      mat* tplist = new mat[np];
      Quad p;
      for (Quadvar pr(quadprimes); pr.ok()&&(ip<np); pr++, ip++)
	{while (n%pr.value()==0) pr++;
	 p=pr.value();
	 cout << "Computing T_p for p = " << p << "\n";
	 tplist[ip] = h.heckeop(p,verbose);
	 cout << "done. " << flush;
	 for (int kp=firstq; kp<nq; kp++)
	   {if (wqlist[kp]*tplist[ip]!=tplist[ip]*wqlist[kp])
	      {cout << "Problem: T_p matrix "<<ip<<" and W_q matrix "<<kp<<" do not commute!" << "\n";
	     }
	  }
	 for (int jp=0; jp<ip; jp++)
	   {if (tplist[ip]*tplist[jp]!=tplist[jp]*tplist[ip])
	      {cout << "Problem: T_p matrices "<<ip<<" and "<<jp<<" do not commute!" << "\n";
	     }
	  }
       }
      delete[] wqlist; delete[] tplist;
    }      // end of if(d>0)

}       // end of while()
exit(0);
       // end of main()

*/

// END OF FILE hecketest.cc
