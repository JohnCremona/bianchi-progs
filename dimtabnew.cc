#include "qidloop.h"
#include "newforms.h"
//#define MODP
//#define LOOPER

int main ()
{
  long d, max(1000);
  cerr << "Enter field (one of "<<valid_fields<<"): " << flush;  cin >> d;
  if (!check_field(d))
    {
      cerr<<"field must be one of: "<<valid_fields<<endl;
      exit(1);
    }
  long ch=0;
#ifdef MODP
  cerr << "Enter characteristic (0 or prime): " << flush;  cin >> ch;
#endif

#ifdef LOOPER
  long firstn, lastn; Quad n;
  int both_conj;
  cerr<<"Both conjugates? (0/1) "; cin >> both_conj;

  cerr<<"Enter first and last norm for Quad loop: ";
  cin >> firstn >> lastn;
  cerr<<endl;
#endif

  Quad::field(d,max);


 cout << "# Table of dimensions of ";
 if (ch) cout<<"mod "<<ch<<" ";
 cout<<"weight 2 Bianchi cuspidal forms and newforms for GL2 over Q(sqrt(-"<<d<<"))" << endl;
 if (Quad::class_group_2_rank>0)
   cout<<"# (with trivial character, and including unramified quadratic twists)"<<endl;
 cout << "# Field\tWeight\tLevel\t";
 cout << "dim(all)\tdim(cuspidal)\tdim(cuspidal, new)" << endl;

#ifdef LOOPER
 Qidealooper loop(firstn, lastn, both_conj, 1); // sorted within norm
 while( loop.not_finished() )
   {
     Qideal N = loop.next();
#else
     Qideal N;
 while(cerr<<"Enter level (ideal label or generator): ", cin>>N, !N.is_zero())
   {
#endif
     cout << "\t"<< d << "\t2\t";                  // field and weight
     cout << ideal_label(N)<<"\t\t"; // level
     homspace hplus(N, 1, 1, 0, ch);  //level, plusflag, cuspidal, verbose
     int dimcusp = hplus.bianchi_form_dimension(1);
     cout << dimcusp << "\t\t";

     int dimcuspnew=0;
     newforms nfdata(N, 0, ch);
     if (nfdata.read_from_file())
       {
         if (Quad::class_group_2_rank == 0)
           {
             dimcuspnew = nfdata.n1ds + nfdata.n2ds;
           }
         else
           {
             dimcuspnew = 2*(nfdata.new1dims[0]+nfdata.new2dims[0]);
             for (auto tdim = nfdata.new1dims.begin()+1; tdim!=nfdata.new1dims.end(); tdim++)
               dimcuspnew += *tdim;
             for (auto tdim = nfdata.new2dims.begin()+1; tdim!=nfdata.new2dims.end(); tdim++)
               dimcuspnew += *tdim;
             dimcuspnew <<= (Quad::class_group_2_rank - 1);
           }
       }
     else
       {
         cout<<"no newforms file! ";
       }
     cout << dimcuspnew << endl;
   }
}
