#include "qidloop.h"
#include "homspace.h"
#define LOOPER
//#define CHECK_CONJUGATE

int main ()
{
  long d, maxpnorm(1000);
  cerr << "Enter field: " << flush;  cin >> d;
  if (!check_field(d))
   {
     cerr<<"field must be one of: "<<valid_fields<<endl;
     exit(1);
   }

 vector<scalar> charlist;
 charlist.push_back(scalar(0));
 cerr << "Enter list of characteristics, ending in 0: " << flush;
 scalar ch(1); // any nonzero
 while(ch!=0)
   {
     cin >> ch;
     if (ch!=0) charlist.push_back(ch);
   }

 Quad::field(d,maxpnorm);
 int plusflag=1, cuspidalflag=1;
 cerr << "Plus space? "; cin>>plusflag;
 cerr << "cuspidal subspace? "; cin>>cuspidalflag;
#ifdef LOOPER
 long firstn, lastn;
 int both_conj;
 cerr<<"Both conjugates? (0/1) "; cin >> both_conj;
 cerr<<"Enter first and last norm for level: ";
 cin >> firstn >> lastn;
 cerr<<endl;
 cout<<"Level";
 for ( const auto& c : charlist)
   cout<<" "<<c;
 cout<<endl;

 Qidealooper loop(firstn, lastn, both_conj, 1); // sorted within norm
 while( loop.not_finished() )
   {
     Qideal N = loop.next();
#else
     Qideal N;
     while(cerr<<"Enter level (ideal label or generator): ", cin>>N, !N.is_zero())
       {
#endif
         cout << ideal_label(N);
         for ( const auto& c : charlist)
           {
             homspace h(N, plusflag, 0, c);  //level, plusflag, verbose, characteristic
             long dim = (cuspidalflag? h.h1cuspdim(): h.h1dim());
             cout << " " << dim;
           }
         cout<<endl;
       }
}
