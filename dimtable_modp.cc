#include "qidloop.h"
#include "homspace.h"
#define LOOPER
//#define CHECK_CONJUGATE

int main ()
{
  long d, max(1000);
  cerr << "Enter field (one of "<<valid_fields<<"): " << flush;  cin >> d;
 if (!check_field(d))
   {
     cerr<<"field must be one of: "<<valid_fields<<endl;
     exit(1);
   }

 vector<long> charlist;
 charlist.push_back(0);
 cerr << "Enter list of characteristics, ending in 0: " << flush;
 long ch=1;
 while(ch)
   {
     cin >> ch;
     if (ch) charlist.push_back(ch);
   }

 Quad::field(d,max);
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
 for (vector<long>::const_iterator ch = charlist.begin(); ch!=charlist.end(); ch++)
   cout<<" "<<(*ch);
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
         for (vector<long>::const_iterator ch = charlist.begin(); ch!=charlist.end(); ch++)
           {
             homspace h(N, plusflag, cuspidalflag, 0, *ch);  //level, plusflag, cuspidal, verbose, characteristic
             long dim = (cuspidalflag? h.h1cuspdim(): h.h1dim());
             cout << " " << dim;
           }
         cout<<endl;
       }
}
