#include "qidloop.h"
#include "homspace.h"
#define LOOPER
//#define CHECK_CONJUGATE

int is_B_smooth(const INT& N, const INT& B)
{
  vector<INT> pr = pdivs(N);
  // cout << "N = "<<N<<" has prime factors "<<pr<<" and is ";
  int ans = std::all_of(pr.begin(), pr.end(), [B](const INT& p){return p<=B;});
  // cout << (ans?"":" NOT ")<< B<<"-smooth"<<endl;
  return ans;
}

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
 cerr << "Include characteristic 0? ";
 int char0 = 0;
 cin >> char0;
 if (char0)
   charlist.push_back(scalar(0));
 cerr << "Enter list of characteristics, ending in either 0: " << flush;
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

 cerr << "Enter smoothness parameter (or 0 for none): ";
 INT B;
 cin >> B;
 // Skip levels whose norm is not B-smooth

 Qidealooper loop(firstn, lastn, both_conj, 1); // sorted within norm
 while( loop.not_finished() )
   {
     Qideal N = loop.next();
     if (!is_B_smooth(N.norm(), B))
       continue;
#else
     Qideal N;
     while(cerr<<"Enter level (ideal label or generator): ", cin>>N, !N.is_zero())
       {
#endif
         cout << label(N);
         for ( const auto& c : charlist)
           {
             scalar modulus = (c==0? default_modulus<scalar>(): c);
             homspace h(N, modulus, plusflag, 0, c);  //level, modulus, plusflag, verbose, characteristic
             long dim = (cuspidalflag? h.h1cuspdim(): h.h1dim());
             cout << " " << dim;
           }
         cout<<endl;
       }
}
