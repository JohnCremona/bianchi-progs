#include "qidloop.h"
#include "newforms.h"
//#define MODP

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

  long firstn=1, lastn; Quad n;
  int both_conj=1;
  //cerr<<"Both conjugates? (0/1) "; cin >> both_conj;

  cerr<<"Enter max norm for Quad loop: ";
  cin >> lastn;
  cerr<<endl;

  Quad::field(d,max);
  int nchi = 1<<(Quad::class_group_2_rank);
  int nchi2 = nchi/2; // only used when nchi is even

  cout << "# Table of dimensions of ";
  if (ch) cout<<"mod "<<ch<<" ";
  cout<<"weight 2 Bianchi cuspidal forms and newforms for GL2 over Q(sqrt(-"<<d<<"))" << endl;
  if (Quad::class_group_2_rank>0)
    cout<<"# (with trivial character, and including unramified quadratic twists)"<<endl;
  cout << "# Field\tWeight\tLevel\t";
  cout << "dim(all)\tdim(cuspidal)\tdim(cuspidal, new)" << endl;

  // dictionary with key-level label, value = list of new dimensions
  // by self-twist character
  map<string, vector<int>> dims;

  Qidealooper loop(firstn, lastn, both_conj, 1); // sorted within norm
  while( loop.not_finished() )
   {
     Qideal N = loop.next();
     string level_label = ideal_label(N);
     vector<Qideal> Ndivisors = alldivs(N);
     // cout<<"Genus characters of divisors of "<<level_label<<":"<<endl;
     // cout << "\t\t"<<Quad::all_disc_factors <<endl;
     // for (auto I=Ndivisors.begin(); I!=Ndivisors.end(); ++I)
     //   {
     //     cout<<ideal_label(*I)<<" --> \t";
     //     for (auto d = Quad::all_disc_factors.begin(); d!=Quad::all_disc_factors.end(); ++d)
     //       cout<<I->genus_character(*d)<<" ";
     //     cout<<endl;
     //   }
     cout << "\t"<< d << "\t2\t";                  // field and weight
     cout << level_label<<"\t\t"; // level

     vector<int>& newdims = dims[level_label];
     newforms nfdata(N, 0, ch);
     if (nfdata.read_from_file())
       {
         newdims.resize(nchi, 0);
         if (nchi == 1)
           {
             newdims[0] = nfdata.n1ds + nfdata.n2ds;
             //cout<<"newdims="<<newdims<<endl;
           }
         else
           {
             newdims[0] = nchi*(nfdata.new1dims[0] + nfdata.new2dims[0]);
             for (int i=1; i<nchi; i++)
               newdims[i] = nchi2*(nfdata.new1dims[i] + nfdata.new2dims[i]);
             //cout<<"newdims="<<newdims<<endl;
           }
         int dimcuspnew = std::accumulate(newdims.begin(), newdims.end(), 0, std::plus<int>());
         int dimcuspold = 0;
         for (auto Di=Ndivisors.begin(); Di!=Ndivisors.end(); ++Di)
           {
             Qideal D = *Di;
             if (N==D)
               break;
             vector<int> olddims = old_multiplicities(D, dims[ideal_label(D)], N);
             int dimcuspold1 = std::accumulate(olddims.begin(), olddims.end(), 0, std::plus<int>());
             // if (dimcuspold1)
             //   cout<<"["<<ideal_label(D)<<":"<<dimcuspold1<<"]";
             dimcuspold += dimcuspold1;
           }
         int dimcusp = dimcuspold + dimcuspnew;
         cout << dimcusp << "\t\t";
         cout << dimcuspnew << endl;
       }
     else
       {
         cout<<"no newforms file for level " << level_label<<"! "<<endl;
         break;
       }
   }
}
