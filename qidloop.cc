// FILE QIDLOOP.CC

#include "qidloop.h"

Qidealooper::Qidealooper(long nmin, long nmax, int both_conjugates, int sorted)
  : n(nmin), maxn(nmax), both(both_conjugates), sort(sorted)
{
  vector<Qideal> Ilist = (sort? Qideal_lists::ideals_with_norm(n): ideals_with_norm(n));
  I_norm_n.insert(I_norm_n.end(), Ilist.begin(), Ilist.end());
  advance();
}

Qideal Qidealooper::next() // must only call when not not_finished
{
  Qideal I = I_norm_n.front();
  I_norm_n.pop_front();
  advance();
  return I;
}

void Qidealooper::advance()
{
  if (! I_norm_n.empty()) return;
  while (I_norm_n.empty() && n<maxn)
    {
      n++;
      vector<Qideal> Ilist = (sort? Qideal_lists::ideals_with_norm(n): ideals_with_norm(n));
      I_norm_n.insert(I_norm_n.end(), Ilist.begin(), Ilist.end());
    }
  // now either n>maxn and I_norm_n is empty, or n<=maxn and it is
  // not.  In the first case, not_finished() will be false
}

vector<Qideal> Quad::class_group;

void Quad::fill_class_group()
{
  long MB = floor(2*sqrt(disc)/PI);
  class_group.push_back(Qideal());
  if(class_number==1) return;
  Qidealooper loop(1, MB, 1, 1);
  while( loop.not_finished() )
    {
      Qideal I = loop.next();
      if (I.is_principal()) continue;
      int new_class=1;
      for (vector<Qideal>::iterator Ji = class_group.begin(); (Ji != class_group.end()) && new_class; ++Ji)
        if (I.is_equivalent(*Ji)) new_class = 0;
      if (new_class)
        {
          //          cout << I << " is in a new ideal class (#" << class_group.size() << ")" << endl;
          class_group.push_back(I);
          Qideal I2 = I.conj();
          if (!I.is_equivalent(I2))
            {
              //              cout << I2 << " is in a new ideal class (#" << class_group.size() << ")" << endl;
              class_group.push_back(I2);
            }
        }
    }
  class_number = class_group.size();
  //  cout << "Class number = " << class_number << " with representatives " << class_group << endl;
}


// END OF FILE QIDLOOP.CC
