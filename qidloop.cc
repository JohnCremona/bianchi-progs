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


// END OF FILE QIDLOOP.CC
