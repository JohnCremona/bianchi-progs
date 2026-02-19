// FILE qidl_labels.cc  -- output a list of labels of ideals in a range of norms

// for use in scripts to with minimal interface

#include <iostream>
#include "qidloop.h"
#include "mat22.h"
#include "primes.h"
int main (int argc, const char *argv[])
{
  if ( (argc < 4) || (argc > 6) ) { puts ("Usage: qidl_labels d min_norm max_norm <n-per-line> <both_conj>"); return 0; }
  long d = atoi(argv[1]);
  long n1 = atoi(argv[2]);
  long n2 = atoi(argv[3]);
  int n_per_line = 0; // default, no newlines except at the end
  if (argc > 4)
    n_per_line = atoi(argv[4]);
  int both_conj = 1; // default, both conjugates
  if (argc > 5)
    both_conj = atoi(argv[5]);

  Quad::field(d);

  Qidealooper loop_both(n1, n2, both_conj, /* sorted */ 1);
  int n=0;
  while( loop_both.not_finished() )
    {
      Qideal I = loop_both.next();
      cout << label(I);
      n++;
      if ((n < n_per_line) || (n_per_line == 0))
        cout << " ";
      else
        {
          cout<<endl;
          n = 0;
        }
    }
  if (n)
    cout << endl;
}


// END OF FILE qidltest.cc
