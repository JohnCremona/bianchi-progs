// FILE qidl_labels.cc  -- output a list of labels of ideals in a range of norms

// for use in scripts to with minimal interface

#include <iostream>
#include "qidloop.h"
#include "mat22.h"
#include "primes.h"
int main (int argc, char *argv[])
{
  if ( argc < 4 ) { puts ("Usage: qidl_labels d min_norm max_norm"); return 0; }
  long d = atoi(argv[1]);
  long n1 = atoi(argv[2]);
  long n2 = atoi(argv[3]);

  Quad::field(d);

  Qidealooper loop_both(n1, n2, /* both conjugates */ 1, /* sorted */ 1);
  while( loop_both.not_finished() )
    {
      Qideal I = loop_both.next();
      cout << ideal_label(I) << " ";
    }
  cout << endl;
}


// END OF FILE qidltest.cc
