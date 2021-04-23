// FILE field.h

//////////////////////////////////////////////////////////////////////
// All data elements of Class field are static.  To initialize,
// include the line
//    Field::init(d);
// at the top of the program, where d is a suitable positive integer.
//////////////////////////////////////////////////////////////////////

#ifndef __FIELD_H__
#define __FIELD_H__
#include <iostream>

class Field {
    static long     h;          // class number (0 if Euclidean)
public:
    static long     d;          // square-free >0
    static long     t;          // trace of w
    static long     n;          // norm of w
    static long     disc;       // discriminant
    static long   nunits;       // number of units
    static Quad fundunit;       // fundamental unit
    static double  rootd;       // sqrt(d)
    static double rootdisc;     // sqrt(-disc)
    static char     name;       // name of w for printing

    static void init(long dd);       // Luiz had a c'tor, but that's bad C++
// next line for debugging
// ~Field() { cerr<< "Destructor ~Field executed!\n" ; }
    static int is_euc() { return (h==0); }
    static int is_princ() { return (h<2); }
    static int not_princ() { return !is_princ(); }
    static long classnum() { if (is_euc()) return 1; else return h; }
    static void display(ostream& s = cout);
};

#endif

// END OF FILE field.h
