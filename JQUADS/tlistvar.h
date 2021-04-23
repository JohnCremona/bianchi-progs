// FILE tlistvar.h

// provides two template classes:
//   1)   Tlist<T>   --  a "list" of objects of type T
//   2)   Tvar<T>    --  a class for looping through a Tlist<T>

// The prototype classes were John's Quadlist and Quadvar classes.
// To use a Quadlist or Quadvar, include the following lines:
//     typedef Tlist<Quad> Quadlist;
//     typedef Tvar<Quad> Quadvar;

#ifndef __TLISTVAR_H__
#define __TLISTVAR_H__

template <class T> class Tlist;
template <class T> class Tvar;

////////////////////////////////////
// Replacement for class Quadlist //
////////////////////////////////////

// Tlist<T> requires there to exist
//   class T with default constructor T::T()
//   void T::operator=(const T&)      // or default assignment operator for T
//   int T::operator==(const T&) const
//   friend ostream& operator<<(ostream&, const T&)

template <class T> class Tlist {
  long length;            // made private by JSB;  access with getlength()
  T* items;  
public:
  long getlength() const {return length;}

  //constructors
  Tlist(long n=0) {if (n<0) n=0; length=n; items=new T[n];}
  ~Tlist() {delete[] items;}   // surely no point in setting length=0
  //Grab an *existing* (newed) array of n t-objs and regard it as a Tlist
  //JSB: don't think this is ever called
  Tlist(long n, T* list) {length=n; items=list; cerr<<"Funny Tlist ctor called!\n";} //NB Does NOT copy!
  Tlist(const Tlist&l) 
    {items=new T[length=l.length]; 
     long n=length; T* x=items, *y=l.items;
     while(n--) {*x++=*y++;}}
  Tlist(const Tlist&l, long n) 
    {length= n<=l.length ? n : l.length;
     items=new T[length]; 
     n=length; T* x=items, *y=l.items;
     while(n--) {*x++=*y++;}}
  Tlist& operator=(const Tlist&l)
    {if (this==&l) return *this;
     delete[] items; items=new T[length=l.length]; 
     long n=length; T* x=items, *y=l.items;
     while(n--) {*x++=*y++;}
     return *this;}

  T get(long i) const {
    if ((i>=0)&&(i<length)) {return items[i];}
    else {cerr<<"bad index "<<i<<" in Tlist!\n"; exit(1); }  // abort program!
  }
  long locate(const T& q)  //return index of q in list (else -1)
    {
      for (long n=0; n<length; n++) if(q==items[n]) return n;
      return -1;
    }
  T& operator[](long i) {
    if ((i>=0)&&(i<length)) {return items[i];}
    else {cerr<<"bad index "<<i<<" in Tlist::[]!\n"; return items[0];}
  }
  T operator()(long i) const {
    if ((i>=0)&&(i<length)) {return items[i];}
    else {cerr<<"bad index "<<i<<" in Tlist::()!\n"; return items[0];}
  }
  void set(long i, const T& val) {
    if ((i>=0)&&(i<length)) items[i]=val; 
    else {cerr<<"bad index "<<i<<" in Tlist::set!\n";}
  }
  void truncate(long newlength) 
    { if (newlength<length) length = newlength;}
  // no need to deallocate memory here because copy constructors do so
  // almost always it is the temporary of a function that calls truncate
  void output(ostream& os) const
    {long i; T *it;
     for(i=0,it=items;i<length;i++,it++)
       {if(i)os<<", ";os<<(*it);}
   }
friend class Tvar<T>;
};

template <class T>
inline ostream& operator<<(ostream& os, const Tlist<T>& l)
{l.output(os);return os;}

///////////////////////////////////
// Replacement for class Quadvar //
///////////////////////////////////

// Tvar<T> requires
//   class T for which Tlist<T> is okay

// also used to need:
//   some conversion  0 -----> T  as the return value if out of range is T(0)
// Now stop with run-time error instead!

template <class T> class Tvar {
public:
        long index;        /* current index */
private:
        T* values;
        long maxindex;     /* max index */
public:
        Tvar(const Tlist<T>& l) 
          {maxindex=l.length-1; index=0; values=l.items;}
        void init(const Tlist<T>& l) 
          {maxindex=l.length-1; index=0; values=l.items;}
        int next() {if ((index++)<maxindex) {values++; return 1;}
                    else return 0;
                   }
/////////////////////////////////////////////////////////////////////////
// JSB: Next comes the prefix increment operator, usage  ++avar
//  For the intended use, shouldn't we define a postfix operator  avar++ ?
//  My reading of "Ellis&Stroustrup: The Annotated C++ Reference Manual"
//  (p338f) suggests that where  avar++  was written, the program had to
//  invoke  ++avar  and would have given the wrong result if, unlike in
//  the usage example below, the return value of  avar++  had ever been used.
//  This is confirmed by experiment: if tmanin (March 1995) is compiled with
//  both operators, only postfix ++ is ever called.  Yee-ha!
/////////////////////////////////////////////////////////////////////////
        T operator++() 
	  { // cout << "Ho-hum.";                  // for debugging
	    if ((index++)<maxindex) {values++; return *values;}
	    else 
	      { cerr << "Error: prefix ++ incrementing out of range!"<<endl;
		exit(1);
	      }
	  }
/////////////////////////////////////////////////////////////////////////
// JSB: Here's my analogous postfix increment operator, usage  avar++
//  It always increments  index,  but only increments  values  if it won't
//  be going out of range.
/////////////////////////////////////////////////////////////////////////
	T operator++(int)
	  { // cout << "Yee-Ha!";                  // for debugging
	    if (ok())
	      if ((index++)<maxindex) { return *(values++); }
	      else { return *values; }
	    else
	      { cerr << "Error: postfix ++ incrementing out of range!"<<endl;
		exit(1);
	      }
	  }
        int ok() const {return index<=maxindex;}
        int more() const {return index<maxindex;}
        T value() const {return *values;}
        operator T() const {return *values;}
};

/* Usage of Tvar<T>: to loop through a Tlist<T> "alist" (which might be empty):
   T a; either of

   for(Tvar<T> avar(alist); avar.ok(); avar++) {a=avar; ... ;}
   for(avar.init(alist); avar.ok(); avar++) {a=avar; ... ;}  
   // second form IFF avar already exists

   Either way, the argument, alist, must be an instance of an existing
   Tlist<T>, and must not be a function returning a Tlist<T> (else
   unflagged error will follow)
   
*/


///////////////////////////////////////////////////////
/// JEC's original Quadlist and Quadvar classes ///////
///////////////////////////////////////////////////////
/******************************************************
class Quadlist {
public:
  long length;        
private:
  Quad *items;  
public:
  Quadlist(long n=0) {length=n; items=new Quad[n];}
  ~Quadlist() {length=0; delete[] items;}
  Quadlist(long n, Quad* list) {length=n; items=list;} //NB Does NOT copy!
  Quadlist(const Quadlist&l) 
    {items=new Quad[length=l.length]; 
     long n=length; Quad* x=items, *y=l.items;
     while(n--) {*x++=*y++;}}
  Quadlist(const Quadlist&l, long n) 
    {length= n<=l.length ? n : l.length;
     items=new Quad[length]; 
     n=length; Quad* x=items, *y=l.items;
     while(n--) {*x++=*y++;}}
  Quadlist& operator=(const Quadlist&l)
    {if (this==&l) return *this;
     delete[] items; items=new Quad[length=l.length]; 
     long n=length; Quad* x=items, *y=l.items;
     while(n--) {*x++=*y++;}
     return *this;}
  Quad get(long i) const {
    if ((i>=0)&&(i<length)) {return items[i];}
    else {cerr<<"bad index "<<i<<" in Quadlist!\n"; return 0;}
  }
  long locate(const Quad& q)  //return index of q in list (else -1)
    {
      for (long n=0; n<length; n++) if(q==items[n]) return n;
      return -1;
    }
  Quad& operator[](long i) {
    if ((i>=0)&&(i<length)) {return items[i];}
    else {cerr<<"bad index "<<i<<" in Quadlist!\n"; return items[0];}
  }
  Quad operator()(long i) const {
    if ((i>=0)&&(i<length)) {return items[i];}
    else {cerr<<"bad index "<<i<<" in Quadlist!\n"; return items[0];}
  }
  void set(long i, const Quad& val) {
    if ((i>=0)&&(i<length)) items[i]=val; 
    else {cerr<<"bad index "<<i<<" in Quadlist!\n";}
  }
  void truncate(long newlength) 
    { if (newlength<length) length = newlength;}
  // no need to deallocate memory here because copy constructors do so
  // almost always it is the temporary of a function that calls truncate
  void output(ostream& os) const
    {long i;Quad *it;
     for(i=0,it=items;i<length;i++,it++)
       {if(i)os<<", ";os<<(*it);}
   }
  friend ostream& operator<<(ostream& os, const Quadlist& l);
  friend class Quadvar;
};

inline ostream& operator<<(ostream& os, const Quadlist& l)
    {l.output(os);return os;}

class Quadvar {
public:
  long index;        // current index
private:
  Quad* values;
  long maxindex;     // max index 
public:
  Quadvar(const Quadlist& l) 
    {maxindex=l.length-1; index=0; values=l.items;}
  void init(const Quadlist& l) 
    {maxindex=l.length-1; index=0; values=l.items;}
  long next() {if ((index++)<maxindex) {values++; return 1;}
  else return 0;
	     }
  Quad operator++() {if ((index++)<maxindex) {values++; return *values;}
  else return 0;
		   }
  long ok() const {return index<=maxindex;}
  long more() const {return index<maxindex;}
  Quad value() const {return *values;}
        operator Quad() const {return *values;}
};

**************************************************/

/* Usage of Quadvar: to loop through a Quadlist "alist" (which might be empty):
   Quad a; either of

   for(Quadvar avar(alist); avar.ok(); avar++) {a=avar; ... ;}
   for(avar.init(alist); avar.ok(); avar++) {a=avar; ... ;}  
   // second form IFF avar already exists

   Either way, the argument, alist, must be an instance of an existing
   Quadlist, and must not be a function returning a Quadlist (else
   unflagged error will follow)
   
*/

#endif

// END OF FILE tlistvar.h
