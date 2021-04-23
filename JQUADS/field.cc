// FILE field.cc

#include "quadarith.h"

//Declaration of static data members of class Field:
long Field::d;
long Field::t;
long Field::n;
long Field::h;
long Field::disc;
long Field::nunits;
Quad Field::fundunit;
double Field::rootd;
double Field::rootdisc;
char Field::name;

Quadlist Quadunits::ulist;

void Field::init(long dd)
{
  d=dd;
  rootd=sqrt(d);
  if ((d+1)%4)
    {
      t=0; disc=-4*d; n=d;       
    }
  else
    {
      t=1;  disc=-d;  n=(d+1)/4; 
    }
  rootdisc=sqrt(-disc);
  switch (d) {
    case 1:  name='i'; nunits=4; fundunit=Quad(0,1); break;
    case 2:  name='t'; nunits=2; fundunit=Quad(-1);  break;
    case 3:  name='w'; nunits=6; fundunit=Quad(0,1); break;
   default:  name='a'; nunits=2; fundunit=Quad(-1); 
  }
 switch (d) {
 case 1: case 2: case 3: case 7: case 11:       // Euclidean
//   h=0; quadgcd=&quadgcd1; quadbezout=&quadbezout1;  
   h=0; quadgcd=&quadgcd2; quadbezout=&quadbezout2;  
   break;
 case 19: case 43: case 67: case 163:           // Non-Euclidean classno. 1
   h=1; quadgcd=&quadgcd2; quadbezout=&quadbezout2;
   break;
 case   5: case   6: case  10: case  13: case  15: case  22: // class number 2
 case  35: case  37: case  51: case  58: case  91: case 115:
 case 123: case 187: case 235: case 267: case 403: case 427:
   h=2;  quadgcd=&quadgcd2; quadbezout=&quadbezout2;
   break;
 default: 
   cerr << "Field not imag quadratic of class number 1 or 2!  Aborting."<<endl;
   exit(1);    //abort program instantly!
 }
  Quadunits::init();
}

void Field::display(ostream& s)
{
  s<<"Field Q(sqrt("<<-d<<"))\tdiscriminant = "<<disc;
  s<<"\tmin poly("<<name<<") = "<<name<<"^2"; if(t) s<<"-"<<name; 
  s<<"+"<<n<<".\n";
  switch (h) {
  case 0:                                     // Euclidean
    s << "Euclidean" << endl;
    break;
  case 1:                                     // Non-Euclidean h=1
    s << "Non-Euclidean Class Number One" << endl;
    break;
  case 2:
    s << "Class Number Two" << endl;
    break;
  default: 
    s << "Class number > 2" << endl;
  }
  s << endl;
}

// END OF FILE field.cc
