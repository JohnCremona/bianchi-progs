// FILE FACE_RELATIONS.H: Declaration of class face_relations

#if     !defined(_FACE_RELATIONS_H)
#define _FACE_RELATIONS_H      1       //flags that this file has been included

#include <eclib/method.h>
#include <eclib/matrix.h>
#include "edge_relations.h"

class face_relations {
public:
  face_relations() {;}
  face_relations(edge_relations*, int plus, int verb=0);
  int ncoords() const {return coord.nrows();}
  vec coords(int i) const {return coord.row(i);}
  long get_denom() const {return denom;}
  long get_rank() const {return rk;}
  long get_hmod() const {return hmod;}
  long gen(int i) const {return pivs[i];}

private:
  edge_relations* ER; // provides coord(i) and symbdata for nsymb, symbol(i), symbops
  int plusflag, verbose;
  long ngens, nsymb, numrel, maxnumrel;
  long hmod, denom, rk;
  vec pivs;

  void make_relations();        // creates relation matrix relmat
  void solve_relations();       // computes kernel of relmat and sets rk, denom, coord[, freegens]

  void add_face_rel(const vector<int>& rel, const vector<int>& types);
  void triangle_relation_0();           // triangle relation for all fields
  void triangle_relation_1_3();         // extra triangle relation for fields 1, 3
  void triangle_relation_2();           // extra triangle relation(s) for fields 19+
  void cyclic_triangle_relation(int i); // generic cyclic triangle relation
  void square_relation_2();             // extra square relation for field 2
  void rectangle_relation_7();          // extra rectangle relation for field 7
  void hexagon_relation_11();           // extra hexagon relation for field 11

  void general_triangle_relation(const vector<int>& tri);                         // generic triangle relation
  void general_square_relation(const vector<int>& squ, const vector<Quad>& xyz);  // generic square relation

#ifdef USE_SMATS
  smat relmat;
#else
  mat relmat;
#endif
protected:
  mat coord;
  friend class newform;
  friend class newforms;
};


#endif