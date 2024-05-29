// FILE FACE_RELATIONS.H: Declaration of class face_relations

#if     !defined(_FACE_RELATIONS_H)
#define _FACE_RELATIONS_H      1       //flags that this file has been included

#include "edge_relations.h"

#define USE_SMATS

class face_relations {
public:
  face_relations() {;}
  face_relations(edge_relations*, int plus, int verb=0, long ch=0);
  long ncoords() const {return coord.nrows();}
  vec coords(int i) const {return coord.row(i);}
  long get_denom() const {return denom;}
  long get_rank() const {return rk;}
  long get_hmod() const {return hmod;}
  long gen(int i) const {return pivs[i];}

private:
  edge_relations* ER; // provides coord(i) and symbdata for nsymb, symbol(i), symbops
  P1N* P1;  // shortcut to ER->P1
  int plusflag, verbose;
  long ngens, nsymb, numrel, maxnumrel;
  long hmod, characteristic, denom, rk;
  vec pivs;

  void make_relations();        // creates relation matrix relmat
  void solve_relations();       // computes kernel of relmat and sets rk, denom, coord[, freegens]

  void add_face_rel(const vector<long>& rel, const vector<int>& types);
  void add_face_rel(const vector<long>& rel, const vector<int>& types, const vector<int>& signs);

  void triangle_relation_0();           // triangle relation for all fields
  void triangle_relation_1_3();         // extra triangle relation for fields 1, 3
  void triangle_relation_2();           // extra triangle relation(s) for fields 19+
  void square_relation_2();             // extra square relation for field 2
  void rectangle_relation_7();          // extra rectangle relation for field 7
  void hexagon_relation_11();           // extra hexagon relation for field 11

  void aaa_triangle_relation(const POLYGON&, int check=1);   // generic aaa-triangle relation
  void aas_triangle_relation(const POLYGON&, int check=1);   // generic aas-triangle relation
  void general_square_relation(const POLYGON&, int check=1);  // generic square relation
  void general_hexagon_relation(const POLYGON&, int check=1);  // generic hexagon relation
  void general_relation(const vector<action>& Mops, const vector<int>& types, const vector<int>& signs,
                        int symmetry=0, int check=1);  // generic relation

#ifdef USE_SMATS
  smat relmat;
  vector<int> relmat_rowdata;
  int maxrowsize;
#else
  mat relmat;
#endif
protected:
  mat coord;
  friend class newform;
  friend class newforms;
};


#endif
