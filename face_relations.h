// FILE FACE_RELATIONS.H: Declaration of class face_relations

#if     !defined(_FACE_RELATIONS_H)
#define _FACE_RELATIONS_H      1       //flags that this file has been included

#include "edge_relations.h"

#define USE_SMATS

class face_relations {
public:
  face_relations() {;}
  face_relations(edge_relations*, scalar mod, int plus, int verb=0, scalar ch=scalar(0));
  int ncoords() const {return coord.nrows();}
  vec coords(int i) const {return coord.row(i);}
  mat get_coord() const {return coord;}
  scalar get_denom() const {return denom;}
  int get_rank() const {return rk;}
  scalar get_hmod() const {return hmod;}
  int gen(int i) const {return pivs[i];}

private:
  edge_relations* ER; // provides coord(i) and symbdata for nsymb, symbol(i), symbops
  int plusflag, verbose, ngens, nsymb, numrel, maxnumrel, rk;
  scalar characteristic, modulus, hmod, denom;
  vec_i pivs;
  action act_with(const mat22& M) {return ER->act_with(M);}
  action act_with(const Quad& a, const Quad& b, const Quad& c, const Quad& d) {return ER->act_with(a,b,c,d);}

  void make_relations();        // creates relation matrix relmat
  void solve_relations();       // computes kernel of relmat and sets rk, denom, coord[, freegens]

  mat22 M_alpha(int j) const {return Quad::SD.Mlist[j];}
  void add_face_rel(const vector<int>& rel, const vector<int>& types);
  void add_face_rel(const vector<int>& rel, const vector<int>& types, const vector<int>& signs);

  void triangle_relation_0();           // triangle relation for all fields
  void triangle_relation_1_3();         // extra triangle relation for fields 1, 3
  // these special relations for class number fields are now included in the general T,Q,H lists
  // void triangle_relation_2();           // extra triangle relation(s) for fields 19+
  // void square_relation_2();             // extra square relation for field 2
  // void rectangle_relation_7();          // extra rectangle relation for field 7
  // void hexagon_relation_11();           // extra hexagon relation for field 11

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
