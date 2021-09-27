alphabet := {@ "a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l",
     "m", "n", "o", "p", "q", "r", "s", "t", "u", "v", "w", "x", "y", "z" @};

// Centres of hemispheres used in the principal covering
C1:=Cusps![w+2,3-w];                         // (5*w-2)/12
C2:=Cusps![w-3,w+2];                         // (5*w-3)/12

C3:=Cusps![2*w-3,4];                         // (2*w-3)/4
C4:=Cusps![2*w-1,4];                         // (2*w-1)/4
C5:=Cusps![2*w+1,4];                         // (2*w+1)/4

C6:=Cusps![w-2,w+1];                         // (3*w+2)/8
C7:=Cusps![3,2-w];                           // (3*w+3)/8
C8:=Cusps![-3,w+1];                          // (3*w-6)/8
C9:=Cusps![w+1,2-w];                         // (3*w-5)/8

C10:=Cusps![w-2,3];                          // (w-2)/3
C11:=Cusps![w+1,3];                          // (w+1)/3

C12:=Cusps![-1,1];                           // -1
C13:=Cusps![0,1];                            // 0
C14:=Cusps![1,1];                            // 1

C15:=Cusps![w-1,2];                          // (w-1)/2
C16:=Cusps![w,2];                            // w/2   

C17:=Cusps![1,0];                            // infinity

EndPoints:=[C1,C2,C4,C6,C9,C11,C13,C15,C16,C17];
np:=#EndPoints-1;

// These are useful matrices in M(2,OK) .

MI:=MK!Matrix([[1,0],[0,1]]);
M2:=MK!Matrix([[0,1],[-1,0]]);
M3:=MK!Matrix([[-w+2,-w-3],[-w-1,w-4]]);
M4:=MK!Matrix([[-w+3,w+2],[-w-2,-w+3]]);
M5:=MK!Matrix([[w+3,-w+4],[-w+2,-w-1]]);
M6:=MK!Matrix([[w-5,2*w+1],[w+2,-w+3]]);
M7:=MK!Matrix([[w+1,3],[-w+2,-w-1]]);
M8:=MK!Matrix([[w+2,-w+3],[-w+3,-w-2]]);
M9:=MK!Matrix([[w-4,w+3],[w+1,-w+2]]);
M10:=MK!Matrix([[1,1],[0,1]]);
M11:=MK!Matrix([[0,1],[-1,1]]);
M12:=MK!Matrix([[-1,0],[-1,-1]]);
M13:=MK!Matrix([[1,0],[1,-1]]);
M14:=MK!Matrix([[w-3,w+2],[w+2,-w+3]]);
M15:=MK!Matrix([[2*w-1,6],[4,-2*w+1]]);
M16:=MK!Matrix([[w-2,3],[w+1,-w+2]]);
M17:=MK!Matrix([[w+1,-w+2],[3,-w-1]]);
M18:=MK!Matrix([[0,1],[1,0]]);
M19:=MK!Matrix([[1,0],[0,-1]]);
M20:=MK!Matrix([[-1,w-1],[0,1]]);
M21:=MK!Matrix([[-1,w],[0,1]]);
M22:=MK!Matrix([[-w+3,-2*w-1],[-w-2,w-5]]);
M23:=MK!Matrix([[-w-1,w-4],[w-2,w+3]]);

// The edge and face relations.  Each is a list of triples <coef,M,k> where:
//
// coef is an integer coefficient (+1 or -1);
// M is a matrix in M(2,OK); and
// k is an index in [1..9] denoting the edge type, as an index into the list EndPoints
//    defined in input.m as [C1,C2,C4,C6,C9,C11,C13,C15,C16,C17].
//    Here C17=oo so does not appear as a type of edge.

rels:=[[<1,M11,4>,<-1,MI,6>,<1,MI,7>],
       [<1,M12,5>,<-1,MI,6>,<1,M10,7>],
       [<1,M13,7>,<1,M10,7>,<-1,MI,7>],
    
       [<-1,MI,1>,<1,M3,3>,<1,MI,4>],
       [<-1,MI,1>,<1,M2,2>,<1,MI,7>],
       [<-1,MI,4>,<1,M2,5>,<1,MI,7>],
      
       [<-1,MI,1>,<1,MI,2>,<1,M4,7>],
       [<1,M2,1>,<-1,MI,2>,<1,MI,7>],

       [<-1,MI,2>,<1,M23,3>,<1,MI,5>],
       [<-1,MI,2>,<1,MI,8>,<-1,M22,8>],
       [<-1,MI,5>,<1,MI,8>,<-1,M7,8>],
      
       [<-1,MI,1>,<1,MI,9>,<-1,M8,9>],
       [<-1,MI,4>,<1,MI,9>,<-1,M3,9>],

       [<1,MI,1>,<1,M8,1>],
       [<1,MI,2>,<1,M14,2>],
       [<1,MI,3>,<1,M15,3>],
       [<1,MI,4>,<1,M16,4>],
       [<1,MI,5>,<1,M7,5>],
       [<1,MI,6>,<1,M17,6>],
       [<1,MI,7>,<1,M18,7>],
       [<1,MI,7>,<-1,M19,7>],
       [<1,MI,8>,<-1,M20,8>],
       [<1,MI,9>,<-1,M21,9>]
];

