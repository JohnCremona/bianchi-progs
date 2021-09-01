alphabet := {@ "a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l",
     "m", "n", "o", "p", "q", "r", "s", "t", "u", "v", "w", "x", "y", "z" @};

// Endpoints of generating edges

EndPoints:=[Cusps![2*w-3,4],
   	    Cusps![2*w-1,4],
            Cusps![9*w-13,20],     // (2*w+1)/(4-w)
            Cusps![9*w+4,20],      // (2*w-3)/(w+3)
            Cusps![w-2,3],
            Cusps![3*w-3,8],       // -3/w
            Cusps![3*w,8],         // 3/(1-w)
            Cusps![w,3],
  	    Cusps![0,1],
	    Cusps![w-1,3],
   	    Cusps![3,2-w],         // (3*w+3)/10
   	    Cusps![-3,w+1],        // (3*w-6)/10
            Cusps![w-1,2],
   	    Cusps![w,2],
	    Cusps![1,0]  	    
];

np:=#EndPoints-1;

MI:=MK!Matrix([[1,0],[0,1]]);
M2:=MK!Matrix([[w-2,w+1],[3,2-w]]);
M3:=MK!Matrix([[1,1],[0,1]]);
M4:=MK!Matrix([[-3,w],[w-1,3]]);
M5:=MK!Matrix([[3-2*w,-2*w-5],[-w-3,w-6]]);
M6:=MK!Matrix([[-2*w+3,-8],[-w-3,2*w-3]]);
M7:=MK!Matrix([[0,1],[-1,0]]);
M8:=MK!Matrix([[-w-5,2*w-7],[w-4,2*w+1]]);
M9:=MK!Matrix([[-w-1,w-2],[-3,w+1]]);
M10:=MK!Matrix([[-1,0],[1,-1]]);
M11:=MK!Matrix([[-1,1],[-1,0]]);
M12:=MK!Matrix([[3,1-w],[-w,-3]]);
M13:=MK!Matrix([[2*w-3,2*w+6],[4,-2*w+3]]);
M14:=MK!Matrix([[2*w-1,8],[4,-2*w+1]]);
M15:=MK!Matrix([[-w,-3],[-3,w-1]]);
M16:=MK!Matrix([[3,-w-1],[2-w,-3]]);
M17:=MK!Matrix([[-3,w-2],[w+1,3]]);
M18:=MK!Matrix([[-1,w-1],[0,1]]);
M19:=MK!Matrix([[-1,w],[0,1]]);
M20:=MK!Matrix([[2*w-3,w+7],[4,1-2*w]]);
M21:=MK!Matrix([[2*w+1,8-w],[4,-2*w+1]]);
M22:=MK!Matrix([[-w-1,-3],[-3,w-2]]);
M23:=MK!Matrix([[1,0],[1,1]]);
M24:=MK!Matrix([[-2*w-1,-8],[w-4,2*w+1]]);
M25:=MK!Matrix([[0,1],[1,0]]);
M26:=MK!Matrix([[1,0],[0,-1]]);
M27:=MK!Matrix([[-3,w-2],[w-2,w+1]]);
M28:=MK!Matrix([[w-2,3],[w+1,2-w]]);
M29:=MK!Matrix([[w-2,w+1],[w+1,3]]);
M30:=MK!Matrix([[0,1],[-1,1]]);

rels:=[[<1,M22,6>,<1,M3,5>,<1,M4,7>],
       [<1,M5,3>,<1,MI,4>,<1,M22,5>],
       [<1,M4,2>,<1,MI,7>,<1,M6,4>],
    

       [<1,M15,6>,<1,MI,8>,<1,M7,9>],
       [<1,M27,5>,<1,MI,11>,<1,M7,9>],
       [<1,M16,7>,<1,MI,11>,<1,M15,10>],

       [<1,M4,7>,<1,M16,8>,<1,MI,11>],
       [<1,M4,7>,<1,M15,9>,<1,MI,8>],
       [<1,M16,7>,<1,MI,11>,<1,M15,10>],

       [<1,M22,6>,<1,M3,5>,<1,M4,7>],
       [<1,M16,8>,<1,MI,11>,<1,M4,7>],
       [<1,M16,9>,<1,MI,11>,<1,M22,5>],

    
       [<1,M6,4>,<1,M8,5>,<1,M3,3>],
       [<1,M8,3>,<1,M9,4>,<1,M3,5>],
       [<1,M5,3>,<1,MI,4>,<1,M22,5>],


       [<1,M27,5>,<1,MI,11>,<1,M7,9>],
       [<1,M28,9>,<1,M3,12>,<1,M16,11>],
       [<1,M23,5>,<1,M3,9>,<1,M29,12>],
       [<1,M30,9>,<1,MI,9>,<1,M11,9>],

       [<1,M29,9>,<1,M3,12>,<1,M22,5>],
       [<1,M16,9>,<1,MI,11>,<1,M22,5>],

      
       [<-1,M14,14>,<1,MI,14>,<1,M14,2>],
       [<1,M4,4>,<1,MI,7>,<1,M14,2>],
       [<1,M4,14>,<1,MI,7>,<-1,MI,14>],


       [<1,M14,13>,<1,MI,2>,<-1,MI,13>],
       [<1,M12,13>,<1,MI,6>,<-1,MI,13>],
       [<1,M12,3>,<1,MI,6>,<1,M14,2>],

       [<1,MI,1>,<1,M13,1>],
       [<1,MI,1>,<1,M20,2>],
       [<1,MI,2>,<1,M14,2>],
       [<1,MI,3>,<1,M24,3>],
       [<1,MI,4>,<1,M6,4>],
       [<1,MI,5>,<1,M2,5>],
       [<1,MI,6>,<1,M12,6>],
       [<1,MI,7>,<1,M4,7>],
       [<1,MI,8>,<1,M15,10>],
       [<1,MI,9>,<1,M7,9>], 
       [<1,MI,9>,<1,M25,9>], 
       [<1,MI,9>,<-1,M26,9>], 

       [<1,MI,11>,<1,M16,11>],
       [<1,MI,12>,<1,M17,12>],

       [<1,MI,13>,<-1,M18,13>],
       [<1,MI,14>,<-1,M19,14>]
];

function actinfinity(M)   // Gives the image of infinity under the matrix M    
   a:=M[1,1]; c:=M[2,1]; 
   if c ne 0 
      then return Cusps![a,c];
      else return Cusps![1,0];
   end if; 
end function;

for i in rels do
if #i eq 3 then
<ActC(i[1][2],EndPoints[i[1][3]]),actinfinity(i[1][2])>;
<ActC(i[2][2],EndPoints[i[2][3]]),actinfinity(i[2][2])>;
<ActC(i[3][2],EndPoints[i[3][3]]),actinfinity(i[3][2])>;
"";
end if;
if #i eq 4 then
<ActC(i[1][2],EndPoints[i[1][3]]),actinfinity(i[1][2])>;
<ActC(i[2][2],EndPoints[i[2][3]]),actinfinity(i[2][2])>;
<ActC(i[3][2],EndPoints[i[3][3]]),actinfinity(i[3][2])>;
<ActC(i[4][2],EndPoints[i[4][3]]),actinfinity(i[4][2])>;
"";
end if;
end for;

for i in rels do
if #i eq 2 then
<ActC(i[1][2],EndPoints[i[1][3]]),actinfinity(i[1][2])>;
<ActC(i[2][2],EndPoints[i[2][3]]),actinfinity(i[2][2])>;
"";
end if;
end for;


