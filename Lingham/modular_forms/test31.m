D := 31;
load "ncore.m";
load "swan.m";

time S1:=InitialList(16);
time S2:=Check(S1);

// D=31 gives (?)

S2:={<-1,1>,<0,1>,<1,1>,<3,1-w>,<-3,w>,<w-2,3>,<w+1,3>,<2*w-3,4>,<2*w-1,4>,<2*w+1,4>,
     <2*w+1,4-w>,<2*w-3,w+3>,<w-6,w+3>,<w+5,4-w>,<5,1-w>,<-5,w>,<2*w-5,w+3>,<2*w+3,4-w>,
     <-3,w+1>,<w,3>,<3,2-w>,<w-1,3>}; 

S2c:=[<-1,1>,<0,1>,<1,1>,<-3,w>,<3,1-w>,<w-2,3>,<w+1,3>,<2*w-3,4>,<2*w-1,4>,<2*w+1,4>,
      <w-6,w+3>,<2*w+1,4-w>,<2*w-3,w+3>,<w+5,4-w>,<-5,w>,<5,1-w>,<2*w+3,4-w>,<2*w-5,w+3>,
      <w,2>,<w-1,2>];

S2a:=[Cusps![x[1],x[2]]:x in S2c];

S2b:=[x[1]/x[2]:x in S2a];

time S3,I1,SP1:=Inter(S2);
assert S2 eq S3;
SP2:=[x:x in SP1];
I2:=Boundary(I1,S3);               // Vertices of f.r.
S4:=[x:x in S3];                   // List of hemispheres forming the f.r.

S5:=[x:x in S4 | translation(Cusps![x[1],x[2]])[2] eq MI]; 

//S5b:=[x:x in S5a | x notin [<2*w-3,4>,<3,1-w>]];

NC:=[Cusps![s[1],s[2]]:s in S5];

// Inversion matrix corresponding to each hemisphere   
IM:=[<c,InversionMatrix(c,NC)>: c in S4];

S6:=[x[1]/x[2]:x in S5] cat SP2;         // List of generating edges
//S6:=[x[1]:x in EndPoints | x[2] ne 0];

VO:=Vorbits(I2);                         // Orbits of vertices under Gamma
RP:=RepPoly(VO,S2c);                     // Representative polyhedra


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

S6:=[x:x in EndPoints | x[2] ne 0];

function Rels(RepPoly)

for p in RepPoly do
<p[1],p[2],Infinity()>;
i,_,M1:=Explode(Image(S2a[p[1]],S2a[p[2]],NC));
a:=Position(S6,i);
i,_,M2:=Explode(Image(S2a[p[2]],Cusps![1,0],NC));
b:=Position(S6,i);
i,_,M3:=Explode(ImageInf(S2a[p[1]],NC));
c:=Position(S6,i);
printf "%o_%o\n",M1,a;
printf "%o_%o\n",M2,b;
printf "%o_%o\n",M3,c;
assert ActC(M1,S6[a]) eq Cusps![S2a[p[1]][1],S2a[p[1]][2]]; 
assert actinfinity(M1)[1] eq S2a[p[2]][1]/S2a[p[2]][2]; 
<ActC(M1,EndPoints[a]),actinfinity(M1)>;
<ActC(M2,EndPoints[b]),actinfinity(M2)>;
<ActC(M3,EndPoints[c]),actinfinity(M3)>;
"";
<p[2],p[3],Infinity()>;
i,_,M1:=Explode(Image(S2a[p[2]],S2a[p[3]],NC));
a:=Position(S6,i);
i,_,M2:=Explode(Image(S2a[p[3]],Cusps![1,0],NC));
b:=Position(S6,i);
i,_,M3:=Explode(ImageInf(S2a[p[2]],NC));
c:=Position(S6,i);
printf "%o_%o\n",M1,a;
printf "%o_%o\n",M2,b;
printf "%o_%o\n",M3,c;
<ActC(M1,EndPoints[a]),actinfinity(M1)>;
<ActC(M2,EndPoints[b]),actinfinity(M2)>;
<ActC(M3,EndPoints[c]),actinfinity(M3)>;
"";
<p[3],p[1],Infinity()>;
i,_,M1:=Explode(Image(S2a[p[3]],S2a[p[1]],NC));
a:=Position(S6,i);
i,_,M2:=Explode(Image(S2a[p[1]],Cusps![1,0],NC));
b:=Position(S6,i);
i,_,M3:=Explode(ImageInf(S2a[p[3]],NC));
c:=Position(S6,i);
printf "%o_%o\n",M1,a;
printf "%o_%o\n",M2,b;
printf "%o_%o\n",M3,c;
<ActC(M1,EndPoints[a]),actinfinity(M1)>;
<ActC(M2,EndPoints[b]),actinfinity(M2)>;
<ActC(M3,EndPoints[c]),actinfinity(M3)>;
"";
"";
end for;
return 0;
end function;


for v in VO do
   vv:=v[1];
   vv;
   "";
   for x in v do
      B:=Boundary(Vertices(x),S3);
      M:=FromPtoQ(x,vv);
      for m in M do
         for b in B do
            act(m,b);
         end for;
         actinfinity(m);
         "";
      end for;
   end for;
   "";
   "";
end for;

//Rels(RP);





