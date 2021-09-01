load "core.m";
load "input.m";
load "levels.m";
load "homology.m"; // loads "input.m"
load "hecke.m";    // loads "homology.m", "region.m"
load "bianchi.m";

P2a,P2b,P3a,P3b:=Explode(PrimeIdeals(3));

N:=P2a*P3b;
assert N eq ideal<OK|w>;
assert  DimH(N,false) eq 1;
data := DimH(N,true);
new  := NewS(N,data);
Eigen(N,data,new,PrimeIdeals(10)); // returns
// [* 1, 1, 0, -1*] -1
// Eigen(N,data,new,PrimeIdeals(50)); // takes ages and returns
// [* 1, 1, 0, -1, -2, 2, -4, -2, 6, 6, 0, -4, -2, 2, 8, 0, -6*] -1

assert (1-2*w)^2 eq -23;
P23:= ideal<OK|1-2*w> ;
N:=P23;

N:=5*P3a;
assert  DimH(N,false) eq 2;
time data := DimH(N,true);
time new  := NewS(N,data);
assert Dimension(new) eq 2;
// Eigen(N,data,new,PrimeIdeals(10)); // fails when dim>1

time Ts:=[HeckeMatrix(P,N,data) : P in PrimeIdeals(30)]; // takes a while
// They commute:
assert &and[A*B eq B*A : A,B in Ts];
[#Eigenvalues(T) : T in Ts]; // returns
// [ 0, 1, 1, 0 ]   so they do not split over Q

D:= Diagonalization(Ts);
eigs:=[[A[i,i] : A in D] : i in [1..2]];
L<a>:=Parent(eigs[1][1]);
assert eigs eq [[a,2,-1,-a],[1-a,2,-1,a-1]];

/*
for N in Levels(10) do
    fl:=IsPrincipal(N);
    d := DimH(N,false);
    if d gt 0 then
	print N,Norm(N),d;
    end if;
end for;
*/
