load "all.m";
load "sorting.m";
P2a,P2b,P3a,P3b:=Explode(PrimeIdeals(3));

N:=P2a*P3b;
KI:=quo<OK|N>;
PKI:=ProjectiveSpace(KI,1);
MKI:=RMatrixSpace(KI,2,2);
PPN := [s[1]^s[2]:s in Factorization(N)];
List(PPN, N, KI, PKI);
