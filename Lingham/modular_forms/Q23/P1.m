Qx<x> := PolynomialRing(Rationals());
K<w> := NumberField(x^2-x+6);
OK := Integers(K);
N := w*OK;
R := quo<OK|N>;
P1 := ProjectiveSpace(R,1);
P1![2,1]; // OK
P1![1,2]; // Runtime error in '!': Projective points may not have all coordinates zero
