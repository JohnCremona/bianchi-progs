Qx<x> := PolynomialRing(Rationals());
K<w> := NumberField(x^2+6);
assert Discriminant(K) eq -24;
OK := Integers(K);
w:= OK.2;
