Qx<x> := PolynomialRing(Rationals());
K<w> := NumberField(x^2-x+8);
assert Discriminant(K) eq -31;
OK := Integers(K);
w:= OK.2;
