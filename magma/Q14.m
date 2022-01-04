Qx<x> := PolynomialRing(Rationals());
K<w> := NumberField(x^2+14);
assert Discriminant(K) eq -4*14;
OK := Integers(K);
w:= OK.2;
