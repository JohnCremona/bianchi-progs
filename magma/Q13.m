Qx<x> := PolynomialRing(Rationals());
K<w> := NumberField(x^2+13);
assert Discriminant(K) eq -4*13;
OK := Integers(K);
w:= OK.2;
