Format of Newspace and Newform data files
=========================================

In directory newspaces/<field-label>/.  Let n2r = 2-rank of class group.

Newspace data file for each level
---------------------------------

- in file newspaces/<field-label>/<level-label>:

- Up to 4 lines if n2r>0 or 2 lines if n2r==0 (odd class number). If
  n=0 on line 1, just 1 line, else 2 or 4 lines.

- field-label level-label n (number of homological newforms)
- list of n dimensions (degrees of principal Hecke fields)
- (if n2r>0) list of n full dimensions (degrees of full Hecke fields)
- (if n2r>0) list of n trivial-character flags

Newform data file for each homological newform
----------------------------------------------

- in file newspaces/<field-label>/<level-label>-<id>

- Only for newforms with trivial character, except for class group C4
  fields where both genus characters are encoded.

- See below for the format of individual data items.

- field-label level-label id
- principal dimension; if n2r>0, full dimension and list of n2r characters
  (these are all 1 if trivial char, and either 1 or -1 if C4).
- Principal Hecke field k (degree d).
- If n2r>0:
  - Full Hecke field K via k mod squares data of rank r.
  - Full Hecke field K as absolute field (degree d*2^r).
  - Embedding of k into K: rational nxn matrix, n=2^r*d.
  - Images in K of the r relative generators (sqrts of the r
    elements in the field-mod-squares data).
  - D (self-twist discriminant dividing field discriminant, or 0).
  - ntw (number of nontrivial unramified quadratic twists up
    to Galois conjugation, at most 2^n2r-1).
  - list of ntw negative discriminants dividing field discriminant.
- base-change code (-1,0,1,2)

- (trivial char only) Atkin-Lehner eigenvalues, one per bad prime:
  - prime-label, +/-1

- aP coefficients ( = T(P) eigenvalues for good P), variable number:
  - prime-label, eigenvalue

Format for individual data items
--------------------------------

Format for a number field: _either_ Q _or_ var coeffs
 e.g. "a [1 -3 -1 1]" for k=Q(a), a^3-3a^2-a-1=0

Format for field elements: d+1 integers, first d coefficients of
numerator w.r.t. power basis, then denominator.

Format for field-mod-squares data defining K as a polyquadratic
extension of k: r (rank, >=0) then r elements of k (nonzero and
multiplicatively linearly independent modulo squares).

Format for a rational nxn matrix: n^2+1 integers, the entries (by
rows) and common denominator.

Format for eigenvalues:
 - (1) field-element
 - (2) root_index (0,+i<r^r) if r>0 (full field != principal field)
 - (3) xf (0,1,-1) if r>0 and complex, where 'complex' is iff the first
   field-mod-squares generator is -1 ( only relevant for C4 data).
