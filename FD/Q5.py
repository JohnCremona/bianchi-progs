from sage.all import QQ, NumberField, polygen

dk = -20
x = polygen(QQ)
if dk%4==1:
    k = NumberField(x**2-x+(1-dk)//4, 'w')
else:
    k = NumberField(x**2-dk//4, 'w')
assert k.discriminant() == dk
w = k.gen()
wbar = w.trace()-w
zero = k(0)
one = k(1)
Ok = k.ring_of_integers()
emb = next(e for e in k.embeddings(CC) if e(w).imag()>0)
rootd=RR(-dk).sqrt()
Ireps=[c.ideal() for c in k.class_group()]

J = Matrix(2,2,[-k(1), k(0), k(0), k(1)])
Smat = Matrix(2,2,[k(0), k(-1), k(1), k(0)])
inf = NFCusp(k,infinity)
tri0 = [NFCusp(k,0), inf, NFCusp(k,1)]

n_alphas = 0
M_alphas = []
alphas = []

add_alpha(0,-1,1,0)   # 0/1
add_alpha(w,2,2,-w)   # w/2
add_two_alphas(2*w, w-4, +1)  # (w-4)/(2*w) = (5+4w)/10 = 1/2+2w/5 (and negative)

alpha_denoms = list(Set(M[1,0] for M in M_alphas))
alpha_denom_dict= dict([(d,[alpha for M,alpha in zip(M_alphas, alphas) if M[1,0]==d]) for d in alpha_denoms])

alpha_inv = [alphas.index(cusp(oo).apply(M.list())) for M in M_alphas]
assert all([alpha_inv[alpha_inv[i]]==i for i in range(n_alphas)])
for i in range(n_alphas):
     r1,s = frac(alphas[i])
     r2,s2 = frac(alphas[alpha_inv[i]])
     assert s==s2
     assert (r1*r2+1)/s in Ok

B = Matrix(2,2,[1+w,2,2,1-w])
assert B.det()==2
assert B*B==2*M_alphas[1]

def make_edge_matrix(a,b):
    return Matrix(2,2,[a.numerator(),b.numerator(),a.denominator(),b.denominator()])

pairs = [(i,j) for i,a in enumerate(alphas) for j,b in enumerate(alphas) if j>i and posdet(make_edge_matrix(a,b)) in alpha_denoms]

sigmas = [cusp(oo), cusp((1+w)/2)]

Tlist, triangles = make_triangles()
print("{} triangles (before eliminating congruences)".format(len(triangles)))
assert len(triangles)==0
T0,t0=reduce_triangles(Tlist, triangles)
print("{} triangles (after eliminating congruences)".format(len(t0)))
assert len(t0)==0

# cyclic triangles for alpha with M_alpha of order 3, i.e. alpha=r/s with r^3=+-1 (mod s) so alpha' = (r+1)/s or (r-1)/s.
# Equivalently, trace(M_alpha) = +/-1
cyclic_triangles = [i for i,M in enumerate(M_alphas) if i>2 and M.trace() in [1,-1]]
assert len(cyclic_triangles)==0
print("no cyclic triangles")

# Triangles (aas so not yet found automatically)
JCT = [tri0]
JCT.append([sigmas[0], alphas[1], sigmas[1]])
JCT.append([sigmas[0], alphas[2], cusp((1-w)/2)])

# Square
JCS = []
JCS.append([sigmas[0], cusp((w-1)/2), alphas[1], sigmas[1]])

# This square has a 2-fold symmetry:
assert apply(M_alphas[1],cycle_poly(JCS[0],2)) == JCS[0]
# and is equivalent to its image under the normaliser matrix B:
assert apply(B, JCS[0]) == cycle_poly(JCS[0], -1)

# Polygons from JB's thesis

# Triangles
JBT = []
JBT.append(tri0)
JBT.append([sigmas[1], cusp((1+w)/3), cusp((2+w)/3)])

# Squares
JBS = []
JBS.append([alphas[0], sigmas[0], sigmas[1], cusp((1+w)/3)])
JBS.append([alphas[0], cusp((1+w)/3), cusp(2*w/5), cusp((w-1)/3)])
JBS.append([sigmas[0], sigmas[1], alphas[1], cusp((w-1)/2)])

assert poly_equiv(apply(B,cycle_poly(JBS[0],2)), JBS[0])

# The 2nd and 3rd of JB's squares are both equivalent to the original square:
assert JBS[1] == apply(Smat,JCS[0])
assert JBS[2] == cycle_poly(reverse_poly(JCS[0]),-1)

# Yasaki Data

# Vertices (indexed from 1)
DYV = [None, cusp(-(1+w)/2), cusp(oo), cusp((1+w)/3), cusp(2/w), cusp(0), cusp(-w/2), cusp(-1), cusp(-(2+w)/3)]

# Triangles
DYT = []
DYT.append([DYV[i] for i in [3,4,5]])
DYT.append([DYV[i] for i in [1,2,6]])
DYT.append([DYV[i] for i in [1,3,8]])
DYT.append([DYV[i] for i in [2,5,7]])

# check that we have all these:
assert poly_equiv(DYT[3], JCT[0])

# Squares
DYT = []
DYT.append([DYV[i] for i in [1,3,5,2]])
DYT.append([DYV[i] for i in [2,5,4,6]])
DYT.append([DYV[i] for i in [6,4,3,1]])
DYT.append([DYV[i] for i in [2,1,8,7]])
DYT.append([DYV[i] for i in [7,8,3,5]])



# missing_triangles = [t for t in DYT if not check_poly_in_list(t,[tri0]+JCT)]
# if missing_triangles:
#     print("Triangles in DY's list not in ours:")
#     print(missing_triangles)

