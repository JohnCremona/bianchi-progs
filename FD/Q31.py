from sage.all import QQ, NumberField, polygen

dk = -31
x = polygen(QQ)
k = NumberField(x**2-x+(1-dk)//4, 'w')
assert k.discriminant() == dk
w = k.gen()
wbar = 1-w
zero = k(0)
one = k(1)
Ok = k.ring_of_integers()
emb = next(e for e in k.embeddings(CC) if e(w).imag()>0)
rootd=RR(-dk).sqrt()

J = Matrix(2,2,[-k(1), k(0), k(0), k(1)])
Smat = Matrix(2,2,[k(0), k(-1), k(1), k(0)])
inf = NFCusp(k,infinity)
tri0 = [NFCusp(k,0), inf, NFCusp(k,1)]

n_alphas = 0
M_alphas = []
all_alphas = []

add_alpha(0,-1,1,0)
add_two_alphas(4, 1+2*w, +1)
add_four_alphas(3, w, 1-w)
add_two_alphas(3, 1+w, -1)

add_two_alphas(w, 3, +1);       # N(s)=8
add_two_alphas(1-w, 3, +1);
add_two_alphas(1+w, 3);         # N(s)=10
add_two_alphas(2-w, 3);
add_two_alphas(3+w, w-6, +1);   # N(s)=20
add_two_alphas(4-w, 5+w, +1);

alpha_denoms = list(Set(M[1,0] for M in M_alphas))
alphas = dict([(d,[alpha for M,alpha in zip(M_alphas, all_alphas) if M[1,0]==d]) for d in alpha_denoms])

alpha_inv = [all_alphas.index(apply(M,cusp(oo))) for M in M_alphas]
assert all([alpha_inv[alpha_inv[i]]==i for i in range(n_alphas)])
for i in range(n_alphas):
     r1,s = frac(all_alphas[i])
     r2,s2 = frac(all_alphas[alpha_inv[i]])
     assert s==s2
     assert (r1*r2+1)/s in Ok

def make_edge_matrix(a,b):
    return Matrix(2,2,[a.numerator(),b.numerator(),a.denominator(),b.denominator()])

pairs = [(i,j) for i,a in enumerate(all_alphas) for j,b in enumerate(all_alphas) if j>i and posdet(make_edge_matrix(a,b)) in alpha_denoms]

sigmas = [cusp(oo), cusp(w/2), cusp((1-w)/2)]

Tlist, triangles = make_triangles()
print("{} triangles (before eliminating congruences)".format(len(triangles)))
assert len(triangles)==44
T0,t0=reduce_triangles(Tlist, triangles)
print("{} triangles (after eliminating congruences)".format(len(t0)))
assert len(t0)==10

# cyclic triangles for alpha with M_alpha of order 3, i.e. alpha=r/s with r^3=+-1 (mod s) so alpha' = (r+1)/s or (r-1)/s.
# Equivalently, trace(M_alpha) = +/-1
cyclic_triangles = [i for i,M in enumerate(M_alphas) if i>2 and M.trace() in [1,-1]]
assert len(cyclic_triangles)==0
print("no cyclic triangles")

# Lingham's vertices: C4, C5, C8, C9, C13, C14, C15, C16, C18, C19, C20, C23, and singular C25, C26
C4 = cusp(to_k(all_alphas[19])-1)
C5 = cusp(to_k(all_alphas[17])+1)
C8 = all_alphas[14]
C9 = all_alphas[15]
C10 = cusp(to_k(all_alphas[14])+1)
C13 = all_alphas[10]
C14 = all_alphas[11]
C15 = cusp(to_k(all_alphas[1])-1)
C16 = cusp(to_k(all_alphas[2])+w)
C18 = cusp(to_k(all_alphas[7])-1)
C19 = all_alphas[6]
C20 = all_alphas[3]
C21 = all_alphas[7]
C23 = all_alphas[0]
C24 = cusp(1)
C25 = cusp((w-1)/2) # = -sigma2
C26 = sigmas[1]

# ML's aaa-triangles (with no singular vertex):
MLT = []
# tetrahedron 4.2(a)
MLT.append([inf, C5, C14])
MLT.append([inf, C5, C21])
MLT.append([inf, C14, C21])
# tetrahedron 4.2(b)
MLT.append([inf, C9, C20])
MLT.append([inf, C9,  C23])
MLT.append([inf, C20, C23])
# tetrahedron 4.2(c)
MLT.append([inf, C9, C14])
MLT.append([inf, C14, C20])
# tetrahedron 4.2(d)
MLT.append([inf, C13, C16])
# tetrahedron 4.2(e)
MLT.append([inf, C14, C16])
# polyhedron 4.2(f)
MLT.append([inf, C9, C10])
MLT.append([inf, C9, C21])
MLT.append([inf, C10, C21])
MLT.append([inf, C10, C24])
MLT.append([inf, C23, C24])
MLT.append([C9, C10, C21])

# check that we have all these:
missing_triangles = [t for t in MLT if not check_poly_in_list(t,[tri0]+t0)]
if missing_triangles:
    print(len(missing_triangles), " Triangles in ML's list not in ours:")
    print(missing_triangles)

# ML's aas-triangles with one singular vertex:
MLST = []
MLST.append([inf, C13, C25])
MLST.append([inf, C16, C25])
MLST.append([inf, C14, C26])
MLST.append([inf, C16, C26])

print("ML has {} singular triangles".format(len(MLST)))

