from sage.all import QQ, NumberField, polygen, oo, NFCusp, CC, RR, Matrix, infinity, Set
from utils import (add_alpha, add_two_alphas, cusp, posdet, frac,
                   reduce_triangles, make_triangles, check_poly_in_list, apply)

dk = -23
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
add_two_alphas(3, 1+w, +1)
add_two_alphas(1+w, 2-w, +1)
add_two_alphas(2-w, 1+w, +1)
add_two_alphas(2+w, w-3, +1)
add_two_alphas(w-3, 2+w, +1)

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
assert len(triangles)==20
T0,t0=reduce_triangles(Tlist, triangles)
print("{} triangles (after eliminating congruences)".format(len(t0)))
assert len(t0)==5

# cyclic triangles for alpha with M_alpha of order 3, i.e. alpha=r/s with r^3=+-1 (mod s) so alpha' = (r+1)/s or (r-1)/s.
# Equivalently, trace(M_alpha) = +/-1
cyclic_triangles = [i for i,M in enumerate(M_alphas) if i>2 and M.trace() in [1,-1]]
assert len(cyclic_triangles)==0
print("no cyclic triangles")

# Lingham's vertices: C1, C2, C10, C13, C17, C19, C20, and singular C21, C22
C1 = all_alphas[12]
C2 = all_alphas[9]
C10 = all_alphas[6]
C13 = all_alphas[7]
C17 = all_alphas[3]
C19 = all_alphas[0]
C20 = cusp(1)
C21 = cusp((w-1)/2) # = -sigma2
C22 = sigmas[1]

# Lingham's triangles:
# with no singular vertex:
MLT = []
MLT.append([inf, C17,  C20])
MLT.append([inf, C20, C19])
MLT.append([inf, C19, C17])
MLT.append([C17, C19, C20])
MLT.append([inf, C1,  C19])
MLT.append([inf, C19, C10])
MLT.append([inf, C10, C1])
MLT.append([C1, C19, C10])
MLT.append([inf, C19, C2])
MLT.append([inf, C2, C1])
MLT.append([C1, C19, C2])
MLT.append([inf, C2, C13])

# check that we have all these:
missing_triangles = [t for t in MLT if not check_poly_in_list(t,[tri0]+t0)]
if missing_triangles:
    print("Triangles in ML's list not in ours:")
    print(missing_triangles)

# with a singular vertex:
MLST = []
MLST.append([inf, C2, C21])
MLST.append([inf, C13, C21])
MLST.append([C2, C13, C21])
MLST.append([inf, C1, C22])
MLST.append([inf, C10, C22])
MLST.append([C1, C10, C22])

print("ML has {} singular triangles".format(len(MLST)))

