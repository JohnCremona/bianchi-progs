from sage.all import QQ, NumberField, polygen

dk = -24
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
add_alpha(1+w,2-w,2,-1-w)   # (1+w)/2
add_two_alphas(2*w,5,+1)

alpha_denoms = list(Set(M[1,0] for M in M_alphas))
alpha_denom_dict= dict([(d,[alpha for M,alpha in zip(M_alphas, alphas) if M[1,0]==d]) for d in alpha_denoms])

alpha_inv = [alphas.index(cusp(oo).apply(M.list())) for M in M_alphas]
assert all([alpha_inv[alpha_inv[i]]==i for i in range(n_alphas)])
for i in range(n_alphas):
     r1,s = frac(alphas[i])
     r2,s2 = frac(alphas[alpha_inv[i]])
     assert s==s2
     assert (r1*r2+1)/s in Ok

def make_edge_matrix(a,b):
    return Matrix(2,2,[a.numerator(),b.numerator(),a.denominator(),b.denominator()])

pairs = [(i,j) for i,a in enumerate(alphas) for j,b in enumerate(alphas) if j>i and posdet(make_edge_matrix(a,b)) in alpha_denoms]

sigmas = [cusp(oo), cusp(w/2)]

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

# Square
JCS = []

