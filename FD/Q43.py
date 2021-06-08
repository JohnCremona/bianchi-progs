from sage.all import QQ, NumberField, polygen, NFCusp, infinity, CC, RR, Matrix

dk = -43
x = polygen(QQ)
k = NumberField(x**2-x+(1-dk)//4, 'w')
assert k.discriminant() == dk
w = k.gen()
zero = k(0)
one = k(1)
Ok = k.ring_of_integers()
emb = next(e for e in k.embeddings(CC) if e(w).imag()>0)
rootd=RR(-dk).sqrt()

alphas = {1: [NFCusp(k,a) for a in [0]],
          2: [NFCusp(k, a/2) for a in [w, 1+w]],
          3: [NFCusp(k, a/3) for a in [w,-w,1-w,w-1,1+w,-(1+w)]]}

M_alphas = [Matrix(2,2,[0,-1,1,0]),
            Matrix(2,2,[-1+w,5, 2,-w]),
            Matrix(2,2,[w,5,2,1-w]),
            Matrix(2,2,[1-w,-4,3,-w]),
            Matrix(2,2,[w-1,-4,3,w]),
            Matrix(2,2,[w,-4,3,w-1]),
            Matrix(2,2,[-w,-4,3,1-w]),
            Matrix(2,2,[1+w,3-w,3,-1-w]),
            Matrix(2,2,[-1-w,3-w,3,1+w])]

all_alphas = [NFCusp(k,-M[1,1]/M[1,0]) for M in M_alphas]
n_alphas = len(all_alphas)

# triangles from EW:

P1T1 = [0,3*w/11,3*(w+1)/13]
P1T2 = [3*w/11,w/3,(5*w+2)/17]
P1T3 = [(w+1)/4, (5*w+2)/17, 4*(w+1)/13]

P2T1 = [0,3*w/11,3*(w-1)/11]
P2T2 = T3 = [infinity, (w-1)/3, w/3]

P3T1 = [(w+1)/3, (5*w+3)/13, 5*(w+1)/13]
P3T2 = T2 = [infinity, w/2, (w+1)/2]

P4T1 = T1 = [0,1,infinity]
P4T2 = [(w+1)/4, 3*(w+1)/13, (4*w+5)/17]
P4T3 = [(w+2)/4, (3*w+7)/13, (4*w+8)/17]
P4T4 = [(w+1)/3, 4*(w+1)/13, (4*w+5)/13]
