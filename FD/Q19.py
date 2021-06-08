from sage.all import QQ, NumberField, polygen, NFCusp

dk = -19
x = polygen(QQ)
k = NumberField(x**2-x+(1-dk)//4, 'w')
assert k.discriminant() == dk
w = k.gen()
Ok = k.ring_of_integers()


alphas = {1: [NFCusp(k,0)],
          2: [NFCusp(k, a/2) for a in [w, 1+w]]}


all_alphas = sum(alphas.values(), [])
