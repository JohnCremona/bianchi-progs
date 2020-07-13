import subprocess

BIANCHI_BIN_DIR = '/home/jec/bianchi-progs/'

def ideal_gen_coeffs(I):
    return " ".join([str(c) for c in list(I.gens_reduced()[0])])

def apdata(E, P):
    ap = E.reduction(P).trace_of_frobenius()
    return " ".join([ideal_gen_coeffs(P), str(ap)])

def check_modularity_modp(E, primes, p=3, verbose=False):
    K = E.base_ring()
    field = K.discriminant().squarefree_part().abs()
    ab = ideal_gen_coeffs(E.conductor())
    np = len(primes)
    input_string = " ".join([str(field), ab, str(p), "1", str(np)] + [apdata(E,P) for P in primes])
    if verbose:
        print("input string: {}".format(input_string))
    cmd = "echo {} | ./modularity_modp".format(input_string)
    if verbose:
        print("command line: {}".format(cmd))
    pipe = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                            shell=True, text=True, cwd=BIANCHI_BIN_DIR)
    if pipe.returncode:
        return None
    outputlines = [str(L) for L in pipe.stdout.readlines()]
    if outputlines:
        return outputlines[0].replace("\n","")
    else:
        return False

# example:

def test1(p=3, verbose=False):
    x = polygen(QQ)
    K = NumberField(x**2+1, 'i')
    i = K.gen()
    E = EllipticCurve([i + 1, i - 1, i + 1, -5*i, 2*i])
    NE = E.conductor()
    primes = [P for P in K.primes_of_bounded_norm(100) if NE.valuation(P)==0]
    res = check_modularity_modp(E,primes,p, verbose)
    if res:
        print("{} matches Bianchi modular form(s) {} (mod {})".format(E, res, p))
    else:
        print("No Bianchi modular form found which matches {} (mod {})".format(E, p))
