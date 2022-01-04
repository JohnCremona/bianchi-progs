import subprocess
    
def ideal_gen_coeffs(I):
    return " ".join([str(c) for c in list(I.gens_reduced()[0])])

def apdata(E, P):
    ap = E.reduction(P).trace_of_frobenius()
    return " ".join([ideal_gen_coeffs(P), str(ap)])

def check_modularity(E, primes, verbose=False):
    K = E.base_ring()
    field = K.discriminant().squarefree_part().abs()
    ab = ideal_gen_coeffs(E.conductor())
    np = len(primes)
    input_string = " ".join([str(field), ab, "1", str(np)] + [apdata(E,P) for P in primes])
    if verbose:
        print("input string: {}".format(input_string))
    cmd = "echo {} | /home/jec/bianchi-progs/modularity".format(input_string)
    if verbose:
        print("command line: {}".format(cmd))
    pipe = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
    if pipe.returncode:
        return None
    return pipe.stdout.readlines()[0].replace("\n","")
    
