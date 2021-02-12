from sage.all import QQ, polygen, NumberField, cartesian_product_iterator, prod

def field(d):
    if not d in [1,2,3,7,11]:
        raise ValueError("only 1,2,3,7,11 are valid")
    t = 0 if d<3 else 1
    n = d if d<3 else (d+1)//4
    x = polygen(QQ)
    K = NumberField(x**2-t*x+n,'-ita---a---a'[d])
    K.n = n
    K.t = t
    return K, K.gen()

def ideal_from_HNF(K,H):
    a,c,d = H
    return K.ideal([a,c+d*K.gen()])

def ideal_HNF(I):
    H = I.pari_hnf().python()
    print(H)
    return [H[0][0],H[0][1],H[1][1]]

def old_ideal_label(I):
    H = ideal_HNF(I)
    N = H[0]*H[2]
    return "[%s,%s,%s]"%(N,H[1],H[2])

import re
ideal_label_regex = re.compile(r'\d+\.\d+\.\d+')

def parse_ideal_label(s):
    if not ideal_label_regex.match(s):
        raise ValueError("invalid label %s"%s)
    #s = re.sub(r']','',re.sub(r'\[','',s))
    return [int(i) for i in s.split(r'.')]

def ideal_from_label(K,s):
    H = parse_ideal_label(s)
    H[0]/= H[2]
    return ideal_from_HNF(K,H)

def ideal_divisor_iterator(I):
    IF = I.factor()
    for ei in cartesian_product_iterator([range(e+1) for p,e in IF]):
        yield prod([p**f for (p,e),f in zip(IF,ei)])

def ideal_divisor_iterator_multi(I):
    IF = I.factor()
    for ei in cartesian_product_iterator([range(e+1) for p,e in IF]):
        yield prod([p**f for (p,e),f in zip(IF,ei)]), prod([(1+e-f) for (p,e),f in zip(IF,ei)])

def read_dimtabeis(d, fname):
    K, a = field(d)
    with open(fname) as infile:
        for L in infile:
            if L.count('eight'):
                continue
            dd, w, label, dimall, dimcusp, dimeis = L.split()
            #H = parse_ideal_label(label)
            I = ideal_from_label(K,label)
            print("%s dimensions: %s %s %s"%(I,dimall,dimcusp,dimeis))

# NB The following function takes as input a dimeistab as produced by
# the C++ program dimeistab.cc, but will only work if this input file
# starts at level norm 1.

def make_dimtabnew(d, fname, min_norm=1, max_norm=None, verbose=False):
    K, a = field(d)
    outfname = fname+'.newdims'
    dimscuspnew={}
    with open(fname) as infile, open(outfname, 'w') as outfile:
        for L in infile:
            if L.count('Table'):
                outfile.write("#"+L)
                continue
            if L.count('Field'):
                #outfile.write(L[:-1]+'\tdim(cuspidal, new)\n')
                outfile.write('#Field\tWeight\tLevel\t\tall\tcusp\tcusp(new)\teis\n')
                continue
            dd, w, label, dimall, dimcusp, dimeis = L.split()
            assert (int(dd)==d) and int(w)==2
            dimcusp = int(dimcusp)
            I = ideal_from_label(K,label)
            Inorm = I.norm()
            if max_norm!=None and Inorm>max_norm:
                break
            dimcuspnew = dimcusp
            if dimcuspnew:
                for J,m in ideal_divisor_iterator_multi(I):
                    if dimcuspnew==0:
                        break
                    if m>1:  # else J=I
                        dimcuspnew -= m*dimscuspnew[J]
            dimscuspnew[I] = dimcuspnew

            if Inorm>=min_norm:
                if verbose:
                    print("%s dimensions: %s %s (%s new) %s"%(old_ideal_label(I),dimall,dimcusp,dimcuspnew,dimeis))
                outfile.write("%s \t%s \t%s \t%s \t%s \t%s \t%s\n"%(dd,w,label,dimall,dimcusp,dimcuspnew,dimeis))

def edit_haluk_output(f, fname):
    outfname = fname+'.new'
    with open(fname) as infile, open(outfname, 'w') as outfile:
        for L in infile:
            if L[0]!="[":
                continue
            L = L.replace(',',' , ')
            L = L.replace(']',' ]')
            lab = range(6)
            lab[0],lab[1],lab[2],lab[3],lab[4],lab[5], dimcusp, dimeis = L.split(None,7)
            dimcusp = int(dimcusp)
            dimeis = int(dimeis.split()[0])
            dimall = dimcusp+dimeis
            label = "".join(lab)
            print(label)
            outfile.write("%s\t2\t%s\t %s\t %s\t %s\n"%(f,label,dimall,dimcusp,dimeis))

