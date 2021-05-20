import os
HOME = os.getenv("HOME")
UPLOAD_DIR = os.path.join(HOME, "bmf-upload")
DATA_DIR = os.path.join(HOME, "bianchi-data")
DIM_DIR = os.path.join(DATA_DIR, "dims")
FORM_DIR = os.path.join(DATA_DIR, "newforms")

from sage.all import QQ, polygen, NumberField, cartesian_product_iterator, prod
from psort import ideal_label

fields = [1,2,3,7,11,19]

def field(d):
    if not d in fields:
        raise ValueError("only d in {} are valid".format(fields))
    t = 0 if d<3 else 1
    n = d if d<3 else (d+1)//4
    x = polygen(QQ)
    gen_name = 'i' if d==1 else 't' if d==2 else 'a'
    K = NumberField(x**2-t*x+n, gen_name)
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

# NB The following function takes as input a dimtabeis as produced by
# the C++ program dimtabeis.cc, but will only work if this input file
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

def read_dimtabeis_new(d, fname):
    """Read a file dimtabeis*newdims as output by make_dimtabnew().
    Ignoring lines starting with #, lines have 7 columns:

    d  field  (1, 2 or abs(disc))
    w  weight (2)
    level_label (HNF-style, e.g. 1.0.1)
    dimall      (dimension of full gl2-space, =dimcusp+dimeis)
    dimcusp     (dimension of cuspidal gl2-space)
    dimcuspnew  (dimension of new cuspidal gl2-space)
    dimeis      (dimension of eisenstein gl2-space)

    NB 1. This script is intended for use with gl2, weight 2, data, so
    the records will not have the keys 'sl2_dims', 'sl2_new_totaldim'
    or 'sl2_cusp_totaldim'.  When uploading we'll need to make sure
    not to overwrite the sl2 data which exists in the database.

    NB 2. The LMFDB only stores cupidal dimensions, so the dimall and
    dimeis columns are ignored.

    """
    K, a = field(d)
    D = K.discriminant().abs()
    field_label = "2.0.{}.1".format(D)
    data = {}
    filename = os.join(DIM_DIR, fname)
    with open(filename) as infile:
        for L in infile:
            if L[0] == '#':
                continue
            dd, w, label, dimall, dimcusp, dimcuspnew, dimeis = L.split()
            assert int(dd) == d
            assert int(w) == 2
            N = ideal_from_label(K,label)
            level_label = ideal_label(N)
            label = "-".join(field_label, level_label)
            level_norm = N.norm()
            gl2_dims = {'2': {'new_dim': dimcuspnew, 'cuspidal_dim': dimcusp}}
            gl2_new_totaldim = dimcuspnew
            gl2_cusp_totaldim = dimcusp
            data['label'] = {
                'field_label': field_label,
                'field_absdisc': D,
                'level_norm': level_norm,
                'level_label': level_label,
                'label': label,
                'gl2_dims': gl2_dims,
                'gl2_new_totaldim': gl2_new_totaldim,
                'gl2_cusp_totaldim': gl2_cusp_totaldim,
            }
    return data

bmf_dims_schema = {
    'gl2_dims': 'jsonb',
    'sl2_dims': 'jsonb',
    'field_absdisc': 'integer',
    'label': 'text',
    'field_label': 'text',
    'level_norm': 'bigint',
    'level_label': 'text',
    'gl2_new_totaldim': 'integer',
    'gl2_cusp_totaldim': 'integer',
    'sl2_new_totaldim': 'integer',
    'sl2_cusp_totaldim': 'integer'
}

int_cols = ['field_absdisc', 'level_norm', 'gl2_cusp_totaldim', 'gl2_new_totaldim', 'sl2_cusp_totaldim', 'sl2_new_totaldim']
str_cols = ['label', 'field_label', 'level_label']
dict_cols = ['gl2_dims', 'sl2_dims']

def encode_col(colname, col=None):
    if 'sl2' in colname:
        return "\N"
    if colname in int_cols:
        return str(col)
    if colname in dict_cols:
        return str(col).replace("'", '"')
    return col

def one_bmf_dims_line(record):
    return "|".join([encode_col(col, record.get(col, None)) for col in bmf_dims_schema])

def write_bmf_dims_upload_file(data, fname):
    """data is a dict as returned by read_dimtabeis_new(), with keys level labels, each value a dict with keys

    'label', 'field_label', 'field_absdisc', 'level_norm',
    'level_label', 'gl2_dims'_dims, 'gl2_new_totaldim'_new_totaldim,
    'gl2_cusp_totaldim'_cusp_totaldim.

    Output columns for the LMFDB bmf_dims table
    """
    filename = os.join(UPLOAD_DIR, fname)
    with open(filename) as infile:

        # Write header lines

        keys = list(bmf_dims_schema.keys())
        keys.remove('label')
        keys = ['label'] + keys
        vals = [bmf_dims_schema[k] for k in keys]
        outfile.write("|".join(keys))
        outfile.write("\n")
        outfile.write("|".join(vals))
        outfile.write("\n\n")

        # Write data lines

        nlines = 0
        for label, record in data.items():
            infile.write(one_bmf_dims_line(record) + "\n")
            nlines +=1
            if nlines%1000 == 0:
                print("{} lines written so far to {}",format(nlines, filename))

    print("{} lines written to {}",format(nlines, filename))
