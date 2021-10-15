import os
import re
from sage.all import QQ, polygen, NumberField, cartesian_product_iterator, prod
from psort import ideal_label, ideal_from_label

whitespace = re.compile(r'\s+')

HOME = os.getenv("HOME")
UPLOAD_DIR = os.path.join(HOME, "bmf-upload")
DATA_DIR = os.path.join(HOME, "bianchi-data")
DIM_DIR = os.path.join(DATA_DIR, "dims")
FORM_DIR = os.path.join(DATA_DIR, "newforms")


Qx = PolynomialRing(QQ,'x')

def split(line):
    return whitespace.split(line.strip())

def field(d):
    t,n = (1, (d+1)//4) if d%4==3 else (0, d)
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
new_ideal_label_regex = re.compile(r'\d+\.\d+')

def parse_ideal_label(s):
    if not ideal_label_regex.match(s):
        raise ValueError("invalid label %s"%s)
    #s = re.sub(r']','',re.sub(r'\[','',s))
    return [int(i) for i in s.split(r'.')]

def ideal_from_IQF_label(K,s):
    try:
        H = parse_ideal_label(s)
        H[0]/= H[2]
        return ideal_from_HNF(K,H)
    except ValueError:
        return ideal_from_label(K,s)

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
            I = ideal_from_IQF_label(K,label)
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
            I = ideal_from_IQF_label(K,label)
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

def parse_dims_line(L, K, d):
    """parse one line from a dimtabeis*newdims as output by
    make_dimtabnew().  Lines have 7 columns:

    d  field  (1, 2 or abs(disc))
    w  weight (2)
    level_label (HNF-style, e.g. 1.0.1)
    dimall      (dimension of full gl2-space, =dimcusp+dimeis)
    dimcusp     (dimension of cuspidal gl2-space)
    dimcuspnew  (dimension of new cuspidal gl2-space)
    dimeis      (dimension of eisenstein gl2-space)

    NB The LMFDB only stores cupidal dimensions, so the dimall and
    dimeis columns are ignored.

    Returns the label and a dict containing all the data.
    """
    dd, w, label, dimall, dimcusp, dimcuspnew, dimeis = L.split()
    assert int(dd) == d
    assert int(w) == 2
    N = ideal_from_IQF_label(K,label)
    level_label = ideal_label(N)
    D = K.discriminant().abs()
    field_label = "2.0.{}.1".format(D)
    label = "-".join([field_label, level_label])
    level_norm = N.norm()
    dimcuspnew = int(dimcuspnew)
    dimcusp = int(dimcusp)
    gl2_dims = {'2': {'new_dim': dimcuspnew, 'cuspidal_dim': dimcusp}}
    gl2_new_totaldim = dimcuspnew
    gl2_cusp_totaldim = dimcusp
    return label, {
        'field_label': field_label,
        'field_absdisc': D,
        'level_norm': level_norm,
        'level_label': level_label,
        'label': label,
        'gl2_dims': gl2_dims,
        'gl2_new_totaldim': gl2_new_totaldim,
        'gl2_cusp_totaldim': gl2_cusp_totaldim,
    }

def read_dimtabeis_new(d, fname):
    """Read a file dimtabeis*newdims as output by make_dimtabnew().
    Ignore lines starting with #.

    NB This script is intended for use with gl2, weight 2, data, so
    the records will not have the keys 'sl2_dims', 'sl2_new_totaldim'
    or 'sl2_cusp_totaldim'.  When uploading we'll need to make sure
    not to overwrite the sl2 data which exists in the database.

    Returns a dict containing all the data, keys all space labels with
    values the data for that space.

    """
    K, a = field(d)
    D = K.discriminant().abs()
    field_label = "2.0.{}.1".format(D)
    data = {}
    filename = os.path.join(DIM_DIR, fname)
    with open(filename) as infile:
        for L in infile:
            if L[0] != '#':
                label, record = parse_dims_line(L, K, d)
                data[label] = record
    return data

def numerify_iso_label(lab):
    from sage.databases.cremona import class_to_int
    return class_to_int(lab.lower())

def parse_newforms_line(line, K):
    r""" Parses one line from a newforms file.  Returns a complete entry
    for the forms collection. This is only for newforms of weight 2.

    Input line fields:

    field_label level_label level_suffix level_ideal wt bc cm sign Lratio AL_eigs pol hecke_eigs

    Sample input line:

    2.0.4.1 65.18.1 a (7+4i) 2 0 ? 1 1 [1,-1] x [0,0,-1,-2,1,-4,0,6,6,-6,2,2,6,0,-4,6,-6,-10,-10,2,2,6,6,-10,8]

    NB We expect 3-component HNF-style labels for ideals (N.c.d) but will convert to LMFDB labels N.i on input.
    """
    data = split(line)
    # base field
    field_label = data[0]
    field_disc = - int(field_label.split(".")[2])
    field_deg = 2
    field_bad_primes = ZZ(field_disc).support()
    # Hecke field degree (=dimension)
    hecke_poly = data[10]
    dimension = 1 if hecke_poly == 'x' else int(Qx(hecke_poly).degree())
    # level
    N = ideal_from_IQF_label(K, data[1])
    level_label = ideal_label(N)
    level_norm = int(level_label.split(".")[0])
    level_ideal = data[3]
    level_gen =  level_ideal[1:-1] # strip (,)
    level_bad_primes = ZZ(level_norm).support()
    label_suffix = data[2]
    label_nsuffix = numerify_iso_label(label_suffix)
    if dimension>1:
        label_suffix = label_suffix+str(dimension)
    short_label = '-'.join([level_label, label_suffix])
    label = '-'.join([field_label, short_label])

    weight = int(data[4])
    assert weight == 2
    bc = data[5]
    if bc!='?': bc=int(bc)
    cm = data[6]
    if cm!='?': cm=int(cm)
    sfe = data[7] # sign
    if sfe!='?': sfe = int(sfe) # sign
    Lratio = data[8]   # string representing rational number
    try:
        AL_eigs = [int(x) for x in data[9][1:-1].split(",")]
    except ValueError:
        AL_eigs = [x for x in data[9][1:-1].split(",")]
    try:
        hecke_eigs = [int(x) for x in data[11][1:-1].split(",")]
    except ValueError:
        hecke_eigs = [x for x in data[11][1:-1].split(",")]

    return label, {
        'label': label,
        'field_label': field_label,
        'field_disc': field_disc,
        'field_deg': field_deg,
        'level_label': level_label,
        'level_norm': level_norm,
        'level_ideal': level_ideal,
        'level_gen': level_gen,
        'label_suffix': label_suffix,
        'label_nsuffix': label_nsuffix,
        'short_label': short_label,
        'dimension': dimension,
        'hecke_poly': hecke_poly,
        'weight': weight,
        'sfe': sfe,
        'Lratio': Lratio,
        'bc': bc,
        'CM': cm,
        'AL_eigs': AL_eigs,
        'hecke_eigs': hecke_eigs,
        'field_bad_primes': field_bad_primes,
        'level_bad_primes': level_bad_primes,
    }

def read_newforms(d, fname):
    """Read a newforms file as output by nflist.cc.
    Ignore lines starting with # or containing "Primes" or "...".

    Returns a dict containing all the data, keys all newform labels with
    values the data for that newform.
    """
    K, a = field(d)
    D = K.discriminant().abs()
    field_label = "2.0.{}.1".format(D)
    data = {}
    filename = os.path.join(FORM_DIR, fname)
    with open(filename) as infile:
        for L in infile:
            if L[0] != '#' and "Primes" not in L and "..." not in L:
                label, record = parse_newforms_line(L, K)
                data[label] = record
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

bmf_dims_schema_no_sl2 = {
    'gl2_dims': 'jsonb',
    'field_absdisc': 'integer',
    'label': 'text',
    'field_label': 'text',
    'level_norm': 'bigint',
    'level_label': 'text',
    'gl2_new_totaldim': 'integer',
    'gl2_cusp_totaldim': 'integer',
}

bmf_forms_schema = {
    'field_disc': 'smallint',
    'sfe': 'smallint',
    'level_gen': 'text',
    'weight': 'smallint',
    'CM': 'smallint',
    'bc': 'smallint',
    'field_deg': 'smallint',
    'label': 'text',
    'label_suffix': 'text',
    'hecke_poly': 'text',
    'field_label': 'text',
    'AL_eigs': 'jsonb', # usually a list of ints, +1 and -1, could be "?"
    'Lratio': 'text',   # string representing a rational number
    'level_norm': 'bigint',
    'hecke_eigs': 'jsonb', # usually a list of ints but could be a list of strings for dim>1
    'short_label': 'text',
    'level_label': 'text',
    'label_nsuffix': 'smallint',
    'dimension': 'smallint',
    'level_ideal': 'text',
    'field_bad_primes': 'integer[]',
    'level_bad_primes': 'integer[]'
}


int_cols = ['field_absdisc', 'level_norm', 'gl2_cusp_totaldim', 'gl2_new_totaldim', 'sl2_cusp_totaldim', 'sl2_new_totaldim',
            'field_disc', 'sfe', 'weight', 'CM', 'bc', 'field_deg', 'label_nsuffix', 'dimension']
str_cols = ['label', 'field_label', 'level_label',
            'level_gen', 'label_suffix', 'hecke_poly', 'Lratio', 'short_label', 'level_ideal']
dict_cols = ['gl2_dims', 'sl2_dims',
             'AL_eigs', 'hecke_eigs']
int_list_cols = ['field_bad_primes', 'level_bad_primes']

def encode_col(colname, col=None):
    if 'sl2' in colname:
        return r"\N"
    if colname in int_cols:
        return str(col)
    if colname in dict_cols:
        return str(col).replace("'", '"').replace(" ", "")
    if colname in int_list_cols:
        return str(col).replace("[", '{').replace("]", '}').replace(" ", "")
    return col

def one_bmf_line(record, table, sl2):
    schema = (bmf_dims_schema_no_sl2 if sl2 else bmf_dims_schema) if table == 'dims' else bmf_forms_schema
    cols = list(schema.keys())
    cols.remove('label')
    cols = ['label'] + cols
    return "|".join([encode_col(col, record.get(col, None)) for col in cols])

def write_bmf_upload_file(data, fname, table, sl2):
    """data is a dict as returned by read_dimtabeis_new() or
    read_newforms(), with keys level labels, each value a dict with
    data for the table.

    fname is the output file name, created in UPLOAD_DIR

    table is 'dims' or 'forms'

    NB For the 'dims' table, assume that sl2_levels holds a (possibly
    empty) list of levels for which the LMFDB table already has sl2
    dimension data which we do not want to over-write.  In this case
    we output two files, one for the sl2_levels and one for the rest.

    e.g. if the file sl2_levels_43 contains a list of the level labels
    for which we have sl2 dimension data, first do

    sage: sl2_levels = read_data("sl2_levels_43", str)

    """
    assert table in ['dims', 'forms']
    schema = (bmf_dims_schema_no_sl2 if sl2 else bmf_dims_schema) if table == 'dims' else bmf_forms_schema
    cols = list(schema.keys())
    cols.remove('label')
    cols = ['label'] + cols
    vals = [schema[k] for k in cols]
    print("cols: {}".format(cols))

    filename = os.path.join(UPLOAD_DIR, fname)
    with open(filename, 'w') as outfile:

        # Write header lines

        outfile.write("|".join(cols))
        outfile.write("\n")
        outfile.write("|".join(vals))
        outfile.write("\n\n")

        # Write data lines

        nlines = 0
        for label, record in data.items():
            if ((table=='forms')
                or (sl2 and record['level_label'] in sl2_levels)
                or (not sl2  and record['level_label'] not in sl2_levels)):
                outfile.write(one_bmf_line(record, table, sl2) + "\n")
                nlines +=1
                if nlines%1000 == 0:
                    print("{} lines written so far to {}".format(nlines, filename))

    print("{} lines written to {}".format(nlines, filename))

# e.g. (assumes directory ~/bmf-upload exists)
# sage: %runfile bianchi.py
# sage: d=43
# sage: N1=1
# sage: N2=1000
# sage: dimdat = read_dimtabeis_new(d, "dimtabeis.{}.all.newdims".format(d))
# sage: sl2_levels = []
# sage: if d in [1,2,3,7,11,19,43,67,163,20]:
# sage:    sl2_levels = read_data("sl2_levels_{}".format(d), str)
# sage:    write_bmf_upload_file(dimdat, "bmf_dims.{}.{}-{}.sl2".format(d,N1,N1), 'dims', sl2=True)
# sage: write_bmf_upload_file(dimdat, "bmf_dims.{}.{}-{}.no_sl2".format(d,N1,N2), 'dims', sl2=False)
# sage: formdat = read_newforms(d, "newforms.{}.{}-{}".format(d,N1,N2))
# sage: write_bmf_upload_file(formdat, "bmf_forms.{}.{}-{}".format(d,N1,N2), 'forms', True)

# Copy the three files to legendre in bmf-upload/.
# Do the upload as follows:
# sage: from lmfdb import db
# sage: d=43
# sage: N1=1
# sage: N2=1000
# sage: db.bmf_forms.copy_from("/scratch/home/jcremona/bmf-upload/bmf_forms.{}.{}-{}".format(d,N1,N2))
# sage: db.bmf_dims.update_from_file("/scratch/home/jcremona/bmf-upload/bmf_dims.{}.{}-{}.sl2".format(d,N1,N2))
# sage: db.bmf_dims.copy_from("/scratch/home/jcremona/bmf-upload/bmf_dims.{}.{}-{}.no_sl2".format(d,N1,N2))
