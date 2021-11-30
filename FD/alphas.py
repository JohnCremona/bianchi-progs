# Precomputed alphas for some fields

def precomputed_alphas(k_or_d):

    from sage.all import ZZ, QQ, polygen, NumberField, NFCusp
    from utils import smallest_ideal_class_representatives

    if k_or_d in ZZ:
        # then d=k>0 should be square-free
        x = polygen(QQ)
        d = k_or_d
        if d%4==3:
            k = NumberField(x**2-x+(1+d)//4, 'w')
        else:
            k = NumberField(x**2+d, 'w')
    else:
        k = k_or_d
        d = -k.discriminant().squarefree_part()

    w = k.gen()
    IC = smallest_ideal_class_representatives(k)
    cusp = lambda a: NFCusp(k, k(a), lreps=IC)

    alphas_dict = {}
    alphas_dict[1] = [cusp(0)]
    alphas_dict[2] = [cusp(0)]
    alphas_dict[3] = [cusp(0)]
    alphas_dict[7] = [cusp(0)]
    alphas_dict[11] = [cusp(0)]
    alphas_dict[19] = [cusp(a) for a in [0, w/2,(w-1)/2]]

    alphas_dict[43] = [cusp(a) for a in [0, w/2, (-1+w)/2, w/3, -w/3,
                                         (1-w)/3, (-1+w)/3, (1+w)/3, (-1-w)/3]]

    alphas_dict[5] = [cusp(a) for a in [0, w/2, (w-4)/(2*w), (4-w)/(2*w)]]

    alphas_dict[6] = [cusp(a) for a in [0, (w+1)/2, 5/(2*w), -5/(2*w)]]

    alphas_dict[10] = [cusp(a) for a in [0, (w+1)/2, w/3, -w/3, (1+w)/3, (-1-w)/3, (1-w)/3,
                                         (-1+w)/3, 9/(2*w), -9/(2*w)]]

    alphas_dict[23] = [cusp(a) for a in [0, (1+2*w)/4, (-1-2*w)/4,
                                         (1+w)/3, (-1-w)/3, (2-w)/(1+w),
                                         (-2+w)/(1+w), (1+w)/(2-w),
                                         (-1-w)/(2-w), (-3+w)/(2+w),
                                         (3-w)/(2+w), (-2-w)/(3-w),
                                         (2+w)/(3-w)]]

    alphas_dict[31] = [cusp(a) for a in [0, (1+2*w)/4, (-1-2*w)/4, w/3,
                                         (-w)/3, (1-w)/3, (-1+w)/3,
                                         (1+w)/3, (-1-w)/3, 3/w, (-3)/w,
                                         3/(1-w), (-3)/(1-w), 3/(1+w),
                                         (-3)/(1+w), 3/(2-w), -3/(2-w),
                                         (-6+w)/(3+w), (6-w)/(3+w),
                                         (5+w)/(4-w), (-5-w)/(4-w)]]

    alphas_dict[67] = [cusp(a) for a in [0, w/2, (-1+w)/2, w/3, (-w)/3,
                                         (1-w)/3, (-1+w)/3, (1+w)/3, (-1-w)/3, w/4, (-w)/4, (-1+w)/4, (1-w)/4,
                                         (1+w)/4, (-1-w)/4, (2-w)/4, (-2+w)/4, (6+w)/(3-w), (-6-w)/(3-w),
                                         (2+w)/(3-w), (-2-w)/(3-w), (7-w)/(2+w), (-7+w)/(2+w), (3-w)/(2+w),
                                         (-3+w)/(2+w)]]

    alphas_dict[163] = [cusp(a) for a in [0, (w)/2, (-1+w)/2, (w)/3,
                                          (-w)/3, (1-w)/3, (-1+w)/3, (1+w)/3, (-1-w)/3, (w)/4, (-w)/4, (-1+w)/4,
                                          (1-w)/4, (1+w)/4, (-1-w)/4, (2-w)/4, (-2+w)/4, (w)/5, (-w)/5,
                                          (-1+w)/5, (1-w)/5, (2*w)/5, (-2*w)/5, (2-2*w)/5, (-2+2*w)/5, (2+w)/5,
                                          (-2-w)/5, (1-2*w)/5, (-1+2*w)/5, (-2+w)/5, (2-w)/5, (2+2*w)/5, (-2-2*w)/5,
                                          (1+w)/5, (-1-w)/5, (1+2*w)/5, (-1-2*w)/5, (w)/6, (-w)/6, (1-w)/6,
                                          (-1+w)/6, (1+w)/6, (-1-w)/6, (-2+w)/6, (2-w)/6, (2+w)/6, (-2-w)/6,
                                          (3-w)/6, (-3+w)/6, (12)/(w), (-12)/(w), (17)/(w), (-17)/(w),
                                          (12)/(1-w), (-12)/(1-w), (17)/(1-w), (-17)/(1-w), (12)/(1+w),
                                          (-12)/(1+w), (-17+w)/(1+w), (17-w)/(1+w), (12)/(2-w), (-12)/(2-w),
                                          (-16-w)/(2-w), (16+w)/(2-w), 7/(2+w), (-7)/(2+w), (18-w)/(2+w),
                                          (-18+w)/(2+w), (-11+w)/(2+w), (11-w)/(2+w), (-16+w)/(2+w),
                                          (16-w)/(2+w), 7/(3-w), (-7)/(3-w), (17+w)/(3-w), (-17-w)/(3-w),
                                          (-10-w)/(3-w), (10+w)/(3-w), (-15-w)/(3-w), (15+w)/(3-w), (3+w)/7,
                                          (-3-w)/7, (-1+2*w)/7, (1-2*w)/7, (1+2*w)/7, (-1-2*w)/7, (3-2*w)/7,
                                          (-3+2*w)/7, (2+3*w)/7, (-2-3*w)/7, (5-w)/(3+w), (-5+w)/(3+w),
                                          (-17+w)/(3+w), (17-w)/(3+w), (4+w)/(4-w), (-4-w)/(4-w), (-16-w)/(4-w),
                                          (16+w)/(4-w)]]

    if d in alphas_dict:
        return alphas_dict[d]

    print("No precomputed alphas for field {}".format(k))
    return None

