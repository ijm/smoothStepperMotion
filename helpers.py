import numpy as np
from numpy.polynomial.polynomial import *
from fractions import Fraction

def floatOrString(s) :
    try:
        s = float(s)
    except:
        pass
    return s
  
def readProfile(fname):
    with open(fname, 'r') as f :
        params = dict()
        rows = []
        for l in f :
            l = l.split('#')[0]
            if l.strip() == "" :
                continue
            bits = l.strip().split()
            if type(floatOrString(bits[0])) == float :
                rows.append( tuple( [ floatOrString(x) for x in bits] ))
            else:
                params[ bits[0] ] = None if len(bits)==1 else (
                    floatOrString(bits[1]) if len(bits)==2 else (
                        tuple( [ floatOrString(x) for x in bits[1:] ] )
                ))
    return params, rows

def reFracPoly(p) :
    return Polynomial ( [ Fraction(x).limit_denominator() for x in p.coef] )


def initPolys(params, rows):
    def normalize(p, x0, x1, r):
        p_min = p(x0) if type(x0)==float else np.min( p(r) )
        p_max = p(x1) if type(x1)==float else np.max( p(r) )
    
        rng = 1.0/(p_max - p_min) if p_max != p_min else 1.0
  
        return rng * (p - p_min)

    polys   = dict()
    integs  = dict()
    polyinx = dict()
    aroots  = dict()
    for (i, (n,r)) in enumerate(params.items()):
        roots = r[2:]
        if len(r) > 2 :
            p_df     = np.sign(r[2]) * Polynomial.fromroots( roots )
            p_f      = normalize( p_df.integ(), r[0], r[1], np.array(roots) )
            p_f_star = p_f.integ()
            p_f_star = p_f_star - p_f_star(-1)
        else:
            p_f = Polynomial(0)
            p_f_star = Polynomial(0)
        aroots[n]   = roots
        polys[n]   = p_f
        integs[n]  = p_f_star
        polyinx[n] = i
  
    return aroots, polys, integs, polyinx

