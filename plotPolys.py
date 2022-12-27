# Plot Motor profile for SL Door and Strike motor

import numpy as np
from argparse import ArgumentParser

import matplotlib as mpl
mpl.use('pgf')
pgf_with_pdflatex = {
    "pgf.texsystem": "pdflatex",
    "pgf.preamble": "\n".join( [
         r"\usepackage[utf8x]{inputenc}",
         r"\usepackage[T1]{fontenc}",
         ] )
}
         #r"\usepackage{cmbright}",
mpl.rcParams.update({'font.size': 9, "font.family": "serif",})
mpl.rcParams.update(pgf_with_pdflatex)
import matplotlib.pyplot as plt

from numpy.polynomial.polynomial import *
from fractions import Fraction

from toLatex import *
from helpers import *

def doArgs() :
    parser = ArgumentParser(description= "Plot Motor profiles")
    parser.add_argument('-p', '--polys', dest='pname',required=True, help='polynomials file')
    return parser.parse_args()

def plotPolys(name, polys, integs, pUsed):
    s=""
    for k in pUsed.keys() :
        n = f"{name}.{k}"

        p_f = polys[k]
        p_f_star = integs[k]
        p_df = p_f.deriv()

        ps = (k,n) + tuple( [ 
                 reFracPoly(p).toLatex("t") for p in [p_f_star, p_f, p_df ]] )
    
        s += """%
\\parbox{{18mm}}{{Polynomial\\newline "{0}"}} &
\\begin{{minipage}}[c]{{0.4\\textwidth}}
\\input{{{1}.pgf}}
\\end{{minipage}} &\\vspace{{10mm}}
\color{{darkgreen}}$x(t)={2}$\\newline\\newline
\color{{blue}}$\dot{{x}}(t)={3}$\\newline\\newline
\color{{red}}$\ddot{{x}}(t)={4}$ \\\\ 
""".format( *ps )

        f, ax = plt.subplots(1,1)
        xs = np.linspace(-1,1,49)
    
        ax.set_ylim( -2.6, 2.6 )
        ax.set_xlim( -1.0, 1.0 )
    
        ax.plot( xs, p_f(xs) )
        ax.plot( xs, p_f_star(xs) )
        ax.plot( xs, p_df(xs) )
    
        #ax.set_xlabel("Normalized Time")
        #ax.set_ylabel("Relative Disp.,\nVel. and Accel.")
        #ax.set_title( "$\\dot{x}\\left(t\\right)="+reFracPoly(p_f).toLatex()+"$",y=1.08 )
        #ax.set_title( k )
        ax.grid(True) 
    
        f.set_size_inches(3, 1.5)

        f.savefig( f"{name}.{k}.pgf", bbox_inches='tight')
        f.savefig( f"{name}.{k}.pdf", bbox_inches='tight')
    

    return f, s

def main():
    args = doArgs()

    bits = args.pname.split(".")
    name = ".".join( bits[:-1] if len(bits)>1 else bits ) 

    aroots, polys, integs, pInx = initPolys( *readProfile(args.pname) )

    fig, tex = plotPolys(name, polys, integs, pInx)

    with open( f"{name}.tex","w") as g:
        g.write(tex)


main() 


