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
import matplotlib.ticker as ticker

from toLatex import *
from helpers import *

colorList = [ 'red', 'green', 'blue', 'magenta', 'yellow', 'cyan' ]

def doArgs() :
    parser = ArgumentParser(description= "Plot Motor profile for SL Door and Strike Motor")
    parser.add_argument('-i', '--input', dest='iname', help='profile name (.in)')
    parser.add_argument('-p', '--polys', dest='polys', help='polynomials file')
    parser.add_argument('-t', '--template', dest='tmpl', help='profile template basename (for .c++tmpl and .textmpl)')
    parser.add_argument('-m', '--motordata', dest='motors',default=None,
            help='motor datasheet torque velocity table file name')
    return parser.parse_args()

def runProfile(params, rows):
    d = dict()
    for k,v in params.items() :
        s = k.split(".")
        if len(s) == 1 :
            d[k] = v
        elif len(s) == 2 :
            e,n = s
            if e not in d :
                d[e] = dict()
            d[e][n] = v
        else :
            raise Exception(f"odd type of environment name ({k})")
    
    nrows = [ x[0:5] + (" ".join(x[5:]),) for x in rows]
    return d, nrows


def checkParam(d, n, k):
    if k not in d :
        raise Exception (f"Parameter {n}{k} is missing and required.")
   
def makeProfile( x,v,t, Dx, Dv, Dt, p_f, p_f_star, env, gl):
    p_s = Polynomial(  [ -2.0*t/Dt-1, 2.0/Dt ] )
    p_ds = p_s.deriv()(0)

    pend   = p_f_star(1)
    alpha = 1 if pend == 0 else abs(Dx / pend) * p_ds
  
    return (x,x+Dx, v,v+Dv, t,t+Dt, p_s, alpha * np.sign(Dx) )

def buildProfiles(params,rows, polys, integs):
    if "global" not in params :
        params["global"] = dict()

    g = params["global"]

    checkParam(g, 'global.', 'pullyDia')

    p = []
    x1 = 0.0
    v1 = 0.0
    t1 = 0.0

    x_min = 0.0
    x_max = 0.0
    v_min = 0.0
    v_max = 0.0

    envused = dict()
    polyused = dict()

    for (Dx, Dv, Dt, prof, env, desc) in rows :
      checkParam(params, '', env)
      checkParam(polys, '', prof)
      envused[env] = 1
      polyused[prof] = 1
      
      x0,x1,v0,v1,t0,t1,p_s,alpha = makeProfile(x1,v1,t1, Dx, Dv, Dt,
                                      polys[prof], integs[prof], params[env], g)
      p.append( (x0,x1, v0,v1, t0,t1,  p_s, alpha, prof, env, desc) )
      if x1 > x_max :
          x_max = x1
      if x1 < x_min :
          x_min = x1
  
    for i,k in enumerate(envused.keys()) :
        envused[k] = i

    for i,k in enumerate(polyused.keys()) :
        polyused[k] = i

    return (p, envused, polyused, x_min, x_max )

def calcLocalPolys(p, polys, integs, params ):
    m = params["global"]["baseMass"]
    dia = params["global"]["pullyDia"]
    q = []
    roots = []
    p_tau0 = 0.
    p_tau1 = 0.
    p_rps0 = 0.
    p_rps1 = 0.
  
    for (x0,x1, v0,v1, t0,t1, p_s, alpha, prof, env, desc) in p :
        e     = params[env]
        p_dx  = alpha * polys [prof](p_s)
        p_x   = alpha * integs[prof](p_s) / (p_s.deriv()(0)) + x0
        p_ddx = p_dx.deriv()
    
        F_spring = e["springK"] * ( p_x - e["springE0"] )
        F_mass   = (m + e["extraMass"]) * p_ddx
        f_F      = F_spring + F_mass
    
        circ = np.pi * dia
        f_tau = f_F * 0.5 * dia
        f_rps = p_dx / circ
    
        roots = p_ddx.deriv().roots()
        p_tau0 = np.min( list(f_tau(roots))+[p_tau0,] )
        p_tau1 = np.max( list(f_tau(roots))+[p_tau1,] )
    
        p_rps0 = np.min( list(f_rps(roots))+[p_rps0,] )
        p_rps1 = np.max( list(f_rps(roots))+[p_rps1,] )
    
        q.append( (p_x, p_dx, p_ddx, F_spring, F_mass, f_F, f_tau, f_rps,
               x0,x1, v0,v1, t0,t1, p_s, alpha, prof, env, desc) )
  
    return q, (p_tau0, p_tau1, p_rps0, p_rps1)

def plotReqProfile(p, eu, x_min, x_max, cl):
    t_min = 0.0
    t_max = p[-1][13]

    f,axs = plt.subplots(3,1)

    for ax in axs :
        ax.set_xlim( t_min, t_max )
        ax.grid(True, linestyle=':')

    ax1b = axs[1].twinx()
    ax1b.set_xlim( t_min, t_max )

    for (p_x, p_dx, p_ddx, F_spring, F_mass, f_F, f_tau, f_rps,
             x0,x1, v0,v1, t0,t1, p_s, alpha, prof, env, desc) in p:

        ts = np.linspace(t0, t1, 19)

        axs[0].plot( t1, x1, 'o', c='black' )
        axs[0].plot( ts, p_x(ts) , c='black')

        axs[1].plot( ts, p_dx(ts) , c='green')
        ax1b.plot( ts, p_ddx(ts) , c='blue')

        axs[2].plot( ts, f_F(ts) , c='red' )
        axs[2].plot( ts, F_spring(ts) , '--', c='m', lw = 1, alpha = 0.66 )
        axs[2].plot( ts, F_mass(ts) , '--', c='brown', lw = 1, alpha = 0.66 )
        for ax in axs :
          ax.axvspan(t0,t1, facecolor = colorList[eu[env]], ec='none', alpha = 0.25 )

    for (n,v,t) in cl:
        axs[0].plot( t, v, '+', c='white', mew=1.5, ms=5 )

    #axs[0].set_xlabel("Time (s)")
    axs[0].set_ylim( x_min*1.1, x_max*1.1 )
    axs[0].yaxis.set_major_locator(ticker.MultipleLocator(0.05))

    axs[0].xaxis.tick_top()
    axs[0].xaxis.set_label_position("top")
    axs[0].xaxis.set_tick_params(top=True, bottom=False)
    plt.setp(axs[0].get_xticklabels(), visible=True)

    axs[1].set_ylim( -0.5, 0.5 ) 
    axs[1].yaxis.set_major_locator(ticker.MultipleLocator(0.2))
    ax1b.set_ylim( -8,8 )
    ax1b.yaxis.set_major_locator(ticker.MultipleLocator(2))
    plt.setp(axs[1].get_xticklabels(), visible=False)

    axs[2].yaxis.set_major_locator(ticker.MultipleLocator(5))
    axs[2].set_xlabel("Time ($\mathrm{s})")
    axs[2].xaxis.tick_bottom()
    axs[2].xaxis.set_tick_params(top=False, bottom=True)
    plt.setp(axs[2].get_xticklabels(), visible=True)


    axs[0].set_ylabel("Displacement ($\mathrm{m}$)")
    axs[1].set_ylabel("Velocity ($\mathrm{m}\cdot\mathrm{s}^{-1}$)")
    ax1b.set_ylabel("Acceleration ($\mathrm{m}\cdot\mathrm{s}^{-2}$)")
    axs[2].set_ylabel("Force ($\mathrm{N}$)")
    #axs[0].set_title("Requested Displacement/Time Profile")

    f.subplots_adjust(hspace=0,wspace=0)
    f.set_size_inches(7.5, 4.5)

    return f

def plotTorqueVelocity(q, m):
    f,ax = plt.subplots(1,1)

    ax.set_ylim( 0.0, 0.8 )
    ax.set_xlim( 0.0, 40.0 )
    ax.xaxis.set_major_locator(ticker.MultipleLocator(5))

    for (p_x, p_dx, p_ddx, F_spring, F_mass, f_F, f_tau, f_rps,
             x0,x1, v0,v1, t0,t1, p_s, alpha, prof, env, desc) in q:
        ts = np.linspace(t0, t1, 29)

        ax.plot( abs(f_rps(ts)) , abs(f_tau( ts )), c='black' )

    i = 0  
    mnames = dict()

    for d,dr in m :
        c = colorList[i]
        mnames[ d["name"] ] = c
        dxs = [ line[int(d["rps"])] for line in dr ]
        dys = [ line[int(d["torque"])] * d["conv"] for line in dr ]

        ax.plot( dxs, dys , "--", c=c )

        i=i+1

    ax.set_ylabel("Torque ($\mathrm{N}\cdot\mathrm{m}$)")
    ax.set_xlabel("Velocity ($\mathrm{s}^{-1}$)")
    ax.grid(True, linestyle=':')

    f.set_size_inches(3.5, 3.5)
  
    return f, mnames

def findZones(p, params, pused, polys, integs):
    Dt = params["smallestDt"]

    v_max = 0
    l_max = 0
    d = dict()
    for k,v in pused.items():
        d[v] = (list( reFracPoly(integs[k]).coef ),
                list( reFracPoly(polys[k]).coef ),
                k)
        if v > v_max :
            v_max = v
        l = len( integs[k] )
        if l > l_max :
            l_max = l
  
    parrs = []
  
    l00 = [Fraction(0)] * l_max
  
    for k in range(v_max+1):
        parrs.append( ( (d[k][0] + l00)[0:l],
                        (d[k][1] + l00)[0:l],
                         d[k][2] ))
  
    inx = dict()
    zones = []
  
    def intChk(x,Dt) :
        xp = x / Dt
        xpi = int(xp)
        if xp != xpi :
            raise Exception( f"Time code {x} isn't on a {Dt} boundry." ) 
        return xpi
  
    for i, (x0,x1, v0,v1, t0,t1, p_s, alpha, prof, env, desc) in enumerate(p) :
        for j in range( intChk(t0,Dt), intChk(t1,Dt)) :
            inx[j] = i
        psf = reFracPoly(p_s).coef
        zones.append(( psf[0].numerator, psf[0].denominator, 
                       psf[1].numerator, psf[1].denominator, 
                       alpha,
                       pused[prof], prof, t0, t1 ) )
  
    inxarray = [ inx[j] for j in range(int(p[-1][5]/Dt)) ] 
  
    return (parrs, zones, inxarray)

def genCCode(name, parrs, zones, inxa, g, tmpl) :
    M = len(parrs[0][0])
    zonefmt = "  {{ {0:.1f}/{1:.1f}, {2:.1f}/{3:.1f}, {4}, {5} }} /* {6}: ({7:.2f} < t < {8:.2f}) */"
  
    return tmpl.format(
        name       = name,
        Dt         = g["smallestDt"],

        N          = len(parrs),
        M          = M,
        L          = len(zones),
        K          = len(inxa),

        kMotor     = 2.0 * g["Nsteps"] / (np.pi * g["pullyDia"]),

        Xcoefs     = ",\n".join( [ "  {"+ ", ".join([ 
                     "{}/{}".format(float(x.numerator), float(x.denominator))
                         for x in nX ]) + "} /* " + cmt + " */" 
                             for (nX,nDx,cmt) in parrs ] ),

        DXcoefs    = ",\n".join( [ "  {"+ ", ".join( [ 
                     "{}/{}".format(float(x.numerator), float(x.denominator))
                         for x in nDx ]) + "} /* " + cmt + " */" 
                             for (nX,nDx,cmt) in parrs ] ),

        zones      = ",\n".join( [ zonefmt.format(*x) for x in zones]),

        inxa       = ", ".join([ str(x) for x in inxa ]),

        spoly      = "("*M + "*s + ".join( [ "c[{0}])".format(i) for i in reversed(range(M)) ] )
        )

def commandAnd(l):
    if len(l) == 0 :
        return ""
    if len(l) == 1 :
        return l[0]
    return ", ".join(l[:-1]) + " and " + l[-1]

def buildDoc(name, mnames, eUsed, mmTauRps, cpr, params, tmpl):
    inx = dict()

    for k,v in eUsed.items() :
        params[k]["zoneName"] = k
        params[k]["zoneN"] = v
        params[k]["zoneColor"] = colorList[v]
        inx[v] = k
  
    maxTau = np.max( (-mmTauRps[0], mmTauRps[1]) ),
    maxRps = np.max( (-mmTauRps[2], mmTauRps[3]) ),
  
    zoneRows = "\\\\\n".join ([ 
      """{zoneN}&{zoneName}&{zoneColor}&${springK}$&
      ${springE0}$&${friction}$&${extraMass}$""".format(**(params[inx[k]]))
      for k in sorted(inx.keys()) ])
  
    def sofl(l) :
        return ", ".join( [ f"{x:.4}" for x in l ]  )
  
    calRows = "\\\\\n".join ([ f"{n} & ${v:.4}$ & {sofl(t)}" for (n, (v,t)) in cpr.items() ] )
  
    return tmpl.format(
        name       = name,
        motorNames = commandAnd([ f"{k} ({v})" for (k,v) in mnames.items()]),
        zoneN      = sorted(inx.keys())[-1]+1,
        zoneRows   = zoneRows,
        calRows    = calRows,
        mass       = params["global"]["baseMass"],
        dia        = params["global"]["pullyDia"],
        maxTau     = np.max( (-mmTauRps[0], mmTauRps[1]) ),
        maxRps     = np.max( (-mmTauRps[2], mmTauRps[3]) ),
        )
  
def calcCalPoints(p, params):
    dt = params["global"]["smallestDt"] * 0.1
    calPoints = params["calpoint"] if "calpoint" in params else dict()

    t_min = 0.0
    t_max = p[-1][13]


    def test(s,e) :
        return [ (n,v) for (n,v) in calPoints.items() if s <= v and v <= e ]

    c = []
    for (p_x, p_dx, p_ddx, F_spring, F_mass, f_F, f_tau, f_rps,
               x0,x1, v0,v1, t0,t1, p_s, alpha, prof, env, desc) in p:
        ts = np.arange(t0,t1, dt)
        xs = p_x(ts)
        for i in range(len(xs)-1):
            l = test( xs[i], xs[i+1] )
            for (n,v) in l :
                c.append( (n, v, ts[i], ts[i+1], p_x, p_dx, p_ddx) )
  
    def solve(f_x, f_dx, f_ddx, x, y0) :
        f_x = f_x - y0
  
        yp = 1.0
  
        while abs(yp) > 1.0e-10 :
            yp = f_x(x)
            dp = f_dx(x)
            x = x - ( yp / ( dp - ( 0.5 * yp * f_ddx(x)/ dp)))
  
        # Newtowns : x1 = x0 - (f_x(x0)-y0)/f_dx(x0)
        return x
  
    c = [ (n,v, solve( p_x, p_dx, p_ddx, (t0+t1)*0.5, v)) for (n,v,t0,t1,p_x,p_dx,p_ddx) in c ]
  
    cpr = dict()
    for n,v in calPoints.items() :
        cpr[n] = (v,[])
  
    for (n,v,t) in c :
        cpr[n][1].append(t)
  
    return cpr, c
 
def main():
    args = doArgs()

    bits = args.iname.split(".")
    iname = ".".join( bits[:-1] if len(bits)>1 else bits ) 

    #*\label{code:\lstname-read-profiles}*
    motors = [ readProfile(args.motors) ] if args.motors else []
    rootsS, polys, integs, pInx = initPolys(*readProfile(args.polys))
    params, row   = runProfile(*readProfile(args.iname))
 
    (p, eUsed, pUsed, x_min, x_max) = buildProfiles(params, row, polys, integs)

    #*\label{code:\lstname-calculate}*
    q, mmTauRps = calcLocalPolys(p, polys, integs, params )
    cpr, cl = calcCalPoints(q, params)

    #*\label{code:\lstname-plot}*
    fig = plotReqProfile(q, eUsed, x_min, x_max, cl)
    fig.savefig( f"{iname}.disp.pgf", bbox_inches = 'tight')
    fig.savefig( f"{iname}.disp.pdf", bbox_inches = 'tight')

    fig, mnames = plotTorqueVelocity(q, motors)
    fig.savefig( f"{iname}.torque.pgf", bbox_inches = 'tight')
    fig.savefig( f"{iname}.torque.pdf", bbox_inches = 'tight')

    (parrs, zones, inxa) = findZones(p, params["global"], pInx, polys, integs)

    #*\label{code:\lstname-gencode}*
    with open( f"{args.tmpl}.c++tmpl", "r" ) as f, open(iname+".c", "w") as g :
        g.write( genCCode( iname, parrs, zones, inxa, params["global"], f.read() ))

    #*\label{code:\lstname-gendoc}*
    with open( f"{args.tmpl}.textmpl", "r" ) as f, open( f"{iname}.section.tex", "w") as g :
        g.write( buildDoc(iname, mnames, eUsed, mmTauRps, cpr, params, f.read()) )

main() 
