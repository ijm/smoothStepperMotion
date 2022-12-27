import numpy as np

from numpy.polynomial.polynomial import *
from fractions import Fraction

def toLatex(x) :
  return x.toLatex() if hasattr(x,'toLatex') else (
    str(int(x)) if type(x) == float and int(x)==x else str(x) )

Fraction.toLatex = lambda self : (
   toLatex(self.numerator) if self.denominator == 1 else
     ("-" if self.numerator < 0 else "") +
     "\\frac{" + str(abs(self.numerator)) + "}{" + str(self.denominator) +"}" )



def polyToLatex(self,n) :
  def signedToLatex(x) :
    return ("" if x<0 else "+") + toLatex(x)

  def coef(x) :
    s = signedToLatex(x)
    return s[0:1] if abs(x) == 1.0 else s

  if n == None :
    n = self.varLabel 

  s = " ".join( reversed( 
    [ signedToLatex(x)                 if i == 0 else (
      coef(x) + n if i == 1 else (
      coef(x) + n + "^{"+toLatex(i)+"}"))
        for i,x in enumerate(self.coef) if x != 0 ]
   ))
  if s == "" :
    s = "0"
  return s[1:] if s[0] == '+' else s

Polynomial.varLabel = "x"
Polynomial.toLatex  = polyToLatex


