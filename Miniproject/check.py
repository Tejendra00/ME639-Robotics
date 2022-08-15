import sympy as sp
from sympy import *

t=sp.symbols("t")
q1=sp.symbols('q1')
q2=sp.symbols('q2')
q1=sp.Function('q1')
q2=sp.Function('q2')
l1=sp.symbols('l1')
l2=sp.symbols('l2')

x2=l1*cos(q1(t))+l2*cos(q2(t))
dx2=x2.diff(t)
ddx2=dx2.diff(t)
print(ddx2)

y2=l1*sin(q1(t))+l2*sin(q2(t))
dy2=y2.diff(t)
ddy2=dy2.diff(t)
print(ddy2)