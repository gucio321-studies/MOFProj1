#!/usr/bin/env python3
import numpy as np
import pint
import numba as nb
from matplotlib import pyplot as plt
from numba.np.arrayobj import default_lt

#si = pint.UnitRegistry()
m = 1 #* si.kg

@nb.njit
def V(x):
    return -np.exp(-((x+2)/3)**2) - 1.2 * np.exp(-((x-2)/3)**2)# * si.J

def devV(xn, delta_x):
    return (V(xn + delta_x) - V(xn))/(2*delta_x)

vs = [0]
xs = [0]
delX = 0.001# * si.m
delT = 0.1# * si.s

def v(n):
   if len(vs) >= n+1:
       return vs[n]
   while len(vs) < n+1:
       last = len(vs) - 1
       vs.append(vs[last] - 1/m * devV(x(last), delX))
   return v(n)

def x(n):
    if len(xs) >= n+1:
        return xs[n]
    while len(xs) < n+1:
        last = len(xs) - 1
        xs.append(xs[last] - 1/m * devV(x(last), delX))
    return x(n)



# ex 1
times = np.linspace(0, 50, int(50/0.1))
# calculate appropiate number of x's and v's
x(len(times)-1)
v(len(times)-1)
# integrate
x_t = [0]
for x in xs:
    x_t.append(x_t[-1] + x)
# delete the first 0
x_t.pop(0)

# do the same for vs
v_t = [0]
for x in xs:
    v_t.append(v_t[-1] + x)
v_t.pop(0)

plt.subplot(1, 2, 1)
plt.plot(times, x_t)
plt.xlabel("Time [s]")
plt.ylabel("x [m]")

plt.subplot(1, 2, 2)
plt.plot(times, v_t)
plt.xlabel("Time [s]")
plt.ylabel("v [m/s]")

plt.show()