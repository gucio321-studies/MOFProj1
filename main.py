#!/usr/bin/env python3
import numpy as np
import pint
import numba as nb
from matplotlib import pyplot as plt
from numba.np.arrayobj import default_lt

#si = pint.UnitRegistry()
m = 1 #* si.kg

class experiment:
    def __init__(self, m, delX, delT):
        self.vs = [0]
        self.xs = [0]
        self.m = m
        self.delX = delX
        self.delT = delT

    def V(self, x):
        return -np.exp(-((x+2)/3)**2) - 1.2 * np.exp(-((x-2)/3)**2)# * si.J

    def devV(self, xn, delta_x):
        return (self.V(xn + delta_x) - self.V(xn))/(2*delta_x)
    def v(self, n):
       if len(self.vs) >= n+1:
           return self.vs[n]
       while len(self.vs) < n+1:
           last = len(self.vs) - 1
           self.vs.append(self.vs[last] - 1/m * self.devV(self.x(last), delX)*self.delT)
       return self.v(n)

    def x(self,n):
        if len(self.xs) >= n+1:
            return self.xs[n]
        while len(self.xs) < n+1:
            last = len(self.xs) - 1
            self.xs.append(self.xs[last] + self.v(last)*self.delT)
        return self.x(n)

    def precalc(self, n):
        """
        makes it so that x(n) is available instantly
        :param n:
        :return:
        """
        self.x(n)
        self.v(n)

    def get_xs(self):
        return self.xs
    def get_vs(self):
        return self.vs
    def plot(self, stop_time):
        """
        makes plots from ex 1 assuming start time = 0 and step delT
        :param stop_time:
        :return:
        """
        times =  np.linspace(0, stop_time, int(stop_time / self.delT))

        self.precalc(len(times) - 1)

        plt.subplot(3, 2, 1)
        plt.plot(times, self.xs)
        plt.xlabel("Time [s]")
        plt.ylabel("x [m]")

        plt.subplot(3, 2, 2)
        plt.plot(times, self.vs)
        plt.xlabel("Time [s]")
        plt.ylabel("v [m/s]")

        plt.subplot(3, 2, 3)
        Eks = [m * v ** 2 / 2 for v in self.vs]
        plt.plot(times, Eks)
        plt.xlabel("Time [s]")
        plt.ylabel("Ek [J]")

        Vs = [self.V(x) for x in self.xs]
        plt.subplot(3, 2, 4)
        plt.plot(times, Vs)
        plt.xlabel("Time [s]")
        plt.ylabel("V [m/s]")

        plt.subplot(3, 2, 5)
        Ls = [Vs[i] + Eks[i] for i in range(len(Vs))]
        plt.plot(times, Ls)

        plt.show()
    def ph(self, stop_time):
        times = np.linspace(0, stop_time, int(stop_time / self.delT))
        self.precalc(len(times) - 1)
        plt.plot(self.xs[:len(times)-1], self.vs[:len(times)-1])


# ex 1
delX = 0.001
exp_duration = 50 # s
# calculate appropiate number of x's and v's
experiment1 = experiment(m, delX, 0.1)
experiment2 = experiment(m, delX, 0.01)

experiment1.plot(exp_duration)
experiment2.plot(exp_duration)

plt.subplot(2,2,1)
experiment1.ph(100)
plt.subplot(2,2,2)
experiment1.ph(1000)
plt.subplot(2,2,3)
experiment2.ph(100)
plt.subplot(2,2,4)
experiment2.ph(1000)
plt.show()