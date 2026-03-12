#!/usr/bin/env python3
import numpy as np
from matplotlib import pyplot as plt
from abc import ABC, abstractmethod
from scipy.optimize import fsolve

m = 1
delX = 0.001

class experimentBase(ABC):
    def v(self, n):
        if len(self.vs) >= n + 1:
            return self.vs[n]

        return self.calc_v_till(n)

    @abstractmethod
    def calc_v_till(self, n):
        """
        This should handle situations when suitable v not in vs vector. It should:
        - do its math
        - append to self.vs
        - return the result
        :param n: index till wihich it should calculate
        :return:
        """
        pass

    def x(self, n):
        if len(self.xs) >= n+1:
            return self.xs[n]
        return self.calc_x_till(n)

    @abstractmethod
    def calc_x_till(self, n):
        """
        This should handle situations when suitable x not in xs vector. It should:
        - do its math
        - append to self.xs
        - return the result
        :param n:
        :return:
        """
        pass

    def __init__(self, m, delX, delT):
        self.vs = [0]
        self.xs = [0]
        self.m = m
        self.delX = delX
        self.delT = delT

    def V(self, x):
        return -np.exp(-((x+2)/3)**2) - 1.2 * np.exp(-((x-2)/3)**2)# * si.J

    def devV(self, xn):
        return (self.V(xn + self.delX) - self.V(xn))/(2*self.delX)

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

class experiment(experimentBase):
    def __init__(self, m, delX, delT, alfa=0):
        super().__init__(m, delX, delT)
        self.alfa = alfa
    def calc_v_till(self, n):
       while len(self.vs) < n+1:
           last = len(self.vs) - 1
           self.vs.append(self.vs[last] - 1/m * self.devV(self.x(last))*self.delT - self.alfa * self.vs[last] * self.delT)
       return self.v(n)

    def calc_x_till(self,n):
        while len(self.xs) < n+1:
            last = len(self.xs) - 1
            self.xs.append(self.xs[last] + self.v(last)*self.delT)
        return self.x(n)

# ex 1
def ex1():
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

# Ex2
def ex2():
    alpha = [.5, 5, 201] # need more!
    for a in alpha:
        experiment1 = experiment(m, delX, 0.01, a)
        plt.title = f"alpha = {a}"
        experiment1.plot(50)

    i = 1
    for a in alpha:
        experiment2 = experiment(m, delX, 0.01, a)
        plt.subplot(2,2,i)
        experiment2.ph(50)
        i += 1
    plt.show()

# This class is an intro for ex3 and ex4 (ex2 is here because I missunderstood the instructions)
class experimentEx2(experimentBase):
    def __init__(self, m, delX, delT, alfa):
        super().__init__(m, delX, delT)
        self.alfa = alfa
        self.insights = [None]
    # arbitrary equations writen in form of f(x) = 0 for fsolve
    def xn_plus_1(self,x, xn, vn_plus_1, vn):
        return x - xn - self.delT/2 * (vn_plus_1 + vn)
    def vn_plus_1(self, v, vn, xn_plus_1, xn):
        return v - vn - self.delT/2 * (-1/m * self.devV(xn_plus_1) - self.alfa * v - 1/m * self.devV(xn) - self.alfa * vn)

    def equations(self, args, last):
        xn_plus_1, vn_plus_1 = args
        f1 = self.xn_plus_1(xn_plus_1, self.xs[last], vn_plus_1, self.vs[last])
        f2 = self.vn_plus_1(vn_plus_1, self.vs[last], xn_plus_1, self.xs[last])
        return [f1, f2]

    def calc_x_till(self, n):
        while len(self.xs) < n+1:
            last = len(self.xs) - 1
            roots, info, ier, msg = fsolve(lambda x : self.equations(x, last), [0,0], full_output=True)
            self.insights.append(info)
            self.xs.append(roots[0])
            self.vs.append(roots[1])
        return self.xs[-1]
    def calc_v_till(self, n):
        self.calc_x_till(n)
        return self.vs[n]


# Ex3
def ex3():
    experimentEx3 = experimentEx2(m, delX, 0.01, 0)
    experimentEx3.precalc(1)
    info = experimentEx3.insights[1] # for the first timestamp
    n_iter = info['nfev']
    accuracy = experimentEx3.equations([experimentEx3.x(1), experimentEx3.v(1)], 0)
    print(f"Number of iterations: {n_iter}. Final accuracy: {accuracy[0]:.20f} for x and {accuracy[1]:.20f} for v")

# ex4
def ex4():
    exp_duration = 50 # s
    # calculate appropiate number of x's and v's
    exp1 = experimentEx2(m, delX, 0.1, 0)
    exp2 = experimentEx2(m, delX, 0.01, 0)

    exp1.plot(exp_duration)
    exp2.plot(exp_duration)

    plt.subplot(2,2,1)
    exp1.ph(100)
    plt.subplot(2,2,2)
    exp1.ph(1000)
    plt.subplot(2,2,3)
    exp2.ph(100)
    plt.subplot(2,2,4)
    exp2.ph(1000)
    plt.show()

    alpha = [.5, 5, 201] # need more!
    for a in alpha:
        ex2experiment1 = experimentEx2(m, delX, 0.01, a)
        plt.title = f"alpha = {a}"
        ex2experiment1.plot(50)

    i = 1
    for a in alpha:
        ex2experiment2 = experimentEx2(m, delX, 0.01, a)
        plt.subplot(2,2,i)
        ex2experiment2.ph(50)
        i += 1
    plt.show()


print("Ex1")
ex1()
print("Ex2")
ex2()
print("Ex3")
ex3()
print("Ex4")
ex4()