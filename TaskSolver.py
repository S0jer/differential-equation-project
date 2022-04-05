from math import sqrt
from xml import dom
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.integrate.odepack import odeint
from scipy.linalg import lu
from numpy import array

E = lambda x: 3.0 if x <= 1.0 else 5.0

n = 100
start = 0
end = 2


def basis(discretization, dom, i, x):
    h = dom / discretization
    hInv = discretization / dom

    center = dom * i / discretization
    left = center - h
    right = center + h

    if x < left or x > right:
        return 0.0
    if x <= center:
        return (x - left) * hInv
    return (right - x) * hInv


def dBasis_dx(discretization, dom, i, x):
    h = dom / discretization
    hInv = discretization / dom

    center = dom * i / discretization
    left = center - h
    right = center + h

    if x < left or x > right:
        return 0.0
    if x <= center:
        return hInv
    return -1*hInv


def L(v):  # , start, end
    return -10*E(0)*v

def drawSolution(U):
    arr = []
    i = 0
    while i < end:
        i += 1 / n
        arr.append(i)
    results = []
    for i in arr:
        results.append(linearCombination(i))
    plt.plot(results)
    plt.show()

def x_i(i):
    return i * (end - start) / n



def linearCombination(x):
    result = 0
    for i in range(n):
        result += U[i] * basis(discretization, dom, i, x)

    return result



if __name__ == '__main__':
    N = int(input("Set N: "))

    dom = 2
    discretization = N

    bMatrix = array([[0.0 for _ in range(discretization)] for _ in range(discretization)])

    for i in range(discretization):
        for j in range(discretization):
            integral = 0
            if abs(i - j) <= 1:
                integrateFrom = dom * max(max(i, j) - 1, 0) / discretization
                integrateTo = dom * min(min(i, j) + 1, discretization) / discretization

                integral = quad(lambda x: E(x) * dBasis_dx(discretization, dom, i, x) \
                                          * dBasis_dx(discretization, dom, j, x), integrateFrom, integrateTo)[0]

            bMatrix[i][j] = -E(0) * basis(discretization, dom, i, 0) * basis(discretization, dom, j, 0) + integral

        #right = np.zeros(discretization + 1)
        right = []

        for i in range(discretization):
            right.append(L(basis(discretization, dom, i, 0)))
            #right[i] = L(basis(discretization, dom, i)(0))  # start, end

        #right[discretization] = 0

        U = lu(bMatrix, right)
        print(U)
        print(bMatrix)

        print("right", right)

    drawSolution(U)