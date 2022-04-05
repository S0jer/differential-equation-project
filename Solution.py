# MES - Problem odkształcenia sprężystego

import math
from numpy import array
import matplotlib.pyplot as plot
from scipy.linalg import solve


# Funkcja E:
def E(x):
    if x <= 1:
        return 3
    return 5


# i-ty punkt podzialu
def x_i(i):
    return i * (end - start) / discretization


# Funkcje bazowe
def basisFunctions(x, i):
    if x_i(i - 1) < x <= x_i(i):
        return discretization / end * (x - x_i(i - 1))
    elif x_i(i) < x < x_i(i + 1):
        return discretization / end * (x_i(i + 1) - x)
    return 0


# Pochodna funkcji bazowej
def dBasisFunctionsDx(x, i):
    if x_i(i - 1) < x <= x_i(i):
        return discretization / end
    elif x_i(i) < x < x_i(i + 1):
        return -discretization / end
    return 0


# Funkcja calkuje kwadratura Gaussa-Legendrea
def integrate(f, a, b):
    middle = (b - a) / 2
    addEnds = (b + a) / 2
    gaussPoint = 1 / math.sqrt(3)
    return middle * (f(middle * gaussPoint + addEnds) + f(middle * (-1) * gaussPoint + addEnds))


# Funkcja oblicza kombinacje liniowa funkcji bazowych
def linearCombination(x):
    result = 0
    for i in range(discretization):
        result = result + uResult[i] * basisFunctions(x, i)
    return result


# Funkcja rysujaca wykres
def drawChart():
    arr, i, results = [start], 0, []

    while i <= end - 1 / discretization:
        i += 1 / discretization
        arr.append(i)
    if end not in arr:
        arr.append(end)
    for i in arr:
        results.append(linearCombination(i))

    plot.plot(arr, results)
    plot.show()


# Rozwiaznie naszego rownania
if __name__ == '__main__':

    # Dane wejściowe:
    discretization, start, end = 100, 0, 2

    # Macierz glowna bMatrix
    bMatrix = array([[0.0 for _ in range(discretization)] for _ in range(discretization)])

    for i in range(discretization):
        for j in range(discretization):
            if abs(i - j) <= 1:
                # Lewa strona sformulowania wariacyjnego:
                fm = max(0, x_i(i - 1), x_i(j - 1))
                t = min(x_i(i + 1), x_i(j + 1), end)

                bMatrix[i][j] = integrate(lambda x: (E(x) * dBasisFunctionsDx(x, i) * dBasisFunctionsDx(x, j)), fm, t) \
                                - E(0) * basisFunctions(0, i) * basisFunctions(0, j)

    # Prawa strona (kolumna)
    lMatrix = array([0 for _ in range(discretization)])

    for i in range(discretization):
        # Prawa strona sformulowania wariacyjnego: (-10) * E(0) * basisFunctions(0, i)
        lMatrix[i] = (-10) * E(0) * basisFunctions(0, i)

    # Rozwiazanie
    uResult = solve(bMatrix, lMatrix)

    drawChart()
