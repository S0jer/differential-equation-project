# MES - ODKSZTALCENIE SPREZYSTE
from math import sqrt
import matplotlib.pyplot as plt
from numpy import array
from scipy.linalg import lu

# Poczatkowe zmienne
n = 100
start = 0
end = 2


# Funkcja E
def E(x):
    if x <= 1.0:
        return 3.0
    return 5.0


# I-ty punkt podzialu
def x_i(i):
    return i * (end - start) / n


# Funkcje bazowe
def e(x, i):
    if x > x_i(i - 1) and x <= x_i(i):
        return n / end * (x - x_i(i - 1))
    elif x > x_i(i) and x < x_i(i + 1):
        return n / end * (x_i(i + 1) - x)
    return 0


# Pochodna funkcji bazowej
def e_dx(x, i):
    if x > x_i(i - 1) and x <= x_i(i):
        return n / end
    elif x > x_i(i) and x < x_i(i + 1):
        return -n / end
    return 0


# Funkcja calkuje kwadratura Gaussa-Legendrea (2 punkty kwadratury)
def integrate(f, a, b):
    middle = (b - a) / 2
    addEnds = (b + a) / 2
    gaussPoint = 1 / sqrt(3)
    return (
            middle * (f(middle * gaussPoint + addEnds) + f(middle * (-1) * gaussPoint + addEnds))
    )


# Lewa strona sformulowania wariacyjnego
def B(i, j):
    f = max(0, x_i(i - 1), x_i(j - 1))
    t = min(x_i(i + 1), x_i(j + 1), end)
    return integrate(lambda x: (E(x) * e_dx(x, i) * e_dx(x, j)), f, t) - E(0) * e(0, i) * e(0, j)


# Prawa strona sformulowania wariacyjnego
def L(i):
    return -10 * E(0) * e(0, i)


# Funkcja ktora rozwiazuje nasze rownanie
def solution():
    # Macierz glowna B
    M = array([[0 for _ in range(n)] for _ in range(n)])
    for i in range(n):
        for j in range(n):
            if abs(i - j) <= 1:
                M[i, j] = B(i - 1, j - 1)

    # Prawa strona (kolumna)
    V = []
    for i in range(n):
        V.append(L(i - 1))

    # Rozwiazanie
    U = lu(M, V)
    return U


# Funkcja oblicza kombinacje liniowa funkcji bazowych
def linearCombination(x):
    result = 0
    for i in range(n):
        result += U[i] * e(x, i - 1)

    return result


# Funkcja rysujaca rozwiazanie rownania
def drawSolution(U):
    print(U)
    # plt.plot(U)
    plt.show()


U = solution()
drawSolution(U)
