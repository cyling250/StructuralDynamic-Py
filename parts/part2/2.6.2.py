"""
本程序为2.4.2节的相关内容
"""
import sympy as sp

print("# ----------(1)----------")
x = sp.Symbol("x")
f_x = x ** 2 + 2 * x - 1
print(sp.solve(f_x, x))

print("# ----------(2)----------")
i = sp.Symbol("i")
f_x2 = sp.summation(i, (i, 1, 10))
print(f_x2)

print("# ----------(3)----------")
x = sp.symbols("x")
f_x3 = sp.sin(x) / x
lim = sp.limit(f_x3, x, 0)
print(lim)

print("# ----------(4)----------")
x = sp.symbols("x")
f_x4 = sp.sin(x)
print(sp.diff(f_x4, x))

print("# ----------(5)----------")
x = sp.symbols("x")
f_x5 = 2 * x + 1
print(sp.integrate(f_x5, x))

print("# ----------(6)----------")
x, y = sp.symbols("x,y")
z = 2 * x ** 2 + 3 * x + 4 * y ** 2 + 5 * y + 6
print(sp.diff(z, x))
print(sp.diff(z, y))