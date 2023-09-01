"""
本程序为2.1.1节的相关内容
"""
print("# ----------(1)----------")
a = 10
b = 10.0
c = complex(10, 10)
print(a, type(a))
print(b, type(b))
print(c, type(c))

print("# ----------(2)----------")
s = "Welcome to dynpy"
print(s)
print(s[0:7])
print(s[::-1])

print("# ----------(3)----------")
l = [1, "welcome", [1, "welcome"], (1, 2), {1: "welcome"}]
print(l)
print(type(l))
print(l[2])
print(l[1:3])

print("# ----------(4)----------")
t = (1, "welcome", [1, "welcome"], (1, 2), {1: "welcome"})
print(t)
print(type(t))
print(t[1])

print("# ----------(5)----------")
d = {1: "number",
     "dynpy": "string"}
print(d)
print(type(d))
print(d["dynpy"])
