"""
本程序为2.1.3节的相关内容
"""
print("# ----------(1)----------")
a = 1
while a < 5:
    print(a)
    a += 1
print(a)

print("# ----------(2)----------")
for i in range(5):
    print(i)

print("# ----------(3)----------")
a = 0
while True:
    a += 1
    if a == 3:
        continue
    elif a > 5:
        break
    print(a)
