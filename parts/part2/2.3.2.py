"""
本程序为2.3.2节的相关内容
"""
import numpy as np

print("# ----------(1)----------")
a = np.zeros((3, 3))
print(a)

print("# ----------(2)----------")
a = np.zeros((2, 3))
print("原数组:\n", a)
a = a.reshape(3, 2)
print("第一次变换形状:\n", a)
a = a.reshape(6, -1)
print("第二次变换形状:\n", a)

print("# ----------(3)----------")
a = np.linspace(1, 9, 9)
a = a.reshape(3, 3)
print(a[1])
print(a[1][1])

print("# ----------(4)----------")
a = np.array([[-3, 3, 2], [-7, 6, -3], [1, -1, 2]])
b = np.array([1, 2, 3])
np.set_printoptions(precision=3, suppress=True)  # 设置矩阵输出格式
print("a矩阵为:\n", a)
print("a矩阵元素之和为:", sum(a))
print("a矩阵的转置为:\n", np.transpose(a))
print("a矩阵与b矩阵的乘积为:\n", np.matmul(a, b))
print("a矩阵的逆为:\n", np.linalg.inv(a))
