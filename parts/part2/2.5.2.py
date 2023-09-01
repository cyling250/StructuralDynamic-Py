"""
本程序为2.4.2节的相关内容
"""

import numpy as np
from scipy.linalg import *
from scipy.fftpack import fft
import matplotlib.pyplot as plt

print("# ----------(1)----------")
a = np.array([[1, 2], [3, 4]])
print("矩阵a的行列式:", det(a))
print("矩阵a的逆:\n", inv(a))

print("# ----------(2)----------")
x = np.linspace(0, 10, 1000)  # 创建横坐标
y = np.sin(2 * np.pi * 2 * x)  # 创建函数
fft_y = fft(y)  # 对函数进行快速傅里叶变换
plt.xlim([0, 10])  # 设置横坐标
plt.ylim(([-1.5, 1.5]))  # 设置纵坐标
plt.xlabel("time(s)")
plt.ylabel("amplitude")
plt.plot(x, y, color="#0080FF")  # 绘制原函数
plt.show()
plt.xlim([0, 0.01])  # 设置横坐标
plt.ylim(([0, 0.5]))  # 设置纵坐标
plt.xlabel("Fs")
plt.ylabel("amplitude")
plt.plot(x / 1000, abs(fft_y / 1000), color="#0080FF")  # 绘制fft后的图像
plt.show()
