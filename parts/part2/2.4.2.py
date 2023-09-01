"""
本程序为2.4.2节的相关内容
"""

import numpy as np
from matplotlib import pyplot as plt

print("# ----------(1)----------")
x = np.linspace(0, np.pi * 2)  # 创建x数组
y = np.sin(x)  # 创建y数组
plt.plot(x, y, color="#0080ff")  # 用plot绘制曲(折)线图
plt.xlim(xmin=0, xmax=2 * np.pi)  # 设置横坐标
plt.ylim([-1.2, 1.2])  # 设置纵坐标
plt.xlabel("x")  # 设置横轴
plt.ylabel("sin(x)")  # 设置纵轴
plt.title("y=sin(x)")
plt.show()  # 显示图像

print("# ----------(2)----------")
x = np.random.normal(50, 20, 500)
# 生成均值为50，标准差为20的500个正态分布数据
plt.xlim(xmin=-20, xmax=120)  # 设置x坐标为[-20,100]
plt.ylim([0, 100])  # 设置y坐标为[0,100]
plt.xlabel("value")
plt.ylabel("rate")
plt.hist(x, bins=15, edgecolor="#000000", color="#0080FF", histtype="bar", alpha=0.5)  # 绘制图像
plt.show()  # 显示图像

print("# ----------(3)----------")
ax = plt.axes(projection='3d')  # 生成3d画布
x = np.arange(-5, 5, 0.5)  # 设置自变量x的值
y = np.arange(-5, 5, 0.5)  # 设置自变量y的值
x, y = np.meshgrid(x, y)  # 网格化
z = -x ** 2 - y ** 2  # 计算
ax.scatter3D(x, y, z, color="#0080FF")  # 绘制3d图像
plt.show()

ax = plt.axes(projection='3d')  # 生成3d画布
x = np.arange(-5, 5, 0.5)  # 设置自变量x的值
y = np.arange(-5, 5, 0.5)  # 设置自变量y的值
x, y = np.meshgrid(x, y)  # 网格化
z = -x ** 2 - y ** 2  # 计算
ax.plot_surface(x, y, z, cmap='viridis')  # 绘制3d图像，cmap设置图像颜色
plt.show()
