"""
本程序是为了统一绘图风格所自定义的绘图库函数
"""
from matplotlib import pyplot as plt
import numpy as np
from matplotlib.pyplot import MultipleLocator
import matplotlib as mpl


def wave_figure(data, comment, x_tick=5, y_tick=0.1, save_file="Figure1.svg"):
    """
    本函数是专门绘制地震波相关图像的函数，可以直接接受read_quake_wave的返回值
    绘图参数如下：
    图框长14cm，宽3.5cm，图框线宽1.5pt。
    纵坐标采用七度标注法，横坐标采用5s标注法，刻度线左边和下边有，线宽1pt。
    Parameters
    ----------
    data 地震波数据

    Returns 直接显示图像
    -------

    """
    fig = plt.figure(figsize=(5.5118, 1.9685))  # 设置图像边框
    fig.patch.set_facecolor('white')  # 设置背景色
    bwith = 1  # 边框宽度为1.5磅
    ax = plt.gca()  # 获取边框
    ax.spines['bottom'].set_linewidth(bwith)  # 图框下边
    ax.spines['left'].set_linewidth(bwith)  # 图框左边
    ax.spines['top'].set_linewidth(bwith)  # 图框上边
    ax.spines['right'].set_linewidth(bwith)  # 图框右边
    # 坐标轴
    plt.xlabel("time/s", fontdict={'family': 'Times New Roman', 'size': 10.5}, labelpad=3)  # labelpad为标题距离刻度线范围
    plt.ylabel(comment, fontdict={'family': 'Times New Roman', 'size': 10.5}, labelpad=3)
    # 刻度标签
    plt.xticks(fontproperties='Times New Roman', size=10.5)
    plt.yticks(fontproperties='Times New Roman', size=10.5)
    plt.tick_params(axis="both", which="major", direction="in", width=1, length=5, pad=5)
    # 刻度范围
    bandwidth = 2 * y_tick

    plt.xlim(0, 45)  # x轴范围设置
    plt.ylim(-bandwidth, bandwidth)  # y轴范围设置
    x_major_locator = MultipleLocator(x_tick)  # x轴刻度线间隔
    y_major_locator = MultipleLocator(y_tick)  # x轴刻度线间隔
    ax.xaxis.set_major_locator(x_major_locator)
    ax.yaxis.set_major_locator(y_major_locator)
    delta_time = data[1]
    x = np.arange(0, len(data[0])) * delta_time
    plt.plot(x, data[0], color="#0080FF")
    plt.gcf().subplots_adjust(left=0.100, right=0.970, top=0.900, bottom=0.220)
    plt.savefig("../../res/" + save_file)
    plt.show()


def response_figure(data, comment, x_tick=5, y_tick=0.005, x_length=48, delta_time=0.02, save_file="Figure1.svg"):
    """
    本函数支持任意条数的地震响应绘制，但是要求响应时长，响应峰值尽量相近
    否则画出来的图像并不美观。
    Parameters
    ----------
    data 地震响应数组
    comment 注释数组，与data相对应[标签,颜色,线性]
    x_tick 横坐标间隔
    y_tick 纵坐标间隔
    save_file 文件保存位置

    Returns
    -------

    """
    mpl.rcParams['font.family'] = 'serif'
    mpl.rcParams['font.serif'] = 'Times New Roman'
    mpl.rcParams['font.style'] = 'normal'
    mpl.rcParams['font.variant'] = 'normal'
    mpl.rcParams['font.weight'] = 'normal'
    mpl.rcParams['font.stretch'] = 'normal'
    mpl.rcParams['font.size'] = 10.5
    fig = plt.figure(figsize=(5.5118, 1.9685))  # 设置图像边框
    fig.patch.set_facecolor('white')  # 设置背景色
    bwith = 1  # 边框宽度为1.5磅
    ax = plt.gca()  # 获取边框
    ax.spines['bottom'].set_linewidth(bwith)  # 图框下边
    ax.spines['left'].set_linewidth(bwith)  # 图框左边
    ax.spines['top'].set_linewidth(bwith)  # 图框上边
    ax.spines['right'].set_linewidth(bwith)  # 图框右边
    # 坐标轴
    plt.xlabel("time/s", fontdict={'family': 'Times New Roman', 'size': 10.5}, labelpad=3)  # labelpad为标题距离刻度线范围
    plt.ylabel("displacement/m", fontdict={'family': 'Times New Roman', 'size': 10.5}, labelpad=3)
    # 刻度标签
    plt.xticks(fontproperties='Times New Roman', size=10.5)
    plt.yticks(fontproperties='Times New Roman', size=10.5)
    plt.tick_params(axis="both", which="major", direction="in", width=1, length=5, pad=5)
    # 刻度范围
    bandwidth = 2 * y_tick

    plt.xlim(0, x_length)  # x轴范围设置
    plt.ylim(-bandwidth, bandwidth)  # y轴范围设置
    x_major_locator = MultipleLocator(x_tick)  # x轴刻度线间隔
    y_major_locator = MultipleLocator(y_tick)  # x轴刻度线间隔
    plt.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))

    ax.xaxis.set_major_locator(x_major_locator)
    ax.yaxis.set_major_locator(y_major_locator)
    for i in range(len(data)):
        x = np.linspace(0, len(data[i]) * delta_time, len(data[i]))
        plt.plot(x, data[i], label=comment[i][0], color=comment[i][1], linestyle=comment[i][2])
    plt.gcf().subplots_adjust(left=0.100, right=0.970, top=0.900, bottom=0.220)
    plt.legend(loc='upper right', prop={'family': 'Times New Roman', 'size': 10.5})
    plt.savefig("../../res/" + save_file)
    plt.show()


def spectrum_figure(data, comment, y_length=1, save_file="Figure1.svg"):
    mpl.rcParams['font.family'] = 'serif'
    mpl.rcParams['font.serif'] = 'Times New Roman'
    mpl.rcParams['font.style'] = 'normal'
    mpl.rcParams['font.variant'] = 'normal'
    mpl.rcParams['font.weight'] = 'normal'
    mpl.rcParams['font.stretch'] = 'normal'
    mpl.rcParams['font.size'] = 10.5
    x = np.arange(0.01, 1, 0.01)
    x = np.hstack((x, np.arange(1.05, 2, 0.05)))
    x = np.hstack((x, np.arange(2.1, 3, 0.1)))
    x = np.hstack((x, np.arange(3.2, 4, 0.2)))
    x = np.hstack((x, np.arange(4.5, 6, 0.5)))
    fig = plt.figure(figsize=(3.1496, 2.3622))  # 设置图像边框
    fig.patch.set_facecolor('white')  # 设置背景色
    bwith = 1  # 边框宽度为1.5磅
    ax = plt.gca()  # 获取边框
    ax.spines['bottom'].set_linewidth(bwith)  # 图框下边
    ax.spines['left'].set_linewidth(bwith)  # 图框左边
    ax.spines['top'].set_linewidth(bwith)  # 图框上边
    ax.spines['right'].set_linewidth(bwith)  # 图框右边
    # 坐标轴
    plt.xlabel("Ts/s", fontdict={'family': 'Times New Roman', 'size': 10.5}, labelpad=3)  # labelpad为标题距离刻度线范围
    plt.ylabel(comment, fontdict={'family': 'Times New Roman', 'size': 10.5}, labelpad=3)
    plt.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    # 刻度标签
    plt.xticks(fontproperties='Times New Roman', size=10.5)
    plt.yticks(fontproperties='Times New Roman', size=10.5)
    plt.tick_params(axis="both", which="major", direction="in", width=1, length=5, pad=5)

    plt.xlim(0, 6)  # x轴范围设置
    plt.ylim(0, y_length)  # y轴范围设置
    x_major_locator = MultipleLocator(1)  # x轴刻度线间隔
    y_major_locator = MultipleLocator(y_length / 5)  # x轴刻度线间隔
    ax.xaxis.set_major_locator(x_major_locator)
    ax.yaxis.set_major_locator(y_major_locator)
    plt.plot(x, data, color="#0080FF")
    plt.gcf().subplots_adjust(left=0.160, right=0.970, top=0.900, bottom=0.220)
    plt.savefig("../../res/" + save_file)
    plt.show()


def hysteresis_figure(data, x_tick=1, y_tick=8, save_file="Figure1.svg"):
    mpl.rcParams['font.family'] = 'serif'
    mpl.rcParams['font.serif'] = 'Times New Roman'
    mpl.rcParams['font.style'] = 'normal'
    mpl.rcParams['font.variant'] = 'normal'
    mpl.rcParams['font.weight'] = 'normal'
    mpl.rcParams['font.stretch'] = 'normal'
    mpl.rcParams['font.size'] = 10.5
    fig = plt.figure(figsize=(3.1496, 2.3622))  # 设置图像边框
    fig.patch.set_facecolor('white')  # 设置背景色
    bwith = 1  # 边框宽度为1.5磅
    ax = plt.gca()  # 获取边框
    ax.spines['bottom'].set_linewidth(bwith)  # 图框下边
    ax.spines['left'].set_linewidth(bwith)  # 图框左边
    ax.spines['top'].set_linewidth(bwith)  # 图框上边
    ax.spines['right'].set_linewidth(bwith)  # 图框右边
    # 坐标轴
    plt.xlabel("displacement/m", fontdict={'family': 'Times New Roman', 'size': 10.5}, labelpad=3)  # labelpad为标题距离刻度线范围
    plt.ylabel("force/N", fontdict={'family': 'Times New Roman', 'size': 10.5}, labelpad=3)
    plt.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    # 刻度标签
    plt.xticks(fontproperties='Times New Roman', size=10.5)
    plt.yticks(fontproperties='Times New Roman', size=10.5)
    plt.tick_params(axis="both", which="major", direction="in", width=1, length=5, pad=5)

    x_length = 3 * x_tick
    y_length = 3 * y_tick

    plt.xlim(-x_length, x_length)  # x轴范围设置
    plt.ylim(-y_length, y_length)  # y轴范围设置
    x_major_locator = MultipleLocator(x_tick)  # x轴刻度线间隔
    y_major_locator = MultipleLocator(y_tick)  # x轴刻度线间隔
    ax.xaxis.set_major_locator(x_major_locator)
    ax.yaxis.set_major_locator(y_major_locator)
    plt.plot(data[0], data[1], color="#0080FF")
    plt.gcf().subplots_adjust(left=0.200, right=0.950, top=0.900, bottom=0.220)
    plt.savefig("../../res/" + save_file)
    plt.show()


def common(data, comment, x_tick=400, y_tick=1, save_file="Figure1.svg"):
    mpl.rcParams['font.family'] = 'serif'
    mpl.rcParams['font.serif'] = 'Times New Roman'
    mpl.rcParams['font.style'] = 'normal'
    mpl.rcParams['font.variant'] = 'normal'
    mpl.rcParams['font.weight'] = 'normal'
    mpl.rcParams['font.stretch'] = 'normal'
    mpl.rcParams['font.size'] = 10.5
    fig = plt.figure(figsize=(3.1496, 2.3622))  # 设置图像边框
    fig.patch.set_facecolor('white')  # 设置背景色
    bwith = 1  # 边框宽度为1.5磅
    ax = plt.gca()  # 获取边框
    ax.spines['bottom'].set_linewidth(bwith)  # 图框下边
    ax.spines['left'].set_linewidth(bwith)  # 图框左边
    ax.spines['top'].set_linewidth(bwith)  # 图框上边
    ax.spines['right'].set_linewidth(bwith)  # 图框右边
    # 坐标轴
    plt.xlabel(comment[0], fontdict={'family': 'Times New Roman', 'size': 10.5}, labelpad=3)  # labelpad为标题距离刻度线范围
    plt.ylabel(comment[1], fontdict={'family': 'Times New Roman', 'size': 10.5}, labelpad=3)
    plt.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    # 刻度标签
    plt.xticks(fontproperties='Times New Roman', size=10.5)
    plt.yticks(fontproperties='Times New Roman', size=10.5)
    plt.tick_params(axis="both", which="major", direction="in", width=1, length=5, pad=5)

    x_length = 5 * x_tick
    y_length = 3 * y_tick

    plt.xlim(0, x_length)  # x轴范围设置
    plt.ylim(-y_length, y_length)  # y轴范围设置
    x_major_locator = MultipleLocator(x_tick)  # x轴刻度线间隔
    y_major_locator = MultipleLocator(y_tick)  # x轴刻度线间隔
    ax.xaxis.set_major_locator(x_major_locator)
    ax.yaxis.set_major_locator(y_major_locator)
    plt.plot(data[0], data[1], color="#0080FF")
    plt.gcf().subplots_adjust(left=0.140, right=0.945, top=0.900, bottom=0.220)
    plt.savefig("../../res/" + save_file)
    plt.show()
