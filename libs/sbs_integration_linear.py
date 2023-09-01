"""
线性逐步积分法计算程序
"""
import numpy as np
from scipy import integrate
import scipy as sc
from multiprocessing import Pool
from scipy.linalg import *
from libs.damping import rayleigh


def duhamel_parse(mass, stiffness, load, delta_time, damping_ratio=0.05, result_length=4000):
    """
    Duhamel积分的计算程序
    选取的应用场景为单自由度有阻尼体系在具有解析表达的外荷载影响下的振动
    注意：本程序只支持解析荷载计算
    Parameters
    ----------
    load 荷载
    mass 质量
    stiffness 刚度
    damping_ratio 阻尼比

    Returns 位移
    -------

    """
    omega_n = (stiffness / mass) ** 0.5  # 无阻尼自由振动周期
    omega_d = omega_n * (1 - damping_ratio ** 2) ** 0.5  # 有阻尼自由振动周期
    dpm = np.zeros(int(result_length))

    def func(tau):
        """
        duhamel被积函数
        """
        return load(tau) * np.e ** (-damping_ratio * omega_n * (t - tau)) * np.sin(omega_d * (t - tau))

    for i in range(len(dpm)):
        t = i * delta_time
        dpm[i] = 1 / omega_d * integrate.quad(func, 0, t)[0]

    return dpm


def duhamel_numerical(mass, stiffness, load, delta_time, damping_ratio=0.05, result_length=0):
    """
    Duhamel积分数值计算程序
    """
    if result_length == 0:
        result_length = int(1.2 * len(load))
    load = np.pad(load, (0, result_length - len(load)))  # 末端补零
    omega_n = (stiffness / mass) ** 0.5  # 无阻尼自由振动周期
    omega_d = omega_n * (1 - damping_ratio ** 2) ** 0.5  # 有阻尼自由振动周期
    dpm = np.zeros(result_length)  # 初始化位移
    for i in range(result_length):
        t = i * delta_time
        dpm_temp = 0
        for j in range(i):
            tau = j * delta_time
            dpm_temp += load[j] * np.e ** (-damping_ratio * omega_n * (t - tau)) * np.sin(
                omega_d * (t - tau)) * delta_time
        dpm[i] = dpm_temp / (omega_d * mass)
    return dpm


def func(i, delta_time, load, damping_ratio, omega_n, t, omega_d):
    dpm_temp = 0
    for j in range(i):
        tau = j * delta_time
        dpm_temp += load[j] * np.e ** (-damping_ratio * omega_n * (t - tau)) * np.sin(
            omega_d * (t - tau)) * delta_time
    return dpm_temp


def duhamel_numerical_multprocess(mass, stiffness, load, delta_time, damping_ratio=0.05):
    """
    Duhamel积分的计算程序
    Parameters
    ----------
    delta_time 积分间隔
    mass 质量
    stiffness 刚度
    load 荷载
    omega 频率
    damping_ratio 阻尼比

    Returns 结构响应
    -------

    """
    length = len(load)
    omega_n = (stiffness / mass) ** 0.5  # 无阻尼自由振动周期
    omega_d = omega_n * (1 - damping_ratio ** 2) ** 0.5  # 有阻尼自由振动周期
    dpm = np.zeros(length)  # 初始化位移
    pool = Pool(8)  # 创建八个线程
    process_list = []
    for i in range(length):
        t = i * delta_time
        process_list.append(
            pool.apply_async(func=func, args=(i, delta_time, load, damping_ratio, omega_n, t, omega_d,)))
    for i in range(length):
        dpm_temp = process_list[i].get()
        dpm[i] = dpm_temp / (omega_d * mass)
        print("\rcalculated time:%.3f dpm = %.6f" % (i * delta_time, dpm[i]), end="")
    pool.close()
    pool.join()
    return dpm


def segmented_parsing(mass, stiffness, load, delta_time,
                      damping_ratio=0.05, dpm_0=0, vel_0=0,
                      result_length=0):
    """
    分段解析法计算程序，分段解析法一般适用于单自由度体系的动力响应求解，
    所以仅考虑单自由度情况下的线性分段解析法。
    Parameters
    ----------
    load 荷载;一维列表
    delta_time 时间步长;浮点数
    mass 质量;浮点数
    stiffness 刚度;浮点数
    damping_ratio 阻尼比;浮点数
    dpm_0 初始位移;浮点数
    vel_0 初始速度;浮点数
    result_length 结果长度;整数

    Returns 位移，速度;二维数组，二维数组
    -------

    """
    # 前期数据准备
    # 为了方便代码阅读和减少重复参数所进行的参数代换
    omega_n = np.sqrt(stiffness / mass)
    omega_d = omega_n * np.sqrt(1 - damping_ratio ** 2)
    temp_1 = sc.e ** (-damping_ratio * omega_n * delta_time)
    temp_2 = damping_ratio / np.sqrt(1 - damping_ratio ** 2)
    temp_3 = 2 * damping_ratio / (omega_n * delta_time)
    temp_4 = (1 - 2 * damping_ratio ** 2) / (omega_d * delta_time)
    temp_5 = omega_n / np.sqrt(1 - damping_ratio ** 2)
    sin = np.sin(omega_d * delta_time)
    cos = np.cos(omega_d * delta_time)

    # 计算所需参数
    A = temp_1 * (temp_2 * sin + cos)
    B = temp_1 * (sin / omega_d)
    C = 1 / stiffness * (temp_3 + temp_1 * (
            (temp_4 - temp_2) * sin - (1 + temp_3) * cos
    ))
    D = 1 / stiffness * (1 - temp_3 + temp_1 * (
            -temp_4 * sin + temp_3 * cos
    ))
    A_prime = -temp_1 * (temp_5 * sin)
    B_prime = temp_1 * (cos - temp_2 * sin)
    C_prime = 1 / stiffness * (-1 / delta_time + temp_1 * (
            (temp_5 + temp_2 / delta_time) * sin + 1 / delta_time * cos
    ))
    D_prime = 1 / (stiffness * delta_time) * (
            1 - temp_1 * (temp_2 * sin + cos)
    )

    # 处理荷载长度
    if result_length == 0:
        result_length = int(1.2 * len(load))
    load = np.pad(load, (0, result_length - len(load)))  # 荷载数据末端补零

    # 初始化位移数组与速度数组
    dpm = np.zeros(result_length)
    vel = np.zeros(result_length)
    acc = np.zeros(result_length)
    dpm[0] = dpm_0
    vel[0] = vel_0

    # 正式开始迭代
    for i in range(result_length - 1):
        dpm[i + 1] = A * dpm[i] + B * vel[i] + C * load[i] + D * load[i + 1]
        vel[i + 1] = A_prime * dpm[i] + B_prime * vel[i] + C_prime * load[i] + D_prime * load[i + 1]
        acc[i + 1] = -2 * damping_ratio * omega_n * vel[i + 1] - stiffness / mass * dpm[i + 1]

    return dpm, vel, acc


def center_difference_single(mass, stiffness, load, delta_time,
                             damping_ratio=0.05, dpm_0=0, vel_0=0,
                             result_length=0):
    """
    中心差分法计算函数。本程序为单自由度计算程序
    Parameters
    ----------
    load 荷载;一维列表
    delta_time 时间步长;浮点数
    mass 质量;浮点数
    stiffness 质量;浮点数
    damping_ratio 阻尼比;浮点数
    dpm_0 初始位移;浮点数
    vel_0 初始速度;浮点数
    result_length 结果长度;整数

    Returns 位移，速度，加速度;一维数组，一维数组，一维数组
    -------

    """

    # 起步条件计算
    damping = damping_ratio * 2 * np.sqrt(stiffness * mass)
    acc_0 = 1 / mass * (load[0] - damping * vel_0 - stiffness * dpm_0)
    dpm_minus1 = dpm_0 - delta_time * vel_0 + delta_time ** 2 * acc_0 / 2

    # 前置参数计算
    equ_k = mass / delta_time ** 2 + damping / (2 * delta_time)  # 等效刚度
    a = (2 * mass) / delta_time ** 2
    b = mass / delta_time ** 2 - damping / (2 * delta_time)

    # 处理荷载长度
    if result_length == 0:
        result_length = int(1.2 * len(load))
    load = np.pad(load, (0, result_length - len(load)))  # 荷载数据末端补零

    # 初始化位移、速度、加速度
    dpm = np.zeros(result_length)
    vel = np.zeros(result_length)
    acc = np.zeros(result_length)
    dpm[0] = dpm_0
    vel[0] = vel_0
    # 迭代开始
    # 起步
    equ_p = load[0] - (stiffness - a) * dpm_0 - b * dpm_minus1
    dpm[1] = equ_p / equ_k

    for i in range(1, result_length - 1):
        # 起步完成，开启后续迭代
        equ_p = load[i] - (stiffness - a) * dpm[i] - b * dpm[i - 1]
        dpm[i + 1] = equ_p / equ_k
        vel[i] = (dpm[i + 1] - dpm[i - 1]) / (2 * delta_time)
        acc[i] = (dpm[i + 1] - 2 * dpm[i] + dpm[i - 1]) / delta_time ** 2

    return dpm, vel, acc


def center_difference_multiple(mass, stiffness, load, delta_time,
                               damping_ratio, dpm_0, vel_0,
                               result_length=0):
    """
    中心差分法计算函数，本程序为多自由度计算程序。
    Parameters
    ----------
    load 时序荷载列阵;二维矩阵
    delta_time 采样间隔;float类型
    mass 质量矩阵;二维矩阵
    stiffness 刚度矩阵;二维矩阵
    damping_ratio 阻尼比;float类型
    dpm_0 初始位移;一维数组
    vel_0 初始速度;一维数组
    result_length 结果长度;int类型

    Returns 位移，速度，加速度;数组列表，数组列表，数组列表
    -------

    """
    # 固有属性计算
    freedom = len(dpm_0.T)  # 计算自由度
    if type(damping_ratio) == np.ndarray:
        damping = damping_ratio
    else:
        damping = rayleigh(mass, stiffness, damping_ratio)  # rayleigh阻尼

    # 起步条件计算
    acc_0 = inv(mass) @ (load[0].T - damping @ vel_0.T - stiffness @ dpm_0.T)
    dpm_minus1 = dpm_0.T - delta_time * vel_0.T + delta_time ** 2 * acc_0 / 2

    # 前置参数计算
    equ_k = mass / delta_time ** 2 + damping / (2 * delta_time)  # 等效刚度矩阵
    a = (2 * mass) / delta_time ** 2
    b = mass / delta_time ** 2 - damping / (2 * delta_time)

    # 荷载长度处理
    if result_length == 0:
        result_length = int(1.2 * len(load))
    load = load + [np.array([0 for i in range(freedom)]) for i in range(result_length - len(load))]  # 荷载末端补零

    # 初始化位移、速度、加速度时程矩阵
    dpm = np.zeros((result_length, freedom))
    vel = np.zeros((result_length, freedom))
    acc = np.zeros((result_length, freedom))
    dpm[0] = dpm_0
    vel[0] = vel_0

    # 迭代开始
    # 起步
    equ_p = load[0].T - (stiffness - a) @ dpm_0.T - b @ dpm_minus1
    dpm[1] = (inv(equ_k) @ equ_p).T

    for i in range(1, result_length - 1):
        # 起步完成，开启后续迭代
        equ_p = load[i].T - (stiffness - a) @ dpm[i].T - b @ dpm[i - 1].T
        dpm[i + 1] = (inv(equ_k) @ equ_p).T
        vel[i] = ((dpm[i + 1] - dpm[i - 1]) / (2 * delta_time))
        acc[i - 1] = ((dpm[i + 1] - 2 * dpm[i] + dpm[i - 1]) / delta_time ** 2)

    return dpm, vel, acc


def newmark_beta_single(mass, stiffness, load, delta_time,
                        damping_ratio=0.05, dpm_0=0, vel_0=0,
                        acc_0=0, beta=0.25, gamma=0.5,
                        result_length=0):
    """
    单自由度体系的Newmark方法，由于Newmark方法可以通过gamma和beta数值
    的调整从而转化为平均加速度法、线性加速度法和中心差分法。所以不编写平均
    加速度法和线性加速度法这两种方法。而将中心差分法和Newmark方法分开写的
    原因是中心差分法是beta取0的Newmark方法，而Newmark方法起步中需要将beta
    作为被除数，要采取另外的一些措施。且Newmark方法是单步法，与中心差分法
    的多步法有差异，需要特殊展示多步法的相关计算过程，故将两者分开写。
    """
    # 基本数据准备和初始条件计算
    if result_length == 0:
        result_length = int(1.2 * len(load))  # 计算持时
    load = np.pad(load, (0, result_length - len(load)))
    dpm = np.zeros(result_length)
    vel = np.zeros(result_length)
    acc = np.zeros(result_length)
    dpm[0] = dpm_0
    vel[0] = vel_0
    acc[0] = acc_0
    a_0 = 1 / (beta * delta_time ** 2)
    a_1 = gamma / (beta * delta_time)
    a_2 = 1 / (beta * delta_time)
    a_3 = 1 / (2 * beta) - 1
    a_4 = gamma / beta - 1
    a_5 = delta_time / 2 * (a_4 - 1)
    a_6 = delta_time * (1 - gamma)
    a_7 = gamma * delta_time
    omega_n = np.sqrt(stiffness / mass)
    damping = 2 * mass * omega_n * damping_ratio
    equ_k = stiffness + a_0 * mass + a_1 * damping  # 计算等效刚度
    # 迭代正式开始
    for i in range(result_length - 1):
        equ_p = load[i + 1] + mass * (
                a_0 * dpm[i] + a_2 * vel[i] + a_3 * acc[i]) + damping * (
                        a_1 * dpm[i] + a_4 * vel[i] + a_5 * acc[i])  # 计算等效荷载
        dpm[i + 1] = equ_p / equ_k  # 计算位移
        acc[i + 1] = a_0 * (dpm[i + 1] - dpm[i]) - a_2 * vel[i] - a_3 * acc[i]  # 计算加速度
        vel[i + 1] = vel[i] + a_6 * acc[i] + a_7 * acc[i + 1]  # 计算速度
    return dpm, vel, acc


def newmark_beta_multiple(mass, stiffness, load, delta_time,
                          damping_ratio, dpm_0, vel_0, acc_0,
                          beta=0.25, gamma=0.5, result_length=0):
    """
    Newmark-beta法计算函数，本程序为多自由度计算程序。
    Parameters
    ----------
    mass 质量矩阵;二维矩阵
    stiffness 刚度矩阵;二维矩阵
    load 时序荷载列阵;二维列表
    delta_time 采样间隔;float类型
    damping_ratio 阻尼比;float类型
    dpm_0 0时刻位移矩阵;一维矩阵
    vel_0 0时刻速度矩阵;一维矩阵
    acc_0 0时刻加速度矩阵;一维矩阵
    beta 计算参数beta;float类型
    gamma 计算参数gamma;float类型
    result_length 计算步数;int类型

    Returns 位移，速度，加速度;数组列表，数组列表，数组列表
    -------

    """
    # 基本数据准备和初始条件计算
    freedom = len(dpm_0.T)  # 计算自由度
    damping = rayleigh(mass, stiffness, damping_ratio)  # 计算阻尼矩阵
    if result_length == 0:
        result_length = int(1.2 * len(load))  # 计算持时
    load = load + [np.array([0 for i in range(freedom)]) for i in range(result_length - len(load))]  # 末端补零
    dpm = np.zeros((result_length, freedom))
    vel = np.zeros((result_length, freedom))
    acc = np.zeros((result_length, freedom))
    dpm[0] = dpm_0
    vel[0] = vel_0
    acc[0] = acc_0
    a_0 = 1 / (beta * delta_time ** 2)
    a_1 = gamma / (beta * delta_time)
    a_2 = 1 / (beta * delta_time)
    a_3 = 1 / (2 * beta) - 1
    a_4 = gamma / beta - 1
    a_5 = delta_time / 2 * (a_4 - 1)
    a_6 = delta_time * (1 - gamma)
    a_7 = gamma * delta_time
    equ_k = stiffness + a_0 * mass + a_1 * damping  # 计算等效刚度
    # 迭代开始
    for i in range(result_length - 1):
        equ_p = load[i + 1].T + mass @ (
                a_0 * dpm[i] + a_2 * vel[i] + a_3 * acc[i]).T + damping @ (
                        a_1 * dpm[i] + a_4 * vel[i] + a_5 * acc[i]).T  # 计算等效荷载
        dpm[i + 1] = (inv(equ_k) @ equ_p).T  # 计算位移
        acc[i + 1] = a_0 * (dpm[i + 1] - dpm[i]) - a_2 * vel[i] - a_3 * acc[i]  # 计算加速度
        vel[i + 1] = vel[i] + a_6 * acc[i] + a_7 * acc[i + 1]  # 计算速度
    return dpm, vel, acc


def wilson_theta_multiple(mass, stiffness, load, delta_time,
                          damping_ratio, dpm_0, vel_0, acc_0,
                          theta=1.37, result_length=0):
    """
    wilson-theta法计算程序，本程序为多自由度计算程序。
    Parameters
    ----------
    mass 质量矩阵;二维矩阵
    stiffness 刚度矩阵;二维矩阵
    load 时序荷载列阵;二维列表
    delta_time 采样间隔;float类型
    damping_ratio 阻尼比;float类型
    dpm_0 0时刻位移矩阵;一维矩阵
    vel_0 0时刻速度矩阵;一维矩阵
    acc_0 0时刻加速度矩阵;一维矩阵
    theta 计算参数theta;flota类型
    result_length 计算步数;int类型

    Returns 位移，速度，加速度;数组列表，数组列表，数组列表
    -------

    """
    # 前置条件计算
    freedom = len(mass)  # 计算自由度
    damping = rayleigh(mass, stiffness, damping_ratio)  # 计算阻尼矩阵
    equ_stiffness = stiffness + 6 / (
            theta * delta_time) ** 2 * mass + 3 / (
                            theta * delta_time) * damping  # 计算等效刚度
    equ_stiffness_inv = inv(equ_stiffness)  # 计算逆矩阵
    # 各项参数初始化
    if result_length == 0:
        result_length = int(1.2 * len(load))  # 荷载末端补零
    load = load + [np.array([0 for i in range(freedom)]) for i in range(result_length - len(load))]
    dpm = np.zeros((result_length, freedom))
    vel = np.zeros((result_length, freedom))
    acc = np.zeros((result_length, freedom))
    dpm[0] = dpm_0
    vel[0] = vel_0
    acc[0] = acc_0
    # 开始迭代求解
    for i in range(result_length - 1):
        # 计算等效荷载
        equ_load = load[i] + theta * (load[i + 1] - load[i])
        mass_timer = 6 / (theta * delta_time) ** 2 * dpm[i] + 6 / (theta * delta_time) * vel[i] + 2 * acc[i]
        damping_timer = 3 / (theta * delta_time) * dpm[i] + 2 * vel[i] + theta * delta_time / 2 * acc[i]
        equ_load += (mass @ mass_timer.T + damping @ damping_timer.T).T
        dpm_temp = (equ_stiffness_inv @ equ_load.T).T
        # 注释掉的代码是没有采用theta=1简化的部分，两者的区别很小，所以采用了theta=1简化计算
        # acc[i + 1] = 6 / (theta ** 3 * delta_time ** 2) * (dpm_temp - dpm[i])
        # acc[i + 1] += -6 / (theta ** 2 * delta_time) * vel[i] + (1 - 3 / theta) * acc[i]
        # vel[i + 1] = (1 - 3 / (2 * theta)) * delta_time * acc[i] + (1 - 3 / theta ** 2) * vel[i]
        # vel[i + 1] += 3 / (theta ** 3 * delta_time) * (dpm_temp - dpm[i])
        # dpm[i + 1] = (delta_time ** 2 / 2 - delta_time ** 2 / (2 * theta)) * acc[i]
        # dpm[i + 1] += (delta_time - delta_time / theta ** 2) * vel[i] + (1 - 1 / theta ** 3) * dpm[i]
        # dpm[i + 1] += 1 / theta ** 3 * dpm_temp
        acc[i + 1] = 6 / (theta ** 3 * delta_time ** 2) * (dpm_temp - dpm[i])
        acc[i + 1] += -6 / (theta ** 2 * delta_time) * vel[i] + (1 - 3 / theta) * acc[i]
        vel[i + 1] = vel[i] + delta_time / 2 * (acc[i + 1] + acc[i])
        dpm[i + 1] = dpm[i] + delta_time * vel[i] + delta_time ** 2 / 6 * (acc[i + 1] + 2 * acc[i])
    return dpm, vel, acc
