import numpy as np


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
    omega_n = np.sqrt(stiffness / mass)
    damping = 2 * mass * omega_n * damping_ratio
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


def newmark_beta_single1(mass, stiffness, load, delta_time,
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
    a_3 = 1 / (2 * beta)
    a_4 = gamma / beta
    a_5 = (a_4 - 2) * delta_time / 2

    omega_n = np.sqrt(stiffness / mass)
    damping = 2 * mass * omega_n * damping_ratio
    equ_stiffness = stiffness + a_0 * mass + a_1 * damping
    # 积分步迭代开始
    for i in range(result_length - 1):
        # 计算荷载增量
        delta_load = load[i + 1] - load[i]
        # 计算等效荷载
        equ_load = delta_load + mass * (a_2 * vel[i] + a_3 * acc[i]) + damping * (a_4 * vel[i] + a_5 * acc[i])
        delta_dpm = equ_load / equ_stiffness  # 计算迭代步位移增量
        # 计算真实位移
        dpm[i + 1] = dpm[i] + delta_dpm
        # 计算真实加速度
        acc[i + 1] = a_0 * delta_dpm - a_2 * vel[i] - (a_3 - 1) * acc[i]
        # 计算速度
        vel[i + 1] = a_1 * delta_dpm + (1 - a_4) * vel[i] + (1 - a_4 / 2) * acc[i] * delta_time

    return dpm, vel, acc


if __name__ == "__main__":
    from read_wave import read_force, read_dpm
    from matplotlib import pyplot as plt

    force1 = np.array(read_force("Shear1_force.txt"))
    force2 = np.array(read_force("Shear2_force.txt"))
    dpm = np.array(read_dpm("Shear1_dpm.txt"))

    plt.plot(dpm, force1)
    plt.show()
