"""
程序描述：本程序封装了阻尼的相关参数和操作，
提供了计算Reyleigh阻尼和Caughey阻尼的接口。
需要注意的是，在本程序中仅适用于计算小型结构的阻尼矩阵，
大型结构计算阻尼矩阵时，应该先用其他高速方法计算出前n阶
振型和周期后代入计算。原因是使用numpy自带的eig方法计
算特征向量和特征值时，对于数据量超过1e6的问题可能出现
计算时间过长。
"""
import numpy as np
from scipy.linalg import *


def rayleigh(mass, stiffness, damping_ratio=0.05, seq=(0, 1)):
    """
    计算Reyleigh阻尼的函数，输入结构的质量矩阵、
    刚度矩阵、阻尼比和‘感兴趣’的振型序列号，输出Reyleigh阻尼
    Parameters
    ----------
    mass 质量矩阵，二维数组
    stiffness 刚度矩阵，二维数组
    damping_ratio 阻尼比，浮点数/一维数组
    seq 序列号列阵，整数元组

    Returns Reyleigh阻尼矩阵，二维数组
    -------

    """
    # 计算w^2
    omega = np.real(eig(stiffness, mass)[0])
    # 加上real()是因为eig计算结果有虚部，但该问题中理论上没有虚数，使用real()去除虚部

    omega = np.sort(omega)  # 使频率按照升序排列
    if len(omega) == 1:
        return 2 * mass * omega * damping_ratio

    # 计算两个频率
    omega_i = np.sqrt(omega[seq[0]])
    omega_j = np.sqrt(omega[seq[1]])
    omega_temp = np.array(
        [[omega_j, -omega_i],
         [-1 / omega_j, 1 / omega_i]]
    )

    # 计算Reyleigh阻尼的a,b
    if type(damping_ratio) == np.array:
        # 提供了变阻尼比的构造方式
        a = 2 * omega_i * omega_j / (omega_j ** 2 - omega_i ** 2) * omega_temp @ damping_ratio
        b = a[1, 0]
        a = a[0, 0]
    else:
        a = 2 * omega_i * omega_j / (omega_i + omega_j) * damping_ratio
        b = 2 * damping_ratio / (omega_i + omega_j)
    if damping_ratio == 0:
        a, b = 0, 0
    # 计算阻尼矩阵
    c = a * mass + b * stiffness
    return c


def caughey(mass, stiffness, damping_ratio, seq=None):
    """
    计算caughey阻尼的函数，输入结构的质量矩阵、
    刚度矩阵、阻尼比列阵、序列号列阵，输出caughey阻尼。
    需要注意的是，caughey阻尼函数必须输入阻尼比，当序列
    号列阵没有输入时，默认采用全部频率的阻尼比构造阻尼矩阵
    Parameters
    ----------
    mass 质量矩阵，二维数组
    stiffness 刚度矩阵，二维数组
    damping_ratio 阻尼比，浮点数
    seq 序列号，整数元组

    Returns caughey阻尼矩阵，二维数组
    -------

    """

    # 如果没有输入阻尼序列，默认使用全部阻尼构造阻尼矩阵
    if not seq:
        seq = np.arange(0, len(mass), 1, dtype='i')

    # 参数初始化
    freedom = len(mass)
    length = len(seq)
    c = np.zeros((freedom, freedom))  # 初始化阻尼矩阵
    damping_ratio = damping_ratio * np.ones(length)
    omega_n = np.zeros((length, length))  # 初始化频率矩阵
    omega = np.real(eig(stiffness, mass)[0])  # 计算w^2
    # 加上real()是因为eig计算结果有虚部，但该问题中理论上没有虚数，使用real()去除虚部
    omega = np.sort(omega)  # 使频率按照升序排列
    omega = np.sqrt(omega)  # 计算真实频率,开方

    # a = 2*(omega_n^-1)*zeta
    # 构造omega_n矩阵
    for i in range(length):
        for j in range(length):
            omega_n[i, j] = omega[seq[i]] ** (2 * (j + 1) - 3)
    a = 2 * inv(omega_n) @ damping_ratio

    # 构造阻尼矩阵
    for i in range(length):
        c = c + a[i] * mass @ (inv(mass) @ stiffness) ** i

    return c


def damping_nonclassical(mass,
                         stiffness,
                         shear_stiffness,
                         child_damping_ratio,
                         child_shear_num,
                         omega_i=0,
                         omega_j=0):
    """
    非比例阻尼计算程序，采用了递归的思路
    需要特别说明的是，child_damping_ratio列表存放的阻尼比应该从结构的高层向低层排序
    Parameters
    ----------
    mass 质量矩阵，二维数组
    stiffness 刚度矩阵，二维数组
    shear_stiffness 各层刚度，一维数组
    child_damping_ratio 子结构阻尼比，一维数组
    child_shear_num 子结构起始与结束层编号，一维数组
    omega_i
    omega_j

    Returns 非比例阻尼
    -------

    """
    # 如果开始没有给出频率，需要自己计算频率
    if omega_i == 0:
        omega = np.real(eig(stiffness, mass)[0])
        omega = np.sort(omega)
        omega_i = np.sqrt(omega[0])
        omega_j = np.sqrt(omega[1])

    if len(child_damping_ratio) == 1:
        # 如果阻尼比列表长度为1，说明当前的结构已经最简化
        a = 2 * omega_i * omega_j / (omega_i + omega_j) * child_damping_ratio[0]
        b = 2 * child_damping_ratio[0] / (omega_i + omega_j)
        return a * mass + b * stiffness

    else:
        # 按照child_shear_num[-1]的指示来分块
        flag = child_shear_num[-1]
        child_mass_1 = mass[:flag, :flag]
        child_mass_4 = mass[flag:, flag:]
        child_stiffness_1 = stiffness[:flag, :flag]
        child_stiffness_2 = stiffness[:flag, flag:]
        child_stiffness_4 = stiffness[flag:, flag:]

        # 迭代
        stiffness_r = np.zeros((len(child_mass_1), len(child_mass_1)))
        stiffness_r[-1, -1] = shear_stiffness[flag]  # 提取出K_r
        damping_1 = damping_nonclassical(child_mass_1,
                                         child_stiffness_1 - stiffness_r,
                                         shear_stiffness,
                                         child_damping_ratio[:-1],
                                         child_shear_num[:-1],
                                         omega_i,
                                         omega_j)  # 递归

        # 计算Reyleigh阻尼的a,b
        a = 2 * omega_i * omega_j / (omega_i + omega_j) * child_damping_ratio[-1]
        b = 2 * child_damping_ratio[-1] / (omega_i + omega_j)
        damping_1 = damping_1 + b * stiffness_r
        damping_2 = b * child_stiffness_2
        damping_3 = damping_2.T
        damping_4 = a * child_mass_4 + b * child_stiffness_4

        # 将阻尼矩阵拼接好
        damping_temp1 = np.vstack((damping_1, damping_3))
        damping_temp2 = np.vstack((damping_2, damping_4))

        return np.hstack((damping_temp1, damping_temp2))
