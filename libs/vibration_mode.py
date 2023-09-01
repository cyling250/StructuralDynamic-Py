"""
程序描述：实用振动分析的计算程序，该文件记录了各种实用振动分析
的计算方法。
需要注意的是，当前的程序只保证了频率运算正确。
当然频率运算正确时，振型一般也正确，但是当前程序没有对振型进行归一化处理。
振型归一化有三种方法，特定坐标归一化、最大值归一化和关于质量矩阵归一化。
当前程序没有进行实现，在以后的修订过程中需要进行实现。
"""
import numpy as np
from scipy.linalg import *


def normalize(mat, psi):
    """
    将psi按照mat进行归一化
    """
    mat_num = (psi.T @ mat @ psi)
    return psi / np.sqrt(mat_num)


def rayleigh_psi(mass, stiffness, psi=None):
    """
    rayleigh振动分析方法
    Parameters
    ----------
    mass 质量矩阵
    stiffness 刚度矩阵
    psi 假设振型

    Returns 结构振动频率
    -------

    """
    if psi is None:
        length = len(mass)  # 获取矩阵的大小
        psi = np.array([1 for i in range(length)]).T
    numerator = psi.T @ stiffness @ psi
    denominator = psi.T @ mass @ psi
    omega = np.sqrt(numerator / denominator)  # 计算频率
    return omega


def rayleigh_ritz(mass, stiffness, psi=1, load=None):
    """
    rayleigh_ritz法计算程序。在外界没有给出假设振型的情况下，程序会自行计算一组合理的
    假设振型，计算方法属于ritz直接法，为荷载相关法。当外界连最基本的荷载都没有给出时，会
    自行在每个自由度生成大小为1的荷载，构成荷载列阵。
    Parameters
    ----------
    mass 质量矩阵
    stiffness 刚度矩阵
    psi 接受两种参数，当phi输入为int类型时，程序自动根据输入的
    值生成一组假设振型。当phi输入为ndarray类型时，程序直接读取假设振型。
    load 用于生成ritz向量的初始荷载

    Returns 一组频率，对应phi
    -------

    """
    if type(psi) == int:
        psi, alpha, beta = load_depended_ritz_vector(mass, stiffness, psi, load)  # 生成假设振型
    # 进行刚度矩阵和质量矩阵缩减
    stiffness_reduce = psi.T @ stiffness @ psi
    mass_reduce = psi.T @ mass @ psi
    # 解矩阵特征值问题
    [omega, gene_crd] = eig(stiffness_reduce, mass_reduce)
    omega = np.sqrt(omega)
    gene_crd = np.array(gene_crd)
    phi = np.zeros((len(mass), psi.shape[1]))
    # 求固有振型
    phi[:, 0] = psi @ gene_crd[:, 0]
    for i in range(1, psi.shape[1]):
        phi_temp = psi @ gene_crd[:, i]
        phi[:, i] = phi_temp
    return omega.real, phi


def load_depended_ritz_vector(mass, stiffness, count=1, load=None):
    """
    荷载相关的ritz向量，本函数将来可以用在振型叠加法中，用ritz向量
    代替结构自振向量。研究表明，振型叠加法计算中忽略了高阶振型的影响，
    会造成高阶振型作用较大时，结构计算收敛慢的现象，荷载相关ritz向量法
    可以很好地解决该问题。
    荷载相关ritz向量法中必须给出荷载列阵。
    Parameters
    ----------
    mass 质量矩阵
    stiffness 刚度矩阵
    count 需要求的ritz向量个数
    load 荷载矩阵

    Returns ritz向量,alpha数组,beta数组
    -------

    """
    if load is None:
        load = np.array([1 for i in range(len(mass))]).T
    stiffness_inv = inv(stiffness)  # 提前求逆，节省计算量
    # 第一阶Ritz向量计算，其定义为静力直接作用在结构上的位移
    dpm = stiffness_inv @ load
    ritz = np.zeros((len(mass), count + 1))
    ritz[:, 0] = normalize(mass, dpm)  # 按照质量矩阵归一化
    alpha = np.zeros(count)
    beta = np.zeros(count)
    # 第n+1阶Ritz向量计算
    for i in range(count):
        alpha_temp = np.zeros(i + 1)
        dpm = stiffness_inv @ mass @ ritz[:, i]  # 计算y[i+1]=K^-1*M*ritz[i]
        for j in range(i + 1):
            alpha_temp[j] = dpm.T @ mass @ ritz[:, j]  # 计算a[i+1,j]=y[i+1]*M*ritz[j]
            dpm -= alpha_temp[j] * ritz[:, j]  # 计算alpha向量的和
        ritz_temp1 = normalize(mass, dpm)  # 关于质量矩阵正交归一化
        ritz[:, i + 1] = ritz_temp1.T  # 整合ritz向量
        alpha[i] = alpha_temp[-1]
        beta[i] = np.sqrt(dpm.T @ mass @ dpm)
    return ritz[:, :-1], alpha, beta


def mat_iterate_base(mass, stiffness, precision=10e-3):
    """
    求解基准模态下的矩阵迭代法
    Parameters
    ----------
    mass 质量矩阵
    stiffness 刚度矩阵
    precision 精度

    Returns 频率，振型
    -------

    """
    # 计算起步条件
    # 计算动力矩阵
    stiffness_inv = inv(stiffness)
    dyn_mat = stiffness_inv @ mass
    # 计算其他条件
    freedom = mass.shape[0]  # 自由度
    phi = np.array([np.random.random() * 2 - 1 for i in range(freedom)]).T  # 随机生成初始振型
    omega_temp = (phi.T @ stiffness @ phi) ** 0.5
    step = 0  # 记录迭代步数
    # 开始迭代
    while True:
        phi_temp = dyn_mat @ phi  # 迭代主要部分
        phi = phi_temp / (phi_temp.T @ mass @ phi_temp) ** 0.5
        omega = (phi.T @ stiffness @ phi) ** 0.5
        if abs(omega - omega_temp) < precision * omega or step > 100:
            break
        omega_temp = omega
        step += 1
    print("迭代%d步收敛。" % step)
    return omega, phi


def mat_iterate_high(mass, stiffness, omega_pre, phi_pre, precision=10e-3):
    """
    求高阶模态的矩阵迭代法
    Parameters
    ----------
    mass 质量矩阵
    stiffness 刚度矩阵
    omega_pre 前n-1阶频率数组
    phi_pre 前n-1阶振型矩阵
    precision 计算精度

    Returns 频率，振型
    -------

    """
    # 计算动力矩阵
    stiffness_inv = inv(stiffness)
    dyn_mat = stiffness_inv @ mass
    # 计算其他条件
    freedom = mass.shape[0]  # 自由度
    phi = np.array([np.random.random() * 2 - 1 for i in range(freedom)]).T  # 随机生成初始振型
    omega_temp = (phi.T @ stiffness @ phi) ** 0.5
    step = 0  # 记录迭代步数

    # 开始迭代
    while True:
        phi_temp = dyn_mat @ phi  # 迭代主要部分
        phi = phi_temp / (phi_temp.T @ mass @ phi_temp) ** 0.5
        # 对高阶振型中的低阶振型进行消除
        alpha_vector = np.array([0.0 for i in range(mass.shape[1])]).T  # 初始化alpha向量
        for i in range(len(omega_pre)):
            alpha = phi.T @ mass @ phi_pre[i]
            alpha_vector += alpha * phi_pre[i]
        phi -= alpha_vector
        omega = (phi.T @ stiffness @ phi) ** 0.5
        if abs(omega - omega_temp) < precision * omega or step > 100:
            break  # 计算达到精度或迭代步数超出设定即跳出迭代
        omega_temp = omega
        step += 1
    print("迭代%d步收敛。" % step)
    return omega, phi


def mat_iterate_highest(mass, stiffness, precision=10e-3):
    """
    求最高阶频率的矩阵迭代法
    Parameters
    ----------
    mass 质量矩阵
    stiffness 刚度矩阵
    precision 计算精度

    Returns 频率，振型
    -------

    """
    # 计算起步条件
    # 计算动力矩阵
    mass_inv = inv(mass)
    dyn_mat = mass_inv @ stiffness
    # 计算其他条件
    freedom = mass.shape[0]  # 自由度
    phi = np.array([np.random.random() * 2 - 1 for i in range(freedom)]).T  # 随机生成初始振型
    omega_temp = (phi.T @ stiffness @ phi) ** 0.5
    step = 0  # 记录迭代步数
    # 开始迭代
    while True:
        phi_temp = dyn_mat @ phi  # 迭代主要部分
        phi = phi_temp / (phi_temp.T @ mass @ phi_temp) ** 0.5
        omega = (phi.T @ stiffness @ phi) ** 0.5
        if abs(omega - omega_temp) < precision * omega or step > 100:
            break
        omega_temp = omega
        step += 1
    print("迭代%d步收敛。" % step)
    return omega, phi


def subspace_iteration(mass, stiffness, psi=1, precision=10e-3):
    """
    子空间迭代法计算程序
    Parameters
    ----------
    mass 质量矩阵
    stiffness 刚度矩阵
    psi 假设振型，如果是int类型，则自动生成该数值的正交向量，否则psi应该直接就是正交向量
    precision 计算精度

    Returns 频率，振型
    -------

    """
    # 生成假设振型
    if type(psi) == int:
        psi, alpha, beta = load_depended_ritz_vector(mass, stiffness, psi)
    # 迭代前准备
    stiffness_inv = inv(stiffness)
    dyn_mat = stiffness_inv @ mass  # 计算动力矩阵
    omega_temp = np.ones(len(mass))  # 初始化精度控制矩阵
    step = 0  # 迭代步数
    # 开始迭代
    while True:
        # 矩阵迭代法迭代
        psi = dyn_mat @ psi
        # Rayleigh-ritz法迭代
        k_temp = psi.T @ stiffness @ psi  # 缩减质量矩阵
        m_temp = psi.T @ mass @ psi  # 缩减刚度矩阵
        [omega, comb_coeff] = eig(k_temp, m_temp)  # 计算特征值和特征向量
        omega = np.sqrt(omega)
        phi = psi @ comb_coeff
        if norm(omega - omega_temp) < norm(precision * omega) or step > 100:
            break
        omega_temp = omega
        step += 1
    return omega.real, phi


def lanczos(mass, stiffness, count=1, load=None):
    """
    兰索斯方法
    Parameters
    ----------
    mass 质量矩阵
    stiffness 刚度矩阵
    count 需要计算的振型数量
    load 荷载

    Returns 频率，振型
    -------

    """
    psi, alpha, beta = load_depended_ritz_vector(mass, stiffness, count, load)
    T = np.diag(alpha) + np.diag(beta[1:], k=1) + np.diag(beta[1:], k=-1)
    eigenvalues, eigenvectors = eigh(T)
    return eigenvalues, psi[:, :count] @ eigenvectors[:, :count]


def dunkerley(mass, stiffness):
    """
    邓克莱方法给出了一种估算体系基频近似值的方法，给出的是结构体系基本频率的下限，
    当其他各阶频率远远高于基频时，利用此法估算基频较为方便。
    (计算出的值相当不准确)
    Parameters
    ----------
    mass 质量矩阵
    stiffness 刚度矩阵

    Returns 频率
    -------

    """
    # 计算动力矩阵
    stiffness_inv = inv(stiffness)
    dyn_mat = stiffness_inv @ mass
    return np.sqrt(1 / dyn_mat.trace())


def jacobi(mass, stiffness, eps=1e-10, max_iter=1000):
    """
    Jacobi迭代法求解实对称矩阵的特征值和特征向量
    Parameters
    ----------
    mass 质量矩阵
    stiffness 刚度矩阵
    eps 精度
    max_iter 最大迭代步数

    Returns 频率,振型
    -------

    """
    # 初始特征向量为单位矩阵
    length = len(mass)
    # 计算动力矩阵
    stiffness_inv = inv(stiffness)
    dyn_mat = stiffness_inv @ mass  # 计算动力矩阵
    v_revlove = np.eye(length)
    dyn_np = np.array(dyn_mat)
    # 迭代计算
    count = 0
    while count < max_iter:
        # 计算非对角元素的最大值和位置
        # 此处这样写是因为matrix格式的数据使用flatten拉直之后
        # 最外层仍为二维数组，所以需要改进
        max_idx = np.argmax(np.abs(np.triu(dyn_np, 1)))  # 找到最大的序号
        if dyn_np.flatten()[max_idx] ** 2 < eps:
            break
        i, j = divmod(max_idx, length)

        # 计算旋转角度
        if dyn_np[i, i] == dyn_np[j, j]:
            theta = np.pi / 4
        else:
            theta = 0.5 * np.arctan(2 * dyn_np[i, j] / (dyn_np[i, i] - dyn_np[j, j]))

        # 构造旋转矩阵
        r_revlove = np.eye(length)
        r_revlove[i, i] = np.cos(theta)
        r_revlove[j, j] = np.cos(theta)
        r_revlove[i, j] = -np.sin(theta)
        r_revlove[j, i] = np.sin(theta)

        # 更新A和V
        dyn_np = r_revlove.T @ dyn_np @ r_revlove
        v_revlove = v_revlove @ r_revlove

        count += 1

    # 提取特征值和特征向量
    eig_values = np.diag(dyn_np)
    eig_vectors = v_revlove.T
    omega = 1 / np.sqrt(eig_values)
    # 返回结果
    return omega, eig_vectors
