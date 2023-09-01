"""
振型叠加法计算程序
"""
import numpy as np
from scipy.linalg import *


def fourier(mass, stiffness, load, delta_time, damping_ratio=0.05):
    """
    傅里叶变换求解地震响应程序
    Parameters
    ----------
    mass 质量
    stiffness 刚度
    damping_ratio 阻尼比
    load 荷载时程
    delta_time 时间步长

    Returns 时程响应
    -------

    """
    load_fft = np.fft.fft(load)  # 对荷载进行傅里叶变换
    omega_n = np.sqrt(stiffness / mass)
    # fs = 1 / delta_time  # 采样频率
    length = len(load_fft)
    result = np.zeros(length, dtype=complex)  # 创建复数矩阵用来存储计算得到的位移谱
    for i in range(length):
        omega = 2 * np.pi * i / (delta_time * length)  # 当前的周期
        complex_reaction = complex(1 - (omega / omega_n) ** 2,
                                   2 * damping_ratio * (omega / omega_n)) * stiffness  # 复频响应函数
        result[i] = load_fft[i] / complex_reaction
    result = np.fft.ifft(result)
    return 2 * result.real  # 只需要取实部并乘2


def modal_superposition(mass, stiffness, load, delta_time, damping_ratio=0.05):
    """
    振型分解法计算程序
    Parameters
    ----------
    mass 质量矩阵
    stiffness 刚度矩阵
    load 时序荷载列阵

    Returns
    -------

    """
    # 数据准备
    freedom = len(mass)
    length = len(load)
    [omega, phi] = eig(stiffness, mass)
    # 按照升序排列
    sort_num = omega.argsort()
    phi = phi.T[sort_num].T
    # 坐标变换
    split_mass = np.diag(np.transpose(phi) @ mass @ phi)
    split_stiffness = np.diag(np.transpose(phi) @ stiffness @ phi)
    load = np.array([np.transpose(phi) @ load[i] for i in range(length)])
    # 计算变换后的坐标
    dpm = np.zeros((length, freedom))
    for i in range(freedom):
        dpm[:, i] = fourier(split_mass[i], split_stiffness[i], load[:, i], delta_time, damping_ratio)
    # 振型坐标转换为位移坐标
    for i in range(length):
        dpm[i] = phi @ dpm[i]
    return dpm


def complex_modal_superposition(mass, stiffness, load, delta_time, damping):
    """
    复模态分析法
    Parameters
    ----------
    mass 质量矩阵
    stiffness 刚度矩阵
    load 荷载列阵
    delta_time 时间间隔

    Returns
    -------

    """
    # 数据准备
    freedom = len(mass)
    length = len(load)
    zero_mat = np.zeros((freedom, freedom))

    unit = np.zeros(freedom)
    unit[-1] = 1
    mass_temp = np.hstack((damping, mass))
    mass_complex = np.hstack((mass, zero_mat))
    mass_complex = np.vstack((mass_temp, mass_complex))  # 复模态质量矩阵
    stiffness_temp = np.hstack((stiffness, zero_mat))
    stiffness_complex = np.hstack((zero_mat, -mass))
    stiffness_complex = np.vstack((stiffness_temp, stiffness_complex))  # 复模态刚度矩阵

    unit_complex = np.hstack((np.zeros(freedom), unit))  # 2N维单位向量
    [lambda_, phi] = eig(stiffness_complex, mass_complex)
    # 这里其实已经计算出了lambda，可以直接使用，但是我们还是再按照书上的方法再计算一边
    real_lambda = lambda_[0::2]
    complex_lambda = lambda_[1::2]
    real_phi = phi[:, 0::2]
    complex_phi = phi[:, 1::2]
    lambda_ = np.hstack((real_lambda, complex_lambda))
    phi = np.hstack((real_phi, complex_phi))
    # 单自由度体系系数lambda,nita
    a = np.zeros(2 * freedom, dtype=complex)
    b = np.zeros(2 * freedom, dtype=complex)
    nita = np.zeros(2 * freedom, dtype=complex)
    for i in range(len(mass_complex)):
        a[i] = phi[:, i].T @ mass_complex @ phi[:, i]
        b[i] = phi[:, i].T @ stiffness_complex @ phi[:, i]
        nita[i] = phi[:, i].T @ mass_complex @ unit_complex / a[i]
    # 开始数值积分
    dpm_temp = np.zeros((length, freedom), dtype=complex)
    dpm = np.zeros((length, freedom), dtype=complex)
    for i in range(2 * freedom):
        quad = 0
        for j in range(length):
            quad += np.e ** (lambda_[i] * j * delta_time) * load[j] * delta_time
            dpm_temp[j] = phi[0:freedom, i] * nita[i] * np.e ** (-lambda_[i] * j * delta_time) * quad
        dpm += dpm_temp
    return np.real(dpm)


def pse_spectrum(quake_wave, delta_time, func):
    """
    反应谱计算函数,这里计算的加速度反应谱是伪加速度反应谱
    Parameters
    ----------
    quake_wave 需要计算反应谱的地震波
    delta_time 时间间隔
    func 计算反应谱的方法

    Returns 位移谱，速度谱，加速度谱
    -------

    """
    ts = np.arange(0.01, 1, 0.01)
    ts = np.hstack((ts, np.arange(1.05, 2, 0.05)))
    ts = np.hstack((ts, np.arange(2.1, 3, 0.1)))
    ts = np.hstack((ts, np.arange(3.2, 4, 0.2)))
    ts = np.hstack((ts, np.arange(4.5, 6, 0.5)))
    length = len(ts)
    spectrum_dpm, spectrum_vel, spectrum_acc = np.zeros(length), np.zeros(length), np.zeros(length)
    for i in range(len(ts)):
        omega = 2 * np.pi / ts[i]
        k = omega ** 2
        print("\r正在计算Ts:%.3fs" % ts[i], end="")
        # func是传进来的不同的计算方法，其传参格式均被规整化了，可以直接调用
        dpm = func(mass=1, stiffness=k, load=quake_wave, delta_time=delta_time, damping_ratio=0.05)
        # 对位移谱进行傅里叶逆变换，得到时域位移，并取响应最大值
        if type(dpm) == tuple:  # 如果是用逐步积分法算出来的
            dpm = dpm[0]
        spectrum_dpm[i] = np.max(np.abs(dpm))
        spectrum_vel[i] = spectrum_dpm[i] * omega
        spectrum_acc[i] = spectrum_vel[i] * omega
    print()
    return spectrum_dpm, spectrum_vel, spectrum_acc


def abs_spectrum(quake_wave, delta_time, func):
    """
    反应谱计算函数,这里计算的加速度反应谱是绝对反应谱
    并且由于绝对反应谱计算时需要计算加速度响应，所以fourier变换
    和duhamel积分适用性不再很好，这里直接采用分段解析法计算
    Parameters
    ----------
    quake_wave 需要计算反应谱的地震波
    delta_time 时间间隔

    Returns 位移谱，速度谱，加速度谱
    -------

    """
    ts = np.arange(0.01, 1, 0.01)
    ts = np.hstack((ts, np.arange(1.05, 2, 0.05)))
    ts = np.hstack((ts, np.arange(2.1, 3, 0.1)))
    ts = np.hstack((ts, np.arange(3.2, 4, 0.2)))
    ts = np.hstack((ts, np.arange(4.5, 6, 0.5)))
    length = len(ts)
    spectrum_dpm, spectrum_vel, spectrum_acc = np.zeros(length), np.zeros(length), np.zeros(length)
    for i in range(len(ts)):
        omega = 2 * np.pi / ts[i]
        k = omega ** 2
        print("\r正在计算Ts:%.3fs" % ts[i], end="")
        # func是传进来的不同的计算方法，其传参格式均被规整化了，可以直接调用
        dpm, vel, acc = func(mass=1, stiffness=k, load=quake_wave, delta_time=delta_time,
                             damping_ratio=0.05)
        # 统计并记录
        spectrum_dpm[i] = np.max(np.abs(dpm))
        spectrum_vel[i] = np.max(np.abs(vel))
        spectrum_acc[i] = np.max(np.abs(acc))
    print()
    return spectrum_dpm, spectrum_vel, spectrum_acc


def design_spectrum(ts, damping_ratio, t_g=0.35, alpha_max=0.08):
    """
    依据抗震规范定义的设计反应谱函数
    Parameters
    ----------
    ts 结构周期
    damping_ratio 结构阻尼比
    t_g 特征周期
    alpha_max 地震影响系数最大值

    Returns 地震影响系数
    -------

    """
    gamma = 0.9 + (0.05 - damping_ratio) / (0.3 + 6 * damping_ratio)
    nita1 = 0.02 + (0.05 - damping_ratio) / (4 + 32 * damping_ratio)
    nita2 = 1 + (0.05 - damping_ratio) / (0.08 + 1.6 * damping_ratio)

    if nita1 < 0:
        nita1 = 0
    if nita2 < 0.55:
        nita2 = 0.55

    if ts < 0.1:
        return 0.45 + alpha_max + ts * (nita2 - 0.45) * alpha_max / 0.1
    elif ts < t_g:
        return nita2 * alpha_max
    elif ts < 5 * t_g:
        return (t_g / ts) ** gamma * nita2 * alpha_max
    else:
        return (nita2 * 0.2 ** gamma - nita1 * (ts - 0.5 * t_g)) * alpha_max


def srss(arr, damping_ratio=None, omega=None):
    return np.sqrt(np.sum(np.square(arr), axis=1))


def cqc(arr, damping_ratio=None, omega=None):
    freedom = len(arr)
    # 计算偶联系数
    result = np.zeros(freedom)
    for i in range(freedom):
        for j in range(freedom):
            # 计算周期之比
            lambda_ = omega[j] / omega[i]  # 周期之比为频率的反比
            if lambda_ > 1:
                # 总是让该系数小于1
                lambda_ = 1 / lambda_
            # print(lambda_)
            temp1 = 8 * (damping_ratio[i] * damping_ratio[j]) ** 0.5 * (
                    damping_ratio[i] + lambda_ * damping_ratio[j]) * lambda_ ** 1.5
            temp2 = (1 - lambda_ ** 2) ** 2 + 4 * damping_ratio[i] * damping_ratio[j] * (
                    1 + lambda_ ** 2) * lambda_ + 4 * (
                            damping_ratio[i] ** 2 + damping_ratio[j] ** 2) * lambda_ ** 2
            rho = temp1 / temp2

            result += rho * arr[:, i] * arr[:, j]

    return np.sqrt(result)


def modal_mass(mass, phi):
    """
    振型参与质量系数
    Parameters
    ----------v
    mass 刚度矩阵
    phi 该结构的振型

    Returns 振型参与质量系数数组
    -------

    """
    freedom = len(mass)
    gamma_n = np.zeros(freedom)  # 振型参与系数
    each_shear_m = np.diagonal(mass)  # 各自由度的质量
    mass_e = np.zeros((freedom, freedom))  # 各振型在各层的参与质量
    result = np.zeros(freedom)  # 振型参与质量系数
    for i in range(freedom):
        # 计算振型参与质量系数
        mass_n = phi[:, i].T @ mass @ phi[:, i]
        gamma_n[i] = (phi[:, i].T @ mass @ np.ones(freedom)) / mass_n
    for i in range(freedom):
        for j in range(freedom):
            mass_e[i, j] = gamma_n[j] * phi[i, j] * each_shear_m[i]
    mass_e = np.sum(mass_e, axis=1)
    mass_sum = np.sum(mass_e)
    for i in range(freedom):
        result[i] = sum(mass_e[:i + 1]) / mass_sum

    return result


def modal_response_spectrum(mass, stiffness, damping_ratio=None, func=srss):
    """
    振型分解反应谱法计算地震响应，本程序计算模态时直接采用了Python中给出的eig函数，
    使用vibration_mode.py文件中的任一计算方法也是可行的。
    此外，因为本文例子的振型数少，所以本函数直接采用了所有振型参与计算，而非严格按照规范讨论前90%
    但给出了计算各振型的振型参与质量系数的函数。
    Parameters
    ----------
    mass 质量矩阵
    stiffness 刚度矩阵
    damping_ratio 阻尼比
    func 计算方法:SRSS或者CQC

    Returns 楼层剪力
    -------

    """
    [omega, phi] = eig(stiffness, mass)  # 计算频率，振型
    sort_num = omega.argsort()  # 计算真实频率，开平方
    # 按照升序排列
    omega = np.sqrt(omega[sort_num]).real
    phi = phi.T[sort_num].T
    freedom = len(mass)  # 自由度
    gamma_n = np.zeros(freedom)  # 振型参与系数
    equ_quake_force = np.zeros((freedom, freedom))  # 各振型的等效地震力
    force = np.zeros((freedom, freedom))  # 真实楼层剪力
    alpha = np.zeros(freedom)  # 各振型的地震影响系数
    modal_damping = np.zeros(freedom)

    for i in range(freedom):
        # 计算振型参与系数
        mass_n = phi[:, i].T @ mass @ phi[:, i]
        gamma_n[i] = (phi[:, i].T @ mass @ np.ones(freedom)) / mass_n
        # 计算等效地震力
        ts = 2 * np.pi / omega[i]
        temp_alpha = 2 * damping_ratio / (omega[0] + omega[1]) * np.array([omega[0] * omega[1], 1])
        modal_damping[i] = temp_alpha[0] / (2 * omega[i]) + temp_alpha[1] * omega[i] / 2
        alpha[i] = design_spectrum(ts, modal_damping[i])
        equ_quake_force[:, i] = alpha[i] * 9.8 * gamma_n[i] * mass @ phi[:, i]
        for j in range(freedom):
            force[j, i] = np.sum(equ_quake_force[j:freedom, i])

    return func(force, modal_damping, omega)
