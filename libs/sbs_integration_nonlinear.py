"""本程序是结构非线性计算的库"""
import numpy as np
from scipy.linalg import *
from libs.damping import rayleigh
from scipy import integrate
from libs.layer_shear import get_k


class Bilinear2:
    def __init__(self, stiffness, stiff_nlr, yeild_force, dpm, vel, force):
        """
        弹性刚度，屈服刚度，屈服应力，初始位移，初始速度，初始应力
        """
        self.stiffness = stiffness  # 弹性刚度
        self.stiff_nlr = stiff_nlr  # 屈服刚度
        self.flag = 0  # 当前标志位;0:加载;1:屈服;2:卸载
        self.yeild_force = yeild_force
        self.tensile_max = [yeild_force / self.stiffness, yeild_force]  # 初始化最大拉伸应力点
        self.compress_max = [-yeild_force / self.stiffness, -yeild_force]  # 初始化最大压缩应力点
        self.dpm_init = 0  # 卸载阶段直线的x轴截距
        self.dpm = dpm  # 上一时刻的位移
        self.vel = vel  # 上一时刻的速度
        self.force = force  # 上一时刻的应力

    def __call__(self, dpm, vel):
        """调用该对象时，需要输出当前的恢复力与刚度"""
        self.flag = self.state_change(dpm, vel)  # 判断当前状态并转换
        if self.flag == 0:
            # 如果是加载阶段
            if dpm >= self.dpm_init:
                # 正向加载
                # stiff = (self.tensile_max[1] - self.force) / (self.tensile_max[0] - self.dpm)
                stiff = self.stiffness  # 要改成不退化的模型，请启用这一行，注释掉上一行
                force = stiff * (dpm - self.dpm) + self.force
                self.dpm_init = -force / self.stiffness + dpm
            else:
                # 反向加载
                # stiff = (self.compress_max[1] - self.force) / (self.compress_max[0] - self.dpm)
                stiff = self.stiffness  # 要改成不退化的模型，请启用这一行，注释掉上一行
                force = stiff * (dpm - self.dpm) + self.force
                self.dpm_init = -force / self.stiffness + dpm
        elif self.flag == 1:
            # 如果是屈服阶段
            if self.force >= 0:
                # 正向屈服
                force = self.stiff_nlr * (dpm - self.tensile_max[0]) + self.tensile_max[1]
                self.dpm_init = -force / self.stiffness + dpm
                stiff = self.stiff_nlr
            else:
                # 反向屈服
                force = self.stiff_nlr * (dpm - self.compress_max[0]) + self.compress_max[1]
                self.dpm_init = -force / self.stiffness + dpm
                stiff = self.stiff_nlr
        else:
            # 如果是卸载阶段
            stiff = self.stiffness
            force = stiff * (dpm - self.dpm_init)

        self.dpm = dpm
        self.vel = vel
        self.force = force
        # 要改成不退化的模型，请启用下面几行
        if self.force > self.yeild_force:
            self.force = self.yeild_force
        elif self.force < -self.yeild_force:
            self.force = -self.yeild_force
        return stiff, self.force

    def state_change(self, dpm, vel):
        """
        状态转换函数，用来判断加载-屈服-卸载的状态转换
        每一个小分支都要return
        """
        if self.flag == 0:
            # 如果是加载阶段
            if vel * self.vel < 0:
                # 如果速度方向改变
                return 2  # 卸载
            elif dpm > self.tensile_max[0] or dpm < self.compress_max[0]:
                return 1  # 屈服
            return 0  # 继续加载
        elif self.flag == 1:
            if vel * self.vel < 0:
                # 更新最大应力点
                if self.force >= 0:
                    # 拉伸状态
                    self.tensile_max = [self.dpm, self.force]
                else:
                    # 压缩状态
                    self.compress_max = [self.dpm, self.force]
                return 2  # 卸载
            return 1  # 继续屈服
        else:
            if vel * self.vel < 0:
                return 0  # 加载
            elif (dpm - self.dpm_init) * vel > 0:
                # 一个强制的拐点处理，可以提高计算精度
                self.dpm = self.dpm_init
                self.force = 0
                return 0  # 加载
            return 2  # 继续卸载


class Bilinear3(Bilinear2):
    def __init__(self, stiffness, stiff_nlr1, stiff_nlr2, crack_force, yeild_force, dpm, vel, force):
        # 请注意，这里使用继承的开裂阶段双折线模型一定是坡顶式的
        super().__init__(stiffness, stiff_nlr1, crack_force, dpm, vel, force)  # 继承父类属性
        self.stiff_nlr2 = stiff_nlr2
        self.yeild_force2 = yeild_force

    def is_yeild(self):
        # 判断是否完全开裂
        if self.flag == 1 and abs(self.force) > abs(self.yeild_force2):
            # 如果处于开裂阶段并且恢复力大于屈服强度，则走完开裂阶段，更新弹性刚度于屈服刚度
            stiff = self.yeild_force2 / (
                    self.yeild_force / self.stiffness + (self.yeild_force2 - self.yeild_force) / self.stiff_nlr)
            self.yeild_force = self.yeild_force2  # 更新屈服强度
            self.stiffness = stiff  # 更新弹性刚度
            self.stiff_nlr = self.stiff_nlr2  # 更新屈服刚度
            if self.force > 0:
                self.force = self.yeild_force
                self.tensile_max[1] = self.force
                self.compress_max[1] = -self.force
            else:
                self.force = -self.yeild_force
                self.tensile_max[1] = -self.force
                self.compress_max[1] = self.force

    def __call__(self, dpm, vel):
        stiff, self.force = Bilinear2.__call__(self, dpm, vel)
        self.is_yeild()
        return stiff, self.force


class Boucwen:
    def __init__(self, stiffness, alpha, A, beta, gamma, n, dpm, vel, z, force):
        self.stiffness = stiffness  # 材料弹性刚度
        self.alpha = alpha  # 屈服前后刚度比
        # 材料参数
        self.A = A
        self.beta = beta
        self.gamma = gamma
        self.n = n
        self.dpm = dpm  # 上一步位移
        self.vel = vel  # 上一步速度
        self.vel_temp = vel  # 这一步速度，迭代变量
        self.z = z  # 上一步滞变位移
        self.force = force  # 上一步恢复力

    def __call__(self, dpm, vel, time):
        self.vel_temp = vel
        self.z = integrate.odeint(self.function, self.z, time)[-1][0]
        temp_force = self.alpha * self.stiffness * dpm + (1 - self.alpha) * self.stiffness * self.z
        temp_stiff = self.stiffness
        if temp_force - self.force != 0:
            temp_stiff = (temp_force - self.force) / (dpm - self.dpm)
        self.force = temp_force
        self.dpm = dpm
        self.vel = vel
        return temp_stiff, self.force

    def function(self, z, t):
        z_dot = self.A * self.vel_temp - self.beta * abs(self.vel_temp) * abs(z) ** (
                self.n - 1) * z - self.gamma * self.vel_temp * abs(z) ** self.n
        return z_dot


def bili_mdof(dpm, vel, states):
    """
    非线性刚度计算
    Parameters
    ----------
    dpm 位移矩阵
    vel 速度矩阵
    states 材料状态

    Returns 刚度矩阵，恢复力矩阵
    -------

    """
    # 初始化相对矩阵
    freedom = len(dpm)
    # 相对位移与绝对位移的基变换矩阵,dpm@trans_mat可以将绝对位移变成相对位移
    trans_mat = np.diag([1 for i in range(freedom)]) + np.diag([-1 for i in range(freedom - 1)], 1)
    # 恢复力乘的矩阵需要特殊注意，需要乘这个plus矩阵而不是上面的那个矩阵
    trans_mat_plus = np.diag([1 for i in range(freedom)]) + np.diag([-1 for i in range(freedom - 1)], -1)
    dpm = dpm @ trans_mat
    vel = vel @ trans_mat
    stiff_new = np.zeros(freedom)
    force = np.zeros(freedom)
    for i in range(0, freedom):
        stiff_new[i], force[i] = states[i](dpm[i], vel[i])
    # 重新组装刚度矩阵
    stiff_new = get_k(stiff_new)
    # 基变换，将以相对位移为基的恢复力转换为绝对位移
    force = force @ trans_mat_plus
    return stiff_new, force


def bouc_mdof(dpm, vel, states, time):
    # 初始化相对矩阵
    freedom = len(dpm)
    # 相对位移与绝对位移的基变换矩阵,dpm@trans_mat可以将绝对位移变成相对位移
    trans_mat = np.diag([1 for i in range(freedom)]) + np.diag([-1 for i in range(freedom - 1)], 1)
    # 恢复力乘的矩阵需要特殊注意，需要乘这个plus矩阵而不是上面的那个矩阵
    trans_mat_plus = np.diag([1 for i in range(freedom)]) + np.diag([-1 for i in range(freedom - 1)], -1)
    dpm = dpm @ trans_mat
    vel = vel @ trans_mat
    force = np.zeros(freedom)
    temp_stiff = np.zeros(freedom)
    for i in range(0, freedom):
        temp_stiff[i], force[i] = states[i](dpm[i], vel[i], time)
    temp_stiff = temp_stiff @ inv(trans_mat)
    # 基变换，将以相对位移为基的恢复力转换为绝对位移
    force = force @ trans_mat_plus

    return temp_stiff, force


def center_bili(mass, stiffness, load,
                delta_time, damping_ratio,
                dpm_0, vel_0, model_parms,
                result_length=0,
                nlr_model=Bilinear2):
    """
    中心差分法多自由度非线性计算函数。
    Parameters
    ----------
    load 时序荷载列阵;二维矩阵
    delta_time 采样间隔;浮点数
    mass 质量矩阵;二维矩阵
    stiffness 刚度矩阵;二维矩阵
    damping_ratio 阻尼比;浮点数
    dpm_0 初始位移;一维数组
    vel_0 初始速度;一维数组
    result_length 结果长度;整数

    Returns 位移，速度，加速度;数组列表，数组列表，数组列表
    -------

    """
    # 固有属性计算
    freedom = len(dpm_0)  # 计算自由度
    damping = rayleigh(mass, stiffness, damping_ratio)  # rayleigh阻尼

    # 起步条件计算
    acc_0 = inv(mass) @ (load[0] - damping @ vel_0 - stiffness @ dpm_0)
    dpm_minus1 = dpm_0 - delta_time * vel_0 + delta_time ** 2 * acc_0 / 2

    # 前置参数计算
    equ_k = mass / delta_time ** 2 + damping / (2 * delta_time)  # 等效刚度矩阵
    a = (2 * mass) / delta_time ** 2
    b = mass / delta_time ** 2 - damping / (2 * delta_time)

    # 荷载长度处理
    if result_length == 0:
        result_length = int(1.2 * len(load))
    load = load + [np.array([0 for i in range(freedom)]) for i in range(result_length - len(load))]

    # 初始化位移、速度、加速度时程矩阵
    dpm = np.zeros((result_length, freedom))
    vel = np.zeros((result_length, freedom))
    acc = np.zeros((result_length, freedom))
    dpm[0] = dpm_0
    vel[0] = vel_0
    # 初始化非线性参数
    models = [nlr_model(i[0], i[1], i[2], i[3], i[4], i[5]) for i in model_parms]
    force = np.zeros((result_length, freedom))

    # 迭代开始
    # 起步
    equ_p = load[0] - (stiffness - a) @ dpm_0 - b @ dpm_minus1
    dpm[1] = inv(equ_k) @ equ_p

    for i in range(1, result_length - 1):
        # print("Step: %d" % i)
        # 计算瞬时恢复力
        delta_stiffness, force[i] = bili_mdof(dpm[i], dpm[i] - dpm[i - 1], models)
        # 计算等效荷载
        equ_p = load[i] - force[i] + a @ dpm[i] - b @ dpm[i - 1]
        # 计算下一步位移
        dpm[i + 1] = inv(equ_k) @ equ_p
        # 计算下一步速度
        vel[i] = (dpm[i + 1] - dpm[i - 1]) / (2 * delta_time)
        # 计算下一步加速度
        acc[i] = (dpm[i + 1] - 2 * dpm[i] + dpm[i - 1]) / delta_time ** 2

    return dpm, vel, acc, force


def center_bouc(mass, stiffness, load,
                delta_time, damping_ratio,
                dpm_0, vel_0, model_parms,
                result_length=0,
                nlr_model=Boucwen):
    """
    中心差分法多自由度非线性计算函数。
    Parameters
    ----------
    load 时序荷载列阵;二维矩阵
    delta_time 采样间隔;浮点数
    mass 质量矩阵;二维矩阵
    stiffness 刚度矩阵;二维矩阵
    damping_ratio 阻尼比;浮点数
    dpm_0 初始位移;一维数组
    vel_0 初始速度;一维数组
    result_length 结果长度;整数

    Returns 位移，速度，加速度;数组列表，数组列表，数组列表
    -------

    """
    # 固有属性计算
    freedom = len(dpm_0)  # 计算自由度
    damping = rayleigh(mass, stiffness, damping_ratio)  # rayleigh阻尼

    # 起步条件计算
    acc_0 = inv(mass) @ (load[0] - damping @ vel_0 - stiffness @ dpm_0)
    dpm_minus1 = dpm_0 - delta_time * vel_0 + delta_time ** 2 * acc_0 / 2

    # 前置参数计算
    equ_k = mass / delta_time ** 2 + damping / (2 * delta_time)  # 等效刚度矩阵
    a = (2 * mass) / delta_time ** 2
    b = mass / delta_time ** 2 - damping / (2 * delta_time)

    # 荷载长度处理
    if result_length == 0:
        result_length = int(1.2 * len(load))
    load = load + [np.array([0 for i in range(freedom)]) for i in range(result_length - len(load))]

    # 初始化位移、速度、加速度时程矩阵
    dpm = np.zeros((result_length, freedom))
    vel = np.zeros((result_length, freedom))
    acc = np.zeros((result_length, freedom))
    dpm[0] = dpm_0
    vel[0] = vel_0
    # 初始化非线性参数
    models = [nlr_model(i[0], i[1], i[2], i[3], i[4], i[5], i[6], i[7], i[8], i[9]) for i in model_parms]
    force = np.zeros((result_length, freedom))

    # 迭代开始
    # 起步
    equ_p = load[0] - (stiffness - a) @ dpm_0 - b @ dpm_minus1
    dpm[1] = inv(equ_k) @ equ_p

    for i in range(1, result_length - 1):
        # 计算瞬时恢复力
        # 在双折线模型中，速度其实并不参与计算而是作为一个评估指标，所以这里采用近似的方法
        temp_stiff, force[i] = bouc_mdof(dpm[i], vel[i - 1], models, np.array([i, i + 1]) * delta_time)
        # 计算等效荷载
        equ_p = load[i] - force[i] + a @ dpm[i] - b @ dpm[i - 1]
        # 计算下一步位移
        dpm[i + 1] = inv(equ_k) @ equ_p
        # 计算下一步速度
        vel[i] = (dpm[i + 1] - dpm[i - 1]) / (2 * delta_time)
        # 计算下一步加速度
        acc[i] = (dpm[i + 1] - 2 * dpm[i] + dpm[i - 1]) / delta_time ** 2

    return dpm, vel, acc, force


def newmark_bili(mass, stiffness, load,
                 delta_time, damping_ratio,
                 dpm_0, vel_0, acc_0,
                 model_parms,
                 beta=0.25, gamma=0.5,
                 result_length=0,
                 nlr_model=Bilinear2):
    """
    MDOF非线性Newmark-beta法计算函数
    Parameters
    ----------
    mass 质量矩阵
    stiffness 刚度矩阵
    load 荷载列阵
    delta_time 时间步长
    damping_ratio 阻尼比
    dpm_0 初始位移
    vel_0 初始速度
    acc_0 初始加速度
    model_parms 非线性模型参数
    beta beta参数
    gamma gamma参数
    result_length 计算长度
    nlr_model 采用的非线性模型

    Returns 位移,速度,加速度,恢复力
    -------

    """
    # 进行计算数据准备
    freedom = len(dpm_0)  # 计算自由度
    if type(damping_ratio) == float:
        damping = rayleigh(mass, stiffness, damping_ratio)  # 计算阻尼矩阵
    else:
        damping = damping_ratio

    # 荷载末端补零
    if result_length == 0:
        result_length = int(1.2 * len(load))  # 计算持时
    load = load + [np.array([0 for i in range(freedom)]) for i in range(result_length - len(load))]

    # 初始化位移、速度、加速度
    dpm = np.zeros((result_length, freedom))
    vel = np.zeros((result_length, freedom))
    acc = np.zeros((result_length, freedom))
    dpm[0] = dpm_0
    vel[0] = vel_0
    acc[0] = acc_0

    # 初始化非线性参数
    models = [nlr_model(i[0], i[1], i[2], i[3], i[4], i[5]) for i in model_parms]
    force = np.zeros((result_length, freedom))

    # 计算中间参数
    a_0 = 1 / (beta * delta_time ** 2)
    a_1 = gamma / (beta * delta_time)
    a_2 = 1 / (beta * delta_time)
    a_3 = 1 / (2 * beta)
    a_4 = gamma / beta
    a_5 = (a_4 - 2) * delta_time / 2
    temp_stiff = stiffness
    # 积分步迭代开始
    for i in range(result_length - 1):
        # 计算荷载增量
        delta_load = load[i + 1] - load[i]
        # 计算等效荷载
        equ_load = delta_load + mass @ (a_2 * vel[i] + a_3 * acc[i]) + damping @ (a_4 * vel[i] + a_5 * acc[i])
        # 计算等效刚度
        equ_stiffness = temp_stiff + a_0 * mass + a_1 * damping
        # newton-raphson迭代
        # 初始化迭代变量
        delta_dpm = dpm[i]  # 迭代初始位移
        delta_force = force[i]  # 迭代初始恢复力
        # print("Step :%d" % i)
        while True:
            # print("等效荷载:", equ_load)
            delta_temp_dpm = inv(equ_stiffness) @ equ_load  # 计算迭代步位移增量
            # print("位移增量:", delta_temp_dpm)
            temp_dpm = delta_dpm + delta_temp_dpm  # 求解迭代后位移
            # print("迭代后位移:", temp_dpm)
            # 计算当前恢复力与刚度
            temp_stiff, temp_force = bili_mdof(temp_dpm, vel[i], models)
            # 计算等效残余力
            # print("上一步恢复力:", delta_force)
            # print("当前恢复力:", temp_force)
            # print("恢复力增量:", temp_force - delta_force)
            # print("当前刚度矩阵:", temp_stiff)
            # print("当前刚度矩阵乘位移增量:", temp_stiff @ delta_temp_dpm)
            # print("当前刚度矩阵乘位移:", temp_stiff @ temp_dpm)
            # print("按照等式右边计算出的等效荷载:", (temp_stiff + a_0 * mass + a_1 * damping) @ delta_temp_dpm)
            equ_load -= temp_force - delta_force  # 这一步恢复力减去上一步恢复力
            equ_stiffness = temp_stiff + a_0 * mass + a_1 * damping
            # if is_dyn:
            # 如果已经减过，则不需要再减
            # print("减去动力项")
            equ_load -= (a_0 * mass + a_1 * damping) @ delta_temp_dpm
            # is_dyn = False

            # 更新位移与恢复力
            delta_dpm = temp_dpm
            delta_force = temp_force
            # print(equ_load)
            if sum(abs(equ_load)) < 1e-5:
                break
        # 计算真实位移
        dpm[i + 1] = temp_dpm
        # 计算真实加速度
        acc[i + 1] = a_0 * (dpm[i + 1] - dpm[i]) - a_2 * vel[i] - (a_3 - 1) * acc[i]
        # 计算速度
        vel[i + 1] = a_1 * (dpm[i + 1] - dpm[i]) + (1 - a_4) * vel[i] + (1 - a_4 / 2) * acc[i] * delta_time
        force[i + 1] = temp_force

    return dpm, vel, acc, force


def newmark_bouc(mass, stiffness, load,
                 delta_time, damping_ratio,
                 dpm_0, vel_0, acc_0,
                 model_parms,
                 beta=0.25, gamma=0.5,
                 result_length=0,
                 nlr_model=Boucwen):
    """
    Newmark-beta法不适用于Boucwen模型
    Parameters
    ----------
    mass 质量矩阵
    stiffness 刚度矩阵
    load 荷载列阵
    delta_time 时间步长
    damping_ratio 阻尼比
    dpm_0 初始位移
    vel_0 初始速度
    acc_0 初始加速度
    model_parms 非线性模型参数
    beta beta参数
    gamma gamma参数
    result_length 计算长度
    nlr_model 采用的非线性模型

    Returns 位移,速度,加速度,恢复力
    -------

    """
    # 进行计算数据准备
    freedom = len(dpm_0.T)  # 计算自由度
    damping = rayleigh(mass, stiffness, damping_ratio)  # 计算阻尼矩阵

    # 荷载末端补零
    if result_length == 0:
        result_length = int(1.2 * len(load))  # 计算持时
    load = load + [np.array([0 for i in range(freedom)]) for i in range(result_length - len(load))]

    # 初始化位移、速度、加速度
    dpm = np.zeros((result_length, freedom))
    vel = np.zeros((result_length, freedom))
    acc = np.zeros((result_length, freedom))
    dpm[0] = dpm_0
    vel[0] = vel_0
    acc[0] = acc_0

    # 初始化非线性参数
    models = [nlr_model(i[0], i[1], i[2], i[3], i[4], i[5], i[6], i[7], i[8], i[9]) for i in model_parms]
    force = np.zeros((result_length, freedom))

    # 计算中间参数
    a_0 = 1 / (beta * delta_time ** 2)
    a_1 = gamma / (beta * delta_time)
    a_2 = 1 / (beta * delta_time)
    a_3 = 1 / (2 * beta)
    a_4 = gamma / beta
    a_5 = (a_4 - 2) * delta_time / 2
    temp_stiff = stiffness
    # 积分步迭代开始
    for i in range(result_length - 1):
        # 计算荷载增量
        delta_load = load[i + 1] - load[i]
        # 计算等效荷载
        equ_load = delta_load + mass @ (a_2 * vel[i] + a_3 * acc[i]) + damping @ (a_4 * vel[i] + a_5 * acc[i])
        # 计算等效刚度
        equ_stiffness = temp_stiff + a_0 * mass + a_1 * damping
        # newton-raphson迭代
        # 初始化迭代变量
        delta_dpm = dpm[i]  # 迭代初始位移
        delta_force = force[i]  # 迭代初始恢复力
        # print("Step :%d" % i)
        # is_dyn = True
        while True:
            # print("等效荷载:", equ_load)
            delta_temp_dpm = inv(equ_stiffness) @ equ_load  # 计算迭代步位移增量
            # print("位移增量:", delta_temp_dpm)
            temp_dpm = delta_dpm + delta_temp_dpm  # 求解迭代后位移
            # print("迭代后位移:", temp_dpm)
            # 计算当前恢复力与刚度
            temp_stiff, temp_force = bouc_mdof(temp_dpm, vel[i], models, np.array([i, i + 1]) * delta_time)
            # 计算等效残余力
            # print("上一步恢复力:", delta_force)
            # print("当前恢复力:", temp_force)
            # print("恢复力增量:", temp_force - delta_force)
            # print("当前刚度矩阵:", temp_stiff)
            # print("当前刚度矩阵乘位移增量:", temp_stiff @ delta_temp_dpm)
            # print("当前刚度矩阵乘位移:", temp_stiff @ temp_dpm)
            # print("按照等式右边计算出的等效荷载:", (temp_stiff + a_0 * mass + a_1 * damping) @ delta_temp_dpm)
            equ_load -= temp_force - delta_force  # 这一步恢复力减去上一步恢复力
            equ_stiffness = temp_stiff + a_0 * mass + a_1 * damping
            # if is_dyn:
            # 如果已经减过，则不需要再减
            # print("减去动力项")
            equ_load -= (a_0 * mass + a_1 * damping) @ delta_temp_dpm
            # is_dyn = False

            # 更新位移与恢复力
            delta_dpm = temp_dpm
            delta_force = temp_force
            # print(equ_load)
            if sum(abs(equ_load)) < 1e-5:
                break
        # 计算真实位移
        dpm[i + 1] = temp_dpm
        # 计算真实加速度
        acc[i + 1] = a_0 * (dpm[i + 1] - dpm[i]) - a_2 * vel[i] - (a_3 - 1) * acc[i]
        # 计算速度
        vel[i + 1] = a_1 * (dpm[i + 1] - dpm[i]) + (1 - a_4) * vel[i] + (1 - a_4 / 2) * acc[i] * delta_time
        force[i + 1] = temp_force

    return dpm, vel, acc, force
