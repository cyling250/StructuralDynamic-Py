import numpy as np
from libs.sbs_integration_nonlinear import Bilinear2, Bilinear3, Boucwen
from libs.figtools import hysteresis_figure, common
from scipy import integrate
from matplotlib import pyplot as plt


def bouc_wen(vel, z_last, time):
    """
    Bouc-Wen模型的等效刚度计算程序
    本程序仅针对单自由度的计算
    Parameters
    ----------
    dpm 当前位移
    vel 当前速度
    z_last 上一时刻滞变位移
    time [上一时刻，本时刻]
    state (刚度，刚度比，A，beta，gamma,n)

    Returns 恢复力
    -------

    """

    def function(z, t):
        """微分方程"""
        nonlocal vel
        A, beta, gamma, n = 1, 1.5, 0.5, 2  # boucwen模型的内置参数
        z_dot = A * vel - beta * abs(vel) * abs(z) ** (n - 1) * z - gamma * vel * abs(z) ** n
        return z_dot

    # 输出滞变位移
    z = integrate.odeint(function, z_last, time)
    return z[-1][0]


if __name__ == "__main__":
    step = 2000
    dpm = np.linspace(0, step / 40, step)
    dpm = np.sin(dpm) * np.linspace(0, 3, step)
    force = np.zeros(step)
    print("----------(1)----------")
    # 拟静力试验
    bilinear2 = Bilinear2(12.35, 0, 15, 0, 0, 0)
    for i in range(1, len(dpm)):
        force[i] = bilinear2(dpm[i], dpm[i] - dpm[i - 1])[1]
    common([np.linspace(0, step, step), dpm], ["Load step", "displacement/m"], save_file="6.2.2_1.svg")
    hysteresis_figure([dpm, force], save_file="6.2.2_2.svg")
    print("----------(2)----------")
    bilinear3 = Bilinear3(12.35, 5, 0, 10, 15, 0, 0, 0)
    for i in range(1, len(dpm)):
        force[i] = bilinear3(dpm[i], dpm[i] - dpm[i - 1])[1]
    hysteresis_figure([dpm, force], save_file="6.2.2_3.svg")
    print("----------(3)----------")
    dpm = np.linspace(0, step / 80, step)
    dpm = np.sin(dpm) * np.linspace(0, 10, step)
    boucwen = Boucwen(12.35, 0.2, 1, 0.5, 0.5, 3, 0, 0, 0, 0)
    for i in range(1, len(dpm)):
        force[i] = boucwen(dpm[i], (dpm[i] - dpm[i - 1]) / 0.02, np.array([i, i + 1]) * 0.02)[1]
    common([np.linspace(0, step, step), dpm], ["Load step", "displacement/m"], y_tick=3, save_file="6.2.2_4.svg")
    hysteresis_figure([dpm, force], x_tick=3.2, y_tick=12, save_file="6.2.2_5.svg")
    # plt.plot(dpm)
    # plt.show()
    # plt.plot(dpm, force)
    # plt.show()
