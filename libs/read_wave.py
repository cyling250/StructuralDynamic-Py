import re
import numpy as np


def read_quake_wave(file_name):
    """
    地震波读取函数
    :return 时程序列，时间间隔
    """
    with open(file_name, 'r') as fp:
        lines = fp.readlines()
    delta_time = eval(re.findall(r"\.[0-9]*", lines[3])[0])  # 提取采样间隔
    quake_wave = []
    for i in range(4, len(lines)):
        temp = lines[i].split()  # 读取elcentro波
        # temp = re.findall(r"[-\.]+[0-9]+E[-+][0-9]+", lines[i])  # 提取加速度数值
        for j in temp:
            quake_wave.append(eval(j))
    return np.array(quake_wave), delta_time


def pga_normal(quake_wave, acc):
    """
    PGA调幅函数
    :param quake_wave:地震波数组
    :param acc: PGA
    :return:
    """
    max_acc = max(quake_wave)
    ratio = acc / max_acc
    return quake_wave * ratio


def length_normal(quake_wave, length, flag=True):
    """
    长度调整函数，为了进行地震波长度的调整
    Parameters
    ----------
    quake_wave 地震波
    length 目标长度，支持两种方式的输入，比例类型或者长度类型
    flag 标志位，当flag为True时，表示length为长度类型；当flag为False时，表示length为比例类型

    Returns 经过长度调整后的地震波
    -------

    """
    if flag:
        length = int(len(quake_wave) * length)
        new_quake = np.zeros(length)
        for i in range(min(length, len(quake_wave))):
            new_quake[i] = quake_wave[i]
    else:
        new_quake = np.zeros(length)
        for i in range(min(length, len(quake_wave))):
            new_quake[i] = quake_wave[i]
    return new_quake


def write_csv(quake, delta_time, file_name):
    with open(file_name, "w") as fp:
        for i in range(len(quake)):
            fp.write(str(delta_time * i) + "," + str(quake[i]) + "\n")
    return


def read_data_sap2000(file_name):
    with open(file_name, "r") as fp:
        lines = fp.readlines()
    dpm = []
    force = []
    for i in range(15,len(lines)):
        dpm.append(float(lines[i].split()[1]))
        force.append(float(lines[i].split()[2]))
    dpm = np.array(dpm)
    force = np.array(force)
    return dpm, force
