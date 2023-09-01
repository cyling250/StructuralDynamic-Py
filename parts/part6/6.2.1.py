import numpy as np
from libs.layer_shear import LayerShear
from libs.read_wave import read_quake_wave, pga_normal
from libs.sbs_integration_linear import newmark_beta_multiple
from libs.figtools import response_figure

if __name__ == "__main__":
    """
    几何非线性计算程序
    该程序一般只做演示使用
    """
    # 参数初始化
    m = np.array([600 for i in range(8)])  # 定义各层质量
    k = np.array([1e5 for i in range(8)])  # 定义各层刚度
    layer_shear = LayerShear(m, k)  # 创建层剪切模型
    temp = np.array([0 for i in range(8)])  # 初始化其余参数
    # 定义荷载
    at2, delta_time = read_quake_wave("../../res/ELCENTRO.DAT")
    at2 = pga_normal(at2, 2)
    load = []
    for i in range(len(at2)):
        load.append(np.array([0, 0, 0, 0, 0, 0, 0, at2[i]]) * 600)
    # 计算结构响应
    new_dpm, new_vel, new_acc = newmark_beta_multiple(layer_shear.m, layer_shear.k, load,
                                                      delta_time, 0.02,
                                                      temp, temp, temp)

    new_dpm_nlr, new_vel_nlr, new_acc_nlr = newmark_beta_multiple(layer_shear.m,
                                                                  layer_shear.k - layer_shear.get_p_delta(), load,
                                                                  delta_time, 0.02,
                                                                  temp, temp, temp)
    response_figure([new_dpm[:, -1], new_dpm_nlr[:, -1]],
                    [["Linear", "#0080FF", "-"], ["P-delta", "#000000", "--"]],
                    x_tick=8, y_tick=0.05, x_length=40,
                    delta_time=delta_time, save_file="../res/6.2.1.svg"
                    )
