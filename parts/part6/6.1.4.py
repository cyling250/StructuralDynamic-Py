if __name__ == "__main__":
    from libs.sbs_integration_linear import newmark_beta_multiple, center_difference_multiple
    from libs.figtools import response_figure
    from libs.layer_shear import LayerShear
    from libs.read_wave import read_quake_wave, pga_normal, length_normal
    import numpy as np

    quake, delta_time = read_quake_wave("../../res/RSN88_SFERN_FSD172.AT2")  # 读取地震波
    quake = pga_normal(quake, 0.35)  # 地震波PGA设置
    quake = length_normal(quake, 2)  # 设置地震波长度为1倍

    layer_shear = LayerShear(np.array([2762, 2760, 2300]), np.array([2.485, 1.921, 1.522]) * 1e4)
    load = []
    for i in range(len(quake)):
        load.append(np.array([0, 0, quake[i] * 2300]))
    center_dpm, center_vel, center_acc = center_difference_multiple(layer_shear.m, layer_shear.k, load, delta_time,
                                                                    0.05, np.array([0, 0, 0]), np.array([0, 0, 0]),
                                                                    len(quake))
    newmark_dpm, newmark_vel, newmark_acc = newmark_beta_multiple(layer_shear.m, layer_shear.k, load, delta_time,
                                                               0.05, np.array([0, 0, 0]), np.array([0, 0, 0]),
                                                               np.array([0, 0, 0]), result_length=len(quake))
    response_figure([center_dpm[:, 0], newmark_dpm[:, 0]],
                    [["Center difference", "#0080FF", "-"], ["Newmark beta", "#000000", "--"]],
                    x_tick=8, y_tick=0.008, x_length=88,
                    delta_time=0.005, save_file="../res/6.1.4_1.svg"
                    )
    response_figure([center_dpm[:, 1], newmark_dpm[:, 1]],
                    [["Center difference", "#0080FF", "-"], ["Newmark beta", "#000000", "--"]],
                    x_tick=8, y_tick=0.015, x_length=88,
                    delta_time=0.005, save_file="../res/6.1.4_2.svg"
                    )
    response_figure([center_dpm[:, 2], newmark_dpm[:, 2]],
                    [["Center difference", "#0080FF", "-"], ["Newmark beta", "#000000", "--"]],
                    x_tick=8, y_tick=0.015, x_length=88,
                    delta_time=0.005, save_file="../res/6.1.4_3.svg"
                    )
