if __name__ == "__main__":
    from libs.read_wave import read_quake_wave, pga_normal, length_normal
    from libs.layer_shear import LayerShear
    import numpy as np
    from libs.sbs_integration_linear import center_difference_multiple
    from libs.modal_analysis import modal_superposition
    from libs.figtools import response_figure

    quake, delta_time = read_quake_wave("../../res/RSN88_SFERN_FSD172.AT2")  # 读取地震波
    quake = pga_normal(quake, 0.35)  # 地震波PGA设置
    quake = length_normal(quake, 2)  # 设置地震波长度为2倍

    layer_shear = LayerShear(np.array([2762, 2760, 2300]), np.array([2.485, 1.921, 1.522]) * 1e4)
    load = []
    for i in range(len(quake)):
        load.append(np.array([0, 0, quake[i] * 2300]))
    center_dpm, center_vel, center_acc = center_difference_multiple(layer_shear.m, layer_shear.k, load, delta_time,
                                                                    0.05, np.array([0, 0, 0]), np.array([0, 0, 0]),
                                                                    len(quake))
    dpm = modal_superposition(layer_shear.m, layer_shear.k, load, delta_time)
    response_figure([center_dpm[:, 0], dpm[:, 0]],
                    [["Central difference", "#0080FF", "-"], ["Fourier", "#000000", "--"]],
                    x_tick=8, y_tick=0.008, x_length=88,
                    delta_time=0.005, save_file="../res/5.1_5.svg"
                    )
    response_figure([center_dpm[:, 1], dpm[:, 1]],
                    [["Central difference", "#0080FF", "-"], ["Fourier", "#000000", "--"]],
                    x_tick=8, y_tick=0.015, x_length=88,
                    delta_time=0.005, save_file="../res/5.1_6.svg"
                    )
    response_figure([center_dpm[:, 2], dpm[:, 2]],
                    [["Central difference", "#0080FF", "-"], ["Fourier", "#000000", "--"]],
                    x_tick=8, y_tick=0.015, x_length=88,
                    delta_time=0.005, save_file="../res/5.1_7.svg"
                    )
