if __name__ == "__main__":
    from libs.modal_analysis import fourier, duhamel_numerical
    from libs.read_wave import read_quake_wave, pga_normal, length_normal
    from libs.sbs_integration_linear import center_difference_single
    from libs.figtools import response_figure

    quake, delta_time = read_quake_wave("../../res/RSN88_SFERN_FSD172.AT2")  # 读取地震波
    quake = pga_normal(quake, 0.35)  # 地震波PGA设置
    quake = length_normal(quake, 1)  # 设置地震波长度.末端补零法

    dpm = fourier(2762, 24850, quake * 2762, delta_time, 0.05)
    dpm_duhamel = duhamel_numerical(2762, 24850, quake * 2762, delta_time, 0.05)
    dpm_center, vel_center, acc_center = center_difference_single(2762, 24850, quake * 2762, delta_time, 0.05,
                                                                  result_length=len(quake))
    response_figure([dpm_center, dpm_duhamel],
                    [["Duhamel", "#0080FF", "-"], ["Central difference", "#000000", "--"]],
                    x_tick=8, y_tick=0.005,
                    delta_time=0.005, save_file="../res/5.1_1.svg")
    response_figure([dpm_center - dpm_duhamel],
                    [["Error", "#0080FF", "-"]],
                    x_tick=8, y_tick=2e-6,
                    delta_time=0.005, save_file="../res/5.1_2.svg")
    response_figure([dpm_center, dpm],
                    [["Central difference", "#0080FF", "-"], ["Fourier", "#000000", "--"]],
                    x_tick=16, y_tick=0.005, x_length=128,
                    delta_time=0.005, save_file="../res/5.1_3.svg")
    response_figure([dpm_center - dpm],
                    [["Error", "#0080FF", "-"]],
                    x_tick=16, y_tick=3e-6, x_length=128,
                    delta_time=0.005, save_file="../res/5.1_4.svg")
