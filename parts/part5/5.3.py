if __name__ == "__main__":
    from libs.modal_analysis import pse_spectrum, fourier
    from libs.sbs_integration_linear import segmented_parsing
    from libs.figtools import spectrum_figure
    from libs.read_wave import read_quake_wave, length_normal, pga_normal

    # 因为工程中使用伪加速度谱比较广泛，这里仅展示伪加速度谱的计算
    quake_wave, delta_time = read_quake_wave("../../res/RSN88_SFERN_FSD172.AT2")
    quake_wave_2 = length_normal(quake_wave, 2)
    quake_wave_2 = pga_normal(quake_wave_2, 0.35)
    quake_wave = pga_normal(quake_wave, 0.35)

    spectrum_dpm_1, spectrum_vel_1, spectrum_acc_1 = pse_spectrum(quake_wave, delta_time, fourier)
    spectrum_dpm_2, spectrum_vel_2, spectrum_acc_2 = pse_spectrum(quake_wave, delta_time,
                                                                  segmented_parsing)

    spectrum_figure(spectrum_acc_1, "acc/m·s$^{-2}$", y_length=2, save_file="fourier_acc_spectrum.svg")
    spectrum_figure(spectrum_vel_1, "vel/m·s$^{-1}$", y_length=0.1, save_file="fourier_vel_spectrum.svg")
    spectrum_figure(spectrum_dpm_1, "dpm/m", y_length=0.08, save_file="fourier_dpm_spectrum.svg")
    spectrum_figure(spectrum_acc_2, "acc/m·s$^{-2}$", y_length=2, save_file="seg_acc_spectrum.svg")
    spectrum_figure(spectrum_vel_2, "vel/m·s$^{-1}$", y_length=0.1, save_file="seg_vel_spectrum.svg")
    spectrum_figure(spectrum_dpm_2, "dpm/m", y_length=0.08, save_file="seg_dpm_spectrum.svg")
