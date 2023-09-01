if __name__ == "__main__":
    from libs.figtools import wave_figure
    from libs.read_wave import read_quake_wave

    at2 = read_quake_wave("../../res/RSN88_SFERN_FSD172.AT2")
    dt2 = read_quake_wave("../../res/RSN88_SFERN_FSD172.DT2")
    vt2 = read_quake_wave("../../res/RSN88_SFERN_FSD172.VT2")
    # wave_figure函数为自定义绘图函数，在libs文件夹的figtools.py文件中
    wave_figure(at2, "acc/g", save_file="acc.svg")
    wave_figure(dt2, "dpm/cm", y_tick=6, save_file="dpm.svg")
    wave_figure(vt2, r"vel/cm·s$^{-1}$", y_tick=6, save_file="vel.svg")
