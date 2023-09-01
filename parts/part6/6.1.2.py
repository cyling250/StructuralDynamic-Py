if __name__ == "__main__":
    import numpy as np
    from libs.sbs_integration_linear import segmented_parsing, duhamel_numerical
    from libs.figtools import response_figure

    x = np.arange(0, 40, 0.02)
    dpm_segmented = segmented_parsing(1, 1, np.sin(x), 0.02)[0]
    dpm_duhamel = duhamel_numerical(1, 1, np.sin(x), 0.02)

    response_figure([dpm_duhamel, dpm_segmented],
                    [["Duhamel", "#0080FF", "-"], ["Segmented", "#000000", "--"]],
                    x_tick=8, y_tick=10, x_length=48,
                    delta_time=0.02, save_file="../res/6.1.2.svg"
                    )
