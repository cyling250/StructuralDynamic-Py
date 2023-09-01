if __name__ == "__main__":
    import numpy as np
    from libs.sbs_integration_linear import duhamel_parse
    from libs.sbs_integration_linear import duhamel_numerical
    from libs.figtools import response_figure

    load = np.sin(np.arange(0, 40, 0.02))
    dpm_parse = duhamel_parse(1, 1, np.sin, 0.02, result_length=48/0.02)
    dpm_numerical = duhamel_numerical(1, 1, load, 0.02)

    response_figure([dpm_parse, dpm_numerical],
                    [["Parse", "#0080FF", "-"], ["Numerical", "#000000", "--"]],
                    x_tick=8, y_tick=10, x_length=48,
                    delta_time=0.02, save_file="../res/6.1.1.svg"
                    )

