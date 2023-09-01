import numpy as np
from libs.sbs_integration_nonlinear import newmark_bili, center_bili, newmark_bouc, center_bouc
from libs.read_wave import read_quake_wave, pga_normal, read_data_sap2000
from libs.layer_shear import LayerShear
from matplotlib import pyplot as plt
from libs.figtools import response_figure, hysteresis_figure

if __name__ == "__main__":
    print("-----------bili----------")
    quake, delta_time = read_quake_wave("../../res/ELCENTRO.DAT")
    m = np.array([2762, 2760, 2300])
    k = np.array([2.485, 1.921, 1.522]) * 1e4
    # m = [41000]
    # k = [15000000]
    quake = pga_normal(quake, 2.2)
    quake[0] = 0  # 动力增量方法的荷载第一步一定必须为0
    freedom = 3
    layer_shear = LayerShear(m, k)
    load = []
    for i in range(len(quake)):
        load.append(np.array([quake[i] * 2762, 0, 0, ]))
    init = np.array([0, 0, 0])
    bili_parms = []
    for i in range(freedom):
        bili_parms.append((k[i], 0 * k[i], 9e2, 0, 0, 0))
    bouc_parms = []
    for i in range(freedom):
        bouc_parms.append((k[i], 0, 0, 0, 0, 0.4, 1, 0.5, 0.5, 3))
    new_dpm, new_vel, new_acc, new_force = newmark_bili(layer_shear.m, layer_shear.k,
                                                        load, delta_time, 0.05,
                                                        init, init, init, bili_parms)
    trans_mat = np.diag([1 for i in range(freedom)]) + np.diag([-1 for i in range(freedom - 1)], 1)
    trans_mat_plus = np.diag([1 for i in range(freedom)]) + np.diag([-1 for i in range(freedom - 1)], -1)
    dpm = new_dpm @ trans_mat
    force = new_force @ np.linalg.inv(trans_mat_plus)
    cen_dpm, cen_vel, cen_acc, cen_force = center_bili(layer_shear.m, layer_shear.k,
                                                       load, delta_time, 0.05,
                                                       init, init, bili_parms)
    dpm1 = cen_dpm @ trans_mat
    force1 = cen_force @ np.linalg.inv(trans_mat_plus)
    # -----------------------------------Bili----------------------------------
    dpm_sap, force_sap = read_data_sap2000("../../res/SAP2000验证/1.txt")
    response_figure([new_dpm[:, 0], dpm_sap], [["Newmark", "#0080FF", "-"], ["SAP2000", "#000000", "--"]],
                    x_tick=8, y_tick=0.03, x_length=40,
                    delta_time=delta_time, save_file="sap_dpm.svg")
    response_figure([force[:, 0], force_sap], [["Newmark", "#0080FF", "-"], ["SAP2000", "#000000", "--"]],
                    x_tick=8, y_tick=500, x_length=40,
                    delta_time=delta_time, save_file="sap_force.svg")
    hysteresis_figure([dpm_sap, force_sap], x_tick=0.02, y_tick=4e2, save_file="sap_hts.svg")
    # -----------------------------------Bouc----------------------------------
    # plt.plot(cen_dpm[:, 0][0:500])
    # plt.plot(dpm[:, 0])
    # plt.legend()
    # plt.show()
    # plt.plot(cen_dpm[:, 0], cen_force[:, 0])
    # plt.show()
