if __name__ == "__main__":
    from libs.vibration_mode import *
    from libs.layer_shear import *

    # 验证矩阵迭代法的程序
    m = np.array([2762, 2760, 2300])
    k = np.array([2.485, 1.921, 1.522]) * 1e4
    layer_shear = LayerShear(m, k)
    [omega, phi] = eig(layer_shear.k, layer_shear.m)
    omega = omega ** 0.5
    phi = phi[:, 2]
    omega_1, phi_1 = mat_iterate_base(layer_shear.m, layer_shear.k)
    omega_2, phi_2 = mat_iterate_high(layer_shear.m, layer_shear.k, [omega_1], [phi_1], 10e-4)
    omega_3, phi_3 = mat_iterate_highest(layer_shear.m, layer_shear.k)
    print("omega精确值:%f,%f,%f" % (omega[0].real, omega[1].real, omega[2].real))
    print("omega迭代值:%f,%f,%f" % (omega_3, omega_2, omega_1))
