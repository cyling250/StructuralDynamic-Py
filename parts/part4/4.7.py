if __name__ == "__main__":
    import numpy as np
    from libs.layer_shear import LayerShear
    from scipy.linalg import eig
    from libs.vibration_mode import dunkerley

    m = np.array([2762, 2760, 2300])
    k = np.array([2.485, 1.921, 1.522]) * 1e4
    layer_shear = LayerShear(m, k)
    [omega, phi] = eig(layer_shear.k, layer_shear.m)
    omega = omega[-1] ** 0.5
    phi = phi[:, 0]
    omega_1 = dunkerley(layer_shear.m, layer_shear.k)
    print("omega精确值:%f" % omega.real)
    print("dunkerley方法计算值:%f" % omega_1)
