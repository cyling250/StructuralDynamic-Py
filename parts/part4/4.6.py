if __name__ == "__main__":
    from libs.layer_shear import LayerShear
    import numpy as np
    from scipy.linalg import eig
    from libs.vibration_mode import lanczos

    m = np.array([2762, 2760, 2300])
    k = np.array([2.485, 1.921, 1.522]) * 1e4
    layer_shear = LayerShear(m, k)
    [omega, phi] = eig(layer_shear.k, layer_shear.m)
    omega = np.sqrt(omega)
    omega_1, phi_1 = lanczos(layer_shear.m, layer_shear.k, 3)
    print("频率精确值为:")
    print(omega.real)
    print("频率计算值为:")
    print(1 / np.sqrt(omega_1))
    print("振型计算值为:")
    print(phi_1)
