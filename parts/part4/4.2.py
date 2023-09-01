if __name__ == "__main__":
    from libs.layer_shear import LayerShear
    import numpy as np
    from scipy.linalg import *
    from libs.vibration_mode import rayleigh_ritz

    m = np.array([2762, 2760, 2300])
    k = np.array([2.485, 1.921, 1.522]) * 1e4
    layer_shear = LayerShear(m, k)
    [omega, phi] = eig(layer_shear.k, layer_shear.m)
    phi_1 = np.array([[0.5, 1.0, 1.5], [-1.0, -1.0, 1.0]]).T  # 这里的振型是任意给出的
    omega_1, phi_1 = rayleigh_ritz(layer_shear.m, layer_shear.k, phi_1)
    print("精确频率为:")
    print(np.sqrt(omega))
    print("频率为:")
    print(omega_1)
    print("振型为:")
    print(phi_1)
