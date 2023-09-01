if __name__ == "__main__":
    from libs.layer_shear import LayerShear
    from scipy.linalg import *
    import numpy as np
    from libs.vibration_mode import subspace_iteration

    # 验证子空间迭代法的程序
    m = np.array([2672, 2760, 2300])
    k = np.array([2.485, 1.921, 1.522]) * 1e4
    layer_shear = LayerShear(m, k)
    [omega, phi] = eig(layer_shear.k, layer_shear.m)
    omega = np.sort(np.sqrt(omega))
    omega_1, phi_1 = subspace_iteration(layer_shear.m, layer_shear.k, 3)
    print("频率精确解为:")
    print(omega.real)
    print("子空间迭代法解为:")
    print(omega_1)
