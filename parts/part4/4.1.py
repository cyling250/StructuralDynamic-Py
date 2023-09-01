if __name__ == "__main__":
    from libs.layer_shear import LayerShear
    from scipy.linalg import *
    import numpy as np
    from libs.vibration_mode import rayleigh_psi

    # 验证rayleigh法的程序
    m = np.array([2762, 2760, 2300])
    k = np.array([2.485, 1.921, 1.522]) * 1e4
    layer_shear = LayerShear(m, k)
    [omega, phi] = eig(layer_shear.k, layer_shear.m)
    omega = omega ** 0.5
    phi1 = phi[:, -1].T
    # 使用了精确振型来计算，得到的结果应该完全一致
    omega1 = rayleigh_psi(
        layer_shear.m,
        layer_shear.k,
        phi1
    )  # 可以在这里调整假设频率的值
    print("精确频率:")
    print(omega)
    print("精确振型:")
    print(phi)
    print("计算得到的基频:")
    print(omega1)
    print("相对误差:{:.2%}".format((omega[-1].real - omega1) / omega[-1].real))
