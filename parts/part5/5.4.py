if __name__ == "__main__":
    from libs.modal_analysis import modal_response_spectrum, cqc, srss, modal_mass
    from libs.layer_shear import LayerShear
    import numpy as np
    from scipy.linalg import eig

    layer_shear = LayerShear(np.array([2762, 2760, 2300]), np.array([2.485, 1.921, 1.522]) * 1e4)
    force_srss = modal_response_spectrum(layer_shear.m, layer_shear.k, 0.05, srss)
    force_cqc = modal_response_spectrum(layer_shear.m, layer_shear.k, 0.05, cqc)
    [omega, phi] = eig(layer_shear.k, layer_shear.m)
    sort_num = np.argsort(omega)
    omega = np.sqrt(omega[sort_num])
    phi = phi[sort_num]
    modal_mass1 = modal_mass(layer_shear.m, phi)
    print("振型参与质量系数：")
    print(modal_mass1)
    print("SRSS方法楼层剪力：")
    print(force_srss)
    print("CQC方法楼层剪力：")
    print(force_cqc)
