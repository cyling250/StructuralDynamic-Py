if __name__ == "__main__":
    from libs.vibration_mode import load_depended_ritz_vector
    import numpy as np
    from libs.layer_shear import LayerShear

    m = np.array([2762, 2760, 2300])
    k = np.array([2.485, 1.921, 1.522]) * 1e4
    layer_shear = LayerShear(each_shear_m=m, each_shear_k=k)
    load = np.array([1 for i in range(3)]).T  # 模拟全横向均布1的荷载
    count = 3
    ritz, alpha, beta = load_depended_ritz_vector(layer_shear.m, layer_shear.k, count=count, load=load)
    print("关于质量矩阵正交的验证:")
    print(ritz.T @ layer_shear.m @ ritz)
