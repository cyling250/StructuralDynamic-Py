if __name__ == "__main__":
    from libs.damping import *
    from libs.layer_shear import LayerShear

    print("----------(1)----------")
    m = np.array([2762, 2760, 2300])
    k = np.array([2.485, 1.921, 1.522]) * 1e4
    layer_shear = LayerShear(m, k)
    c_r = rayleigh(layer_shear.m, layer_shear.k, 0.05)
    c_c = caughey(layer_shear.m, layer_shear.k, 0.05, [0, 1])
    print("Rayleigh阻尼矩阵为:")
    print(c_r)
    print("Caughey阻尼矩阵为:")
    print(c_c)

    print("----------(2)----------")
    c_n = damping_nonclassical(layer_shear.m, layer_shear.k, k, [0.05, 0.02], [1, 2])
    print("非经典阻尼矩阵为:")
    print(c_n)
