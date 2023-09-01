if __name__ == "__main__":
    import numpy as np
    from libs.layer_shear import LayerShear
    from libs.beam6 import Structure

    print("# ----------(1)----------")
    m = np.array([2762, 2760, 2300])
    k = np.array([2.485, 1.921, 1.522]) * 1e4
    layer_shear = LayerShear(m, k)
    print(layer_shear.m)
    print(layer_shear.k)

    print("# ----------(2)----------")
    structure = Structure()  # 创建结构对象
    L = 1
    structure.add_element((0, 0), (0, 0.5 * L))  # 添加单元
    structure.add_element((0, 0.5 * L), (0, L))  # 添加单元
    structure.get_K_C()  # 组装未经过约束的刚度矩阵
    structure.add_restrain([0, 1, 2])
    structure.get_C()
    print("整体刚度矩阵为:")
    print(structure.K)
    print("整体质量矩阵为:")
    print(structure.M)
    print("整体阻尼矩阵为:")
    print(structure.C)
