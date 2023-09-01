"""
程序说明：这是层剪切模型的程序，也称葫芦串模型。
程序封装了层剪切模型的基本参数和相关操作。
"""

import numpy as np


class LayerShear:
    def __init__(self, each_shear_m, each_shear_k, h=3):
        """
        层剪切模型
        Parameters
        ----------
        each_shear_m 各层质量，一维数组，长度为层数
        each_shear_k 各层刚度，一维数组
        h 结构高度，浮点数
        """
        self.each_shear_k = each_shear_k
        self.each_shear_m = each_shear_m
        self.h = h

        self.freedom = len(each_shear_k)  # 计算当前结构自由度
        self.m = np.diag(each_shear_m)  # 构造质量矩阵，即对角化each_shear_m
        self.k = self.get_k(each_shear_k)  # 构造刚度矩阵

    def get_k(self, each_shear_k):
        """
        刚度矩阵构造函数
        Parameters
        ----------
        each_shear_k 各层刚度，一维数组

        Returns 刚度矩阵，二维数组
        -------

        """
        # 全零初始化刚度矩阵
        k = np.zeros((self.freedom, self.freedom))

        for i in range(self.freedom - 1):
            k[i, i] = each_shear_k[i] + each_shear_k[i + 1]
            k[i, i + 1] = -each_shear_k[i + 1]
            k[i + 1, i] = -each_shear_k[i + 1]
        k[self.freedom - 1, self.freedom - 1] = each_shear_k[self.freedom - 1]

        return k

    def get_p_delta(self):
        """
        P-delta矩阵构造函数
        Parameters
        ----------

        Returns 考虑p-delta效应的刚度矩阵，二维数组
        -------

        """
        # 初始化
        k_p = np.zeros((self.freedom, self.freedom))
        w = np.zeros(self.freedom)

        for i in range(self.freedom):
            w[i] = np.sum(self.each_shear_m[:self.freedom - i]) * 9.8

        # 开始构造p_delta矩阵
        for i in range(self.freedom - 1):
            k_p[i, i] = w[i] + w[i + 1]
            k_p[i, i + 1] = -w[i + 1]
            k_p[i + 1, i] = -w[i + 1]
        k_p[self.freedom - 1, self.freedom - 1] = w[self.freedom - 1]

        return k_p / self.h


def get_k(each_shear_k):
    """
    摆脱类的刚度矩阵构造函数
    Parameters
    ----------
    each_shear_k 各层刚度，一维数组

    Returns 刚度矩阵，二维数组
    -------

    """
    # 初始化
    freedom = len(each_shear_k)
    k = np.zeros((freedom, freedom))

    # 开始构造
    for i in range(freedom - 1):
        k[i, i] = each_shear_k[i] + each_shear_k[i + 1]
        k[i, i + 1] = -each_shear_k[i + 1]
        k[i + 1, i] = -each_shear_k[i + 1]
    k[freedom - 1, freedom - 1] = each_shear_k[freedom - 1]

    return k


def get_each_shear_k(k):
    """
    通过层剪切模型的刚度矩阵，提取各层刚度
    Parameters
    ----------
    k 刚度矩阵，二维数组

    Returns 各层刚度，一维数组
    -------

    """
    # 初始化
    freedom = len(k)
    each_shear_k = np.zeros(freedom)  # 初始化层刚度
    each_shear_k[-1] = k[-1, -1]

    for i in range(freedom - 2, -1, -1):
        each_shear_k[i] = k[i, i] - each_shear_k[i + 1]

    return each_shear_k
