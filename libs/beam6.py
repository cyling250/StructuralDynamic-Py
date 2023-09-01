"""
程序描述：这是单根等截面六自由度梁单元程序。
程序封装了单个等截面六自由度梁单元的相关参数与操作
"""
import numpy as np
from libs.damping import rayleigh


class Beam6:
    """
    扩展的Beam6单元类型
    """

    def __init__(self,
                 node_1,
                 node_2,
                 freedom,
                 elastic=1,
                 inertia=1,
                 area=1,
                 rho=1):
        self.nodes = (node_1, node_2)  # 梁的方向始终由node_1指向node_2
        self.elastic = elastic  # 弹性模量
        self.inertia = inertia  # 截面惯性矩
        self.area = area  # 横截面积
        self.rho = rho  # 密度
        self.beam_length = self.get_length()  # 长度
        self.freedom = freedom  # 用False代表自由，用True代表约束
        # 提供两种质量矩阵，集中质量矩阵和一致质量矩阵
        self.loca_m = np.diag([0.5, 0.5, 0, 0.5, 0.5, 0]) * self.rho * self.area * self.beam_length  # 此为集中质量矩阵
        self.loca_m_ = self.get_loca_m_()  # 此为一致质量矩阵
        self.loca_k = self.get_loca_k()  # 局部单元刚度矩阵
        self.transpose = self.get_transpose()  # 坐标变换矩阵
        self.k = self.transpose.T @ self.loca_k @ self.transpose  # 整体坐标下的单元刚度矩阵
        self.m = self.transpose.T @ self.loca_m @ self.transpose  # 整体坐标下的集中质量矩阵
        self.m_ = self.transpose.T @ self.loca_m_ @ self.transpose  # 整体坐标下的一致质量矩阵

    def get_length(self):
        """计算单元长度"""
        length = (self.nodes[0][0] - self.nodes[1][0]) ** 2 + (
                self.nodes[0][1] - self.nodes[1][1]) ** 2
        return length ** 0.5

    def get_loca_k(self):
        """组装局部坐标下的单元刚度矩阵"""
        L = self.beam_length
        E = self.elastic
        I = self.inertia
        A = self.area
        k = np.array([
            [E * A / L, 0, 0, -E * A / L, 0, 0],
            [0, 12 * E * I / L ** 3, 6 * E * I / L ** 2, 0, -12 * E * I / L ** 3, 6 * E * I / L ** 2],
            [0, 6 * E * I / L ** 2, 4 * E * I / L, 0, -6 * E * I / L ** 2, 2 * E * I / L],
            [-E * A / L, 0, 0, E * A / L, 0, 0],
            [0, -12 * E * I / L ** 3, -6 * E * I / L ** 2, 0, 12 * E * I / L ** 3, -6 * E * I / L ** 2],
            [0, 6 * E * I / L ** 2, 2 * E * I / L, 0, -6 * E * I / L ** 2, 4 * E * I / L]
        ])
        return k

    def get_loca_m_(self):
        """组装局部坐标下的一致质量矩阵"""
        rho = self.rho
        A = self.area
        L = self.beam_length
        m_ = np.array([
            [140, 0, 0, 70, 0, 0],
            [0, 156, 22 * L, 0, 54, -13 * L],
            [0, 22 * L, 4 * L ** 2, 0, 13 * L, -3 * L ** 2],
            [70, 0, 0, 140, 0, 0],
            [0, 54, 13 * L, 0, 156, -22 * L],
            [0, -13 * L, -3 * L ** 2, 0, -22 * L, 4 * L ** 2]
        ]) * (rho * A * L / 420)
        return m_

    def get_transpose(self):
        """计算坐标变换矩阵"""
        cosa = (self.nodes[0][0] - self.nodes[1][0]) / self.beam_length
        sina = (self.nodes[0][1] - self.nodes[1][1]) / self.beam_length
        transpose = np.array([
            [cosa, sina, 0, 0, 0, 0],
            [-sina, cosa, 0, 0, 0, 0],
            [0, 0, 1, 0, 0, 0],
            [0, 0, 0, cosa, sina, 0],
            [0, 0, 0, -sina, cosa, 0],
            [0, 0, 0, 0, 0, 1]])
        return transpose


class Structure:
    """结构类"""

    def __init__(self):
        self.nodes = []
        self.elements = []  # 单元列表，下标就是单元编号
        self.K = []  # 整体刚度矩阵
        self.M = []  # 整体质量矩阵
        self.C = []  # 整体阻尼矩阵

    def add_element(self, node_1, node_2):
        """添加单元"""
        if node_1 not in self.nodes:
            self.nodes.append(node_1)
        if node_2 not in self.nodes:
            self.nodes.append(node_2)
        node_id_1 = self.nodes.index(node_1)
        node_id_2 = self.nodes.index(node_2)
        freedom = np.array([node_id_1 * 3, node_id_1 * 3 + 1, node_id_1 * 3 + 2,
                            node_id_2 * 3, node_id_2 * 3 + 1, node_id_2 * 3 + 2])
        self.elements.append(Beam6(node_1, node_2, freedom))

    def add_restrain(self, restrain):
        """添加约束"""
        self.M = np.delete(self.M, restrain, axis=0)
        self.M = np.delete(self.M, restrain, axis=1)
        self.K = np.delete(self.K, restrain, axis=0)
        self.K = np.delete(self.K, restrain, axis=1)

    def get_K_C(self):
        """求解整体刚度矩阵、整体质量矩阵"""
        size = len(self.nodes)
        self.K = np.array(np.zeros((size * 3, size * 3)))
        self.M = np.array(np.zeros((size * 3, size * 3)))
        for element in self.elements:
            for i in range(6):
                for j in range(6):
                    self.K[element.freedom[i], element.freedom[j]] += element.k[i, j]
                    self.M[element.freedom[i], element.freedom[j]] += element.m[i, j]
        return

    def get_C(self):
        self.C = rayleigh(self.M, self.K, 0.05)
