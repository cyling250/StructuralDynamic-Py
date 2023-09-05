from sklearn.model_selection import train_test_split, cross_val_score
from sklearn.datasets import load_iris
from sklearn.linear_model import LogisticRegression
# 加载数据集
iris = load_iris()
X = iris.data
y = iris.target
# 划分数据集为训练集和测试集
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=7)
# 定义模型
model = LogisticRegression()
# 进行交叉验证
cv_scores = cross_val_score(model, X_train, y_train, cv=10)
# 输出交叉验证的平均得分
print("交叉验证平均得分：", cv_scores.mean())
