from sklearn import datasets
from sklearn import svm
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score
from sklearn.metrics import precision_score
from sklearn.metrics import recall_score
iris = datasets.load_iris()  # 导入鸢尾花数据集
iris_X = iris.data
iris_y = iris.target
X_train, X_test, y_train, y_test = train_test_split(iris_X, iris_y, test_size=0.3, random_state=1)  # 划分训练集和测试集
model = svm.SVC()  # 搭建机器学习算法svm
model.fit(X_train, y_train)  # 训练svm
y_pred = model.predict(X_test)  # 预测结果
print("预测准确率：", accuracy_score(y_test, y_pred))
print("预测查准率：", precision_score(y_test, y_pred,average='macro'))
print("预测召回率：", recall_score(y_test, y_pred,average='macro'))
