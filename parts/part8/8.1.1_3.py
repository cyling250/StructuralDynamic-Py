from sklearn import datasets
import xgboost as xgb
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error
from sklearn.metrics import mean_absolute_error
from sklearn.metrics import r2_score
import numpy as np
import pandas as pd
dataset=pd.read_csv("res/boston_housing_data.csv")
dataset=dataset.dropna()
data = dataset.iloc[:, 0:13]
target = dataset.iloc[:, 13]
X_train, X_test, y_train, y_test = train_test_split(data, target, test_size=0.3, random_state=3)  # 划分训练集和测试集
model = xgb.XGBRegressor()  # 搭建机器学习算法svm
model.fit(X_train, y_train)  # 训练svm
y_pred = model.predict(X_test)  # 预测结果
print("平均绝对误差：", mean_absolute_error(y_test, y_pred))
print("均方误差：", mean_squared_error(y_test, y_pred))
print("决定系数：", r2_score(y_test, y_pred))
results=np.stack([y_test, y_pred])
np.savetxt('results/8.1.1_3_pred_results.csv',results,delimiter=',')
