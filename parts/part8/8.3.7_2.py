import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split
from sklearn.model_selection import cross_val_score, GridSearchCV, RandomizedSearchCV
from sklearn.metrics import mean_squared_error, r2_score, mean_absolute_error
from sklearn.preprocessing import MinMaxScaler
import matplotlib.pyplot as plt


def read_data(path):
    dataset = pd.read_csv(path, encoding='utf-8')
    data = dataset.iloc[:, 0:9]
    data = pd.get_dummies(data)

    target = dataset.iloc[:, 9]
    train_X, test_X, train_y, test_y = train_test_split(data, target, test_size=0.3, random_state=7)
    return train_X, test_X, train_y, test_y


model = RandomForestRegressor()
path = "res/data_norm_8.3.csv"
train_X, test_X, train_y, test_y = read_data(path)
params = {
    'n_estimators': [100, 200, 300, 400, 500],
    'max_depth': [3, 5, 7, 9, 11],
    'min_samples_split': [2,3, 4,5, 6,7, 8,9, 10],
    'min_samples_leaf': [1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
}
# 随机搜索迭代1000轮调优
random_search = RandomizedSearchCV(model, params, cv=10, scoring='neg_mean_squared_error', n_iter=100, n_jobs=-1)
random_search.fit(train_X, train_y)
# 获取最佳超参数组合和对应的模型权重系数
best_params = random_search.best_params_
best_estimator = random_search.best_estimator_
# 输出最佳超参数组合
print(best_params)
pred_y = best_estimator.predict(test_X)
# 性能评估
mse = mean_squared_error(test_y, pred_y)
r2 = r2_score(test_y, pred_y)
mae = mean_absolute_error(test_y, pred_y)
print('预测MSE值为：',mse)
print('预测R2值为：',r2)
print('预测MAE值为：',mae)
results1 = np.stack([test_y, pred_y])
# 输出预测结果
np.savetxt('results/rf_pred_results.csv', results1, delimiter=',')
# 绘制结果
plt.figure(figsize=(10, 10), dpi=600)
plt.scatter(test_y, pred_y)
plt.show()