import pandas as pd
import numpy as np
df = pd.DataFrame({'A': [1, 2, np.nan, 4], 'B': [5, np.nan, np.nan, 1], 'C': [1, 2, 3, 2]})
print("填充整个数据框")
new_df1 = df.fillna(value=0)
print(new_df1)
new_df2 = df.fillna(method='ffill')  # 使用前一个值填充
print("使用前一个值填充")
print(new_df2)
new_df3 = df.fillna(method='bfill')  # 使用后一个值填充
print("使用后一个值填充")
print(new_df3)
