import pandas as pd
import numpy as np
df = pd.DataFrame({'A': [1, 2, np.nan], 'B': [5, np.nan, np.nan], 'C': [1, 2, 3]})
new_df1 = df.dropna(axis=1)  # 删除列
print("删除存在NaN值的列")
print(new_df1)
new_df2 = df.dropna(axis=0)  # 删除行
print("删除存在NaN值的行")
print(new_df2)
new_df3 = df.dropna(thresh=2)  # 删除含有少于2个非NaN值的行
print("删除含有少于2个非NaN值的行")
print(new_df3)
