import pandas as pd
data = {'name': ['Ana', 'Bob', 'Tom', 'Ana', 'David'],
        'age': [25, 30, 35, 25, 40]}
df = pd.DataFrame(data)
# 使用duplicated()函数查找重复值
duplicates = df.duplicated()
# 打印重复值的位置
print(duplicates)
drop_duplicates = df.drop_duplicates()
# 打印新的数据框
print(drop_duplicates)
