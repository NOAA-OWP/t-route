import pandas as pd

df = pd.read_csv("qlats_croton.csv", index_col=0)

print(df)

df1 = df.iloc[:, 0:12]

df2 = df.iloc[:, 12:26]

df1.to_csv("qlats_croton_loop1.csv")

df2.to_csv("qlats_croton_loop2.csv")
