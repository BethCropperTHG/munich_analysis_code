import glob
import pandas as pd
import numpy as np

fs = np.sort(glob.glob("*.pkl"))

dfs = []

for f in fs:
    print(f)
    df = pd.read_pickle(f)
    df2 = df[["EASSIGN", "AREA", "sAREA"]]
    print(df2)


