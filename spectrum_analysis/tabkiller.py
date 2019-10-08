import pandas as pd

fname = '116Cd_d,p_117Cd_10degrees_bite_3_peaks.txt'

df = pd.read_table(fname, index_col = False)
df.to_csv(fname, sep = ' ', index = False)
