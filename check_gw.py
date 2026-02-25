import pandas as pd

df = pd.read_csv('./data/gw_ecolik12.csv') 
coords = df['sgRNA_loc'].str.split(':', expand=True)[1].str.split('-', expand=True)

df['loc_start'] = pd.to_numeric(coords[0])
df['loc_end'] = pd.to_numeric(coords[1])

range_start, range_end = 2001, 420000
# range_start, range_end = 80001, 120000
# range_start, range_end = 160001, 200000
# range_start, range_end = 240001, 280000
# range_start, range_end = 320001, 360000

filtered_df = df[(df['loc_start'] >= range_start) & (df['loc_end'] <= range_end)]
count = len(filtered_df)
print(count)
