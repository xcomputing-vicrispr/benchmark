import pandas as pd

input_file = 'gw_data/CHOPCHOP320_360.tsv'
output_file = 'gw_data/CHOPCHOP320_360.csv'

df = pd.read_csv(input_file, sep='\t')

df.to_csv(output_file, index=False, encoding='utf-8')

print(f"Đã chuyển đổi thành công {input_file} sang {output_file}")