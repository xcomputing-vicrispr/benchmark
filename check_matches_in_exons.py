import pandas as pd

df_reject = pd.read_csv('mm10dbRejected.tsv', sep='\t', header=None, names=['sequence_tsv', 'reason'])
filtered_reject = df_reject[df_reject['reason'].str.strip() == "Multiple matches in exons"]

df_data = pd.read_csv('data500k_filter.csv')
df_data['seq_20_csv'] = df_data['sequence'].str[:20]

#gop
result = pd.merge(
    df_data, 
    filtered_reject, 
    left_on='seq_20_csv', 
    right_on='sequence_tsv'
)

final_output = result[['sequence', 'reason']]


final_output.to_csv('merged_results.csv', index=False)