import sys, math
import pandas as pd
from sklearn.metrics import precision_score, recall_score, f1_score, confusion_matrix
from experiment_functions_500 import fold_rna


#Doench
GC_MIN = 25
GC_MAX = 75
#MFE_THRESHOLD = -24.5
#MFE_THRESHOLD = -26
#MFE_THRESHOLD = -22.5
MFE_THRESHOLD = -22
scaffold_default = "GUUUUAGAGCUAGAAAUAGCAAGUUAAAAUAAGGCUAGUCCGUUAUCAACUUGAAAAAGUGGCACCGAGUCGGUGCUUUU"


def extract_30mer_with_ngg(file_path):
    
    df = pd.read_csv(file_path, skiprows=1)# tai dong dau tien ghi supplementary
    matches = []
    column_name = 'Expanded Sequence'
    if column_name in df.columns:
        matches = df[column_name].tolist()            
    return matches

def process_sgrna_data(file_path, output_path):

    df = pd.read_csv(file_path, skiprows=1)
    df['Gene % Rank'] = pd.to_numeric(df['Gene % Rank'], errors='coerce')
    
    def label_top_20(group):
        group = group.sort_values(by='Gene % Rank', ascending=False)
        
        num_rows = len(group)
        n_top = math.ceil(num_rows * 0.2)
        
        group['Is_Top_20'] = [1] * n_top + [0] * (num_rows - n_top)
        
        return group

    final_df = df.groupby('Transcript', group_keys=False).apply(label_top_20)
    
    final_df.to_csv(output_path, index=False)
    print(f"File saved at: {output_path}")



def filter_and_save_csv(results, output_file="./vicrispr_doench_all.csv"):

    df1 = pd.read_csv("./Doench_data_processed.csv")
    y_true_list = df1.iloc[:, -1].astype(int).tolist()
    y_true_rank = df1.iloc[:, -3].astype(float).tolist()

    data_rows = []
    all_sgrnas = []
    
    for i in range(len(results)):
        sgrna = results[i]
        x = sgrna[4:23]
        gc_count = x.count('G') + x.count('C')
        gc_content = (gc_count / len(x)) * 100 if len(x) > 0 else 0


        ss, mfe_score = fold_rna(x + scaffold_default)

        is_valid = (
            GC_MIN <= gc_content <= GC_MAX and
            mfe_score >= MFE_THRESHOLD
        )

        if is_valid:
            row = {
                "seq": x,
                "flank_sequence": results[i],
                "GC_Content": round(gc_content, 2),
                "MFE Score": mfe_score,
            }
            data_rows.append(row)

        row = {
            "Sequence": x,
            "GC_Content": gc_content,
            "MFE Score": mfe_score,
            "Vicrispr_rate": is_valid,
            "Doench_rate": y_true_list[i],
            "Gene % rank": y_true_rank[i]
        }

        all_sgrnas.append(row)

    df = pd.DataFrame(data_rows)

    df.to_csv(output_file, index=False, encoding='utf-8-sig')

    dff = pd.DataFrame(all_sgrnas)
    dff.to_csv(output_file, index=False, encoding='utf-8-sig')


def calculate_detailed_metrics(file1_path, file2_path):
    df1 = pd.read_csv(file1_path)
    df2 = pd.read_csv(file2_path)

    y_true = df1.iloc[:, -1].astype(int) 
    y_pred = df2.iloc[:, -3].astype(int)

    tn, fp, fn, tp = confusion_matrix(y_true, y_pred).ravel()

    precision = precision_score(y_true, y_pred)
    recall = recall_score(y_true, y_pred)
    f1 = f1_score(y_true, y_pred)

    npv = tn / (tn + fn) if (tn + fn) > 0 else 0

    print("="*40)
    print("DETAILED PERFORMANCE STATISTICS")
    print("="*40)
    print(f"Total samples compared:    {len(y_true)}")
    print("-" * 40)
    print(f"Ground Truth counts (y_true):")
    print(f"  - Positive (1):          {sum(y_true)}")
    print(f"  - Negative (0):          {len(y_true) - sum(y_true)}")
    print("-" * 40)
    print(f"Prediction counts (y_pred):")
    print(f"  - Predicted Positive (1): {sum(y_pred)}")
    print(f"  - Predicted Negative (0): {len(y_pred) - sum(y_pred)}")
    print("-" * 40)
    print(f"Confusion Matrix breakdown:")
    print(f"  - True Positives (TP):   {tp} (Pred 1, Actual 1)")
    print(f"  - True Negatives (TN):   {tn} (Pred 0, Actual 0)")
    print(f"  - False Positives (FP):  {fp} (Pred 1, Actual 0)")
    print(f"  - False Negatives (FN):  {fn} (Pred 0, Actual 1)")
    print("-" * 40)
    print(f"  - Total Predicted Pos:   {tp + fp}")
    print(f"  - Precision:             {precision:.4f}")
    print(f"  - Recall:                {recall:.4f}")
    print(f"  - F1-score:              {f1:.4f}")
    print(f"  - NPV (Neg. Pred. Value): {npv:.4f}")
    print("="*40)




results = extract_30mer_with_ngg("./Doench_data.csv")
process_sgrna_data('Doench_data.csv', 'Doench_data_processed.csv')
print(f"Found {len(results)} sequences.")
filter_and_save_csv(results)
calculate_detailed_metrics("./Doench_data_processed.csv", "./vicrispr_doench_all.csv")