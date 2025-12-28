
import pandas as pd
df = pd.read_csv('data/processed/liver_proteome.csv')
min_tpm = df['liver_tpm'].min()
count_below_1 = len(df[df['liver_tpm'] < 1.0])
print(f"Min TPM: {min_tpm}")
print(f"Values < 1.0: {count_below_1}")
