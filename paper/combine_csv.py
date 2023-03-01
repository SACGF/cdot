#!/bin/env python3

columns = ['date', 'details', 'server', 'hgvs', 'valid_hgvs', 'can_resolve']
df_list = []
for filename in glob.glob("validate*.csv"):
    df = pd.read_csv(filename)
    df = df[columns]
    df_list.append(df)
df_combined.sort_values("date").to_csv("combo.csv", index=False)


