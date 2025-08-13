#!/usr/bin/env python3
"""
DiRT Candidate Index Gene Finder

This script identifies the top candidate index genes for each target gene
based on normalized dispersion (std/mean) of expression ratios
in control samples.

Author: Your Name
Date: 2025-08-13
"""

import pandas as pd
import numpy as np

# --- Load data ---
# Test.csv must have:
#   - First column header: "Geneid"
#   - CPM values for 10,000 genes x 18 samples (C1..C9 = control, T1..T9 = treated)
df = pd.read_csv("Test.csv")

# Avoid divide-by-zero in ratios
EPS = 1e-12

def tencand(idx,
             control_start="C1", control_end="C9",
             all_start="C1", all_end="T9",
             top_n=10):
    """
    Identify top_n candidate index genes for a target gene at row index `idx`,
    ranked by normalized dispersion (std/mean) of control ratios.
    """
    target_geneid = df.at[idx, 'Geneid']
    target_cdata = df.loc[idx, control_start:control_end]
    target_all_data = df.loc[idx, all_start:all_end]

    # Normalized dispersion across controls
    def ndiv(x):
        ratios = target_cdata / (x[control_start:control_end] + EPS)
        return np.std(ratios, ddof=1) / np.mean(ratios)

    # Full ratio across all columns
    def ratio_all(x):
        return target_all_data / (x[all_start:all_end] + EPS)

    # Candidate IDs as "target/candidate"
    def id_pair(x):
        return f"{target_geneid}/{x['Geneid']}"

    # Build DataFrame
    mdf = pd.DataFrame({
        "ID": df.apply(id_pair, axis=1),
        "ndiv": df.apply(ndiv, axis=1)
    })

    # Append ratio columns
    ratio_df = df.apply(ratio_all, axis=1)
    mdf = pd.concat([mdf, ratio_df], axis=1)

    # Remove self-pair
    mdf = mdf[~mdf["ID"].str.endswith(f"/{target_geneid}")]

    # Sort and return top_n
    return mdf.sort_values(by="ndiv", ascending=True).head(top_n)

if __name__ == "__main__":
    # --- Test run on first 10 genes ---
    df_C = pd.concat([tencand(i) for i in range(10)], axis=0, ignore_index=True)
    df_C.to_csv("DiRT_test.csv", index=False)

    print("Wrote: DiRT_test.csv")
    print(df_C.head())

    # --- Notes for full run ---
    # Full run on 10,000 genes will be slow (~24h).
    # Split across cores:
    #   Half 1: indices [0, 5000]
    # df_1 = pd.concat([tencand(i) for i in range(0, 5000)], axis=0, ignore_index=True)
    # df_1.to_csv('DiRT_test1.csv', index=False)
    #
    #   Half 2: indices [5000, 10000]
    # df_2 = pd.concat([tencand(i) for i in range(5000, 10000)], axis=0, ignore_index=True)
    # df_2.to_csv('DiRT_test2.csv', index=False)
    # Merge:
    # pd.concat([df_1, df_2], axis=0, ignore_index=True).to_csv('DiRT_full.csv', index=False)
