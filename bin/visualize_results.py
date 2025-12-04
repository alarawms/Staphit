#!/usr/bin/env python3

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import sys
import os

def generate_visualizations(summary_csv, output_dir):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    if summary_csv.endswith('.tsv'):
        df = pd.read_csv(summary_csv, sep='\t')
    else:
        df = pd.read_csv(summary_csv)
    
    # 1. MLST Distribution
    if 'mlst_scheme' in df.columns and 'mlst_sequence_type' in df.columns:
        plt.figure(figsize=(10, 6))
        sns.countplot(x='mlst_sequence_type', data=df)
        plt.title('MLST Sequence Type Distribution')
        plt.savefig(f"{output_dir}/mlst_distribution.png")
        plt.close()

    # 2. AMR Genes Heatmap (Binary)
    # Assuming AMR genes are in columns or we need to parse them
    # This depends on the aggregator output format. 
    # For now, we'll look for columns that might represent AMR classes or genes
    # Or if there's a column 'amr_genes' with comma-separated values
    
    if 'amrfinder_genes' in df.columns:
        # Expand the list
        amr_series = df['amrfinder_genes'].fillna('').astype(str).str.split(';').apply(lambda x: x if isinstance(x, list) else [])
        all_genes = sorted(list(set([item for sublist in amr_series for item in sublist])))
        
        if all_genes:
            matrix = pd.DataFrame(0, index=df['sample_id'], columns=all_genes)
            for idx, row in df.iterrows():
                genes = row['amrfinder_genes'].split(';') if isinstance(row['amrfinder_genes'], str) else []
                for gene in genes:
                    if gene in all_genes:
                        matrix.loc[row['sample_id'], gene] = 1
            
            plt.figure(figsize=(12, 8))
            sns.heatmap(matrix, cmap="Blues", cbar=False, linewidths=.5)
            plt.title('AMR Gene Presence/Absence')
            plt.tight_layout()
            plt.savefig(f"{output_dir}/amr_heatmap.png")
            plt.close()

    # 3. Virulence Genes Heatmap
    if 'virulence_genes' in df.columns:
         # Expand the list
        vir_series = df['virulence_genes'].fillna('').astype(str).str.split(';').apply(lambda x: x if isinstance(x, list) else [])
        all_vir = sorted(list(set([item for sublist in vir_series for item in sublist])))
        
        if all_vir:
            matrix = pd.DataFrame(0, index=df['sample_id'], columns=all_vir)
            for idx, row in df.iterrows():
                genes = row['virulence_genes'].split(';') if isinstance(row['virulence_genes'], str) else []
                for gene in genes:
                    if gene in all_vir:
                        matrix.loc[row['sample_id'], gene] = 1
            
            plt.figure(figsize=(12, 8))
            sns.heatmap(matrix, cmap="Reds", cbar=False, linewidths=.5)
            plt.title('Virulence Gene Presence/Absence')
            plt.tight_layout()
            plt.savefig(f"{output_dir}/virulence_heatmap.png")
            plt.close()

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: visualize_results.py <summary.csv> <output_dir>")
        sys.exit(1)
    
    summary_csv = sys.argv[1]
    output_dir = sys.argv[2]
    
    generate_visualizations(summary_csv, output_dir)
