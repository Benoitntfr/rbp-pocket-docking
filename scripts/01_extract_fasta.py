#!/usr/bin/env python3
"""
Étape 1: Extraire séquences FASTA depuis RBPbase.csv

Input:  data/RBPbase/RBPbase.csv
Output: data/RBPfasta/all_rbp.fasta
"""
import csv
from pathlib import Path

# Config
INPUT_FILE = Path("data/RBPbase/RBPbase.csv")
OUTPUT_FILE = Path("data/RBPfasta/all_rbp.fasta")
MAX_LENGTH = 600  # ESMFold limit


def extract_fasta():
    OUTPUT_FILE.parent.mkdir(parents=True, exist_ok=True)
    
    stats = {'total': 0, 'kept': 0, 'filtered': 0}
    
    with open(INPUT_FILE, 'r') as f_in, open(OUTPUT_FILE, 'w') as f_out:
        reader = csv.DictReader(f_in)
        for row in reader:
            stats['total'] += 1
            seq = row['protein_sequence']
            seq_len = len(seq)
            
            if seq_len > MAX_LENGTH:
                stats['filtered'] += 1
                continue
            
            protein_id = row['protein_ID']
            phage_id = row['phage_ID']
            
            f_out.write(f">{protein_id}|{phage_id}|{seq_len}\n")
            for i in range(0, seq_len, 80):
                f_out.write(seq[i:i+80] + "\n")
            
            stats['kept'] += 1
    
    print(f"Input:  {INPUT_FILE}")
    print(f"Output: {OUTPUT_FILE}")
    print(f"Total: {stats['total']}")
    print(f"Kept: {stats['kept']}")
    print(f"Filtered (>{MAX_LENGTH} aa): {stats['filtered']}")


if __name__ == "__main__":
    extract_fasta()