#!/usr/bin/env python3
"""
Étape 2: Prédire structures 3D via ESMFold API

Input:  data/RBPfasta/all_rbp.fasta
Output: data/RBPstructures/<protein_id>.pdb

Note: ~2-3h pour 93 séquences (API rate-limited)
"""
import requests
import time
from pathlib import Path
from Bio import SeqIO

# Config
INPUT_FILE = Path("data/RBPfasta/all_rbp.fasta")
OUTPUT_DIR = Path("data/RBPstructures")
DELAY = 3.0  # secondes entre requêtes
ESMFOLD_API = "https://api.esmatlas.com/foldSequence/v1/pdb/"


def predict_structure(sequence: str, protein_id: str, retry: int = 3) -> bool:
    output_file = OUTPUT_DIR / f"{protein_id}.pdb"
    
    if output_file.exists():
        print(f"  [SKIP] {protein_id} - already exists")
        return True
    
    for attempt in range(retry):
        try:
            response = requests.post(
                ESMFOLD_API,
                data=sequence,
                headers={"Content-Type": "text/plain"},
                timeout=300
            )
            
            if response.status_code == 200:
                with open(output_file, 'w') as f:
                    f.write(response.text)
                return True
            elif response.status_code == 429:
                wait = 60 * (attempt + 1)
                print(f"  [RATE LIMIT] Waiting {wait}s...")
                time.sleep(wait)
            else:
                print(f"  [ERROR] {protein_id}: HTTP {response.status_code}")
                return False
                
        except requests.exceptions.Timeout:
            print(f"  [TIMEOUT] {protein_id} - attempt {attempt + 1}")
            time.sleep(30)
        except Exception as e:
            print(f"  [ERROR] {protein_id}: {e}")
            return False
    
    return False


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    
    sequences = list(SeqIO.parse(INPUT_FILE, "fasta"))
    
    print(f"Input:  {INPUT_FILE}")
    print(f"Output: {OUTPUT_DIR}")
    print(f"Found {len(sequences)} sequences")
    print()
    
    success = 0
    failed = []
    
    for i, record in enumerate(sequences):
        protein_id = record.id.split("|")[0]
        seq_len = len(record.seq)
        
        print(f"[{i+1}/{len(sequences)}] {protein_id} ({seq_len} aa)")
        
        if predict_structure(str(record.seq), protein_id):
            success += 1
        else:
            failed.append(protein_id)
        
        time.sleep(DELAY)
    
    print(f"\n=== SUMMARY ===")
    print(f"Success: {success}/{len(sequences)}")
    if failed:
        print(f"Failed: {failed}")


if __name__ == "__main__":
    main()