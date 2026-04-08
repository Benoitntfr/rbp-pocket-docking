#!/usr/bin/env python3
"""
Étape 3: Prédire pockets de binding via P2Rank

Prérequis:
    - Java 8+ installé
    - P2Rank téléchargé: https://github.com/rdk/p2rank/releases

Usage:
    python 03_predict_pockets.py --input data/structures --output data/pockets --p2rank /path/to/p2rank

Output:
    - data/pockets/<protein_id>/
        - <protein_id>.pdb_predictions.csv  (liste des pockets avec scores)
        - <protein_id>.pdb_residues.csv     (résidus de chaque pocket)
        - visualizations/                    (fichiers PyMOL)
"""
import subprocess
import argparse
from pathlib import Path
import os
import shutil

def check_java():
    """Vérifie que Java est installé"""
    try:
        result = subprocess.run(["java", "-version"], capture_output=True, text=True)
        return result.returncode == 0
    except FileNotFoundError:
        return False

def run_p2rank(p2rank_dir: str, pdb_file: Path, output_dir: Path):
    """Exécute P2Rank sur un fichier PDB"""
    p2rank_script = Path(p2rank_dir) / "prank"
    
    if not p2rank_script.exists():
        # Windows
        p2rank_script = Path(p2rank_dir) / "prank.bat"
    
    cmd = [
        str(p2rank_script),
        "predict",
        "-f", str(pdb_file),
        "-o", str(output_dir)
    ]
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    return result.returncode == 0, result.stderr

def main(structures_dir: str, output_dir: str, p2rank_dir: str):
    if not check_java():
        print("ERROR: Java not found. Install Java 8+")
        return
    
    structures_path = Path(structures_dir)
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    pdb_files = list(structures_path.glob("*.pdb"))
    print(f"Found {len(pdb_files)} PDB files")
    
    success = 0
    failed = []
    
    for i, pdb_file in enumerate(pdb_files):
        protein_id = pdb_file.stem
        print(f"[{i+1}/{len(pdb_files)}] {protein_id}")
        
        protein_output = output_path / protein_id
        
        ok, error = run_p2rank(p2rank_dir, pdb_file, protein_output)
        
        if ok:
            success += 1
        else:
            failed.append(protein_id)
            print(f"  [ERROR] {error[:200] if error else 'Unknown'}")
    
    print(f"\n=== SUMMARY ===")
    print(f"Success: {success}/{len(pdb_files)}")
    if failed:
        print(f"Failed: {failed}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, help="Directory with PDB files")
    parser.add_argument("--output", default="data/pockets")
    parser.add_argument("--p2rank", required=True, help="Path to P2Rank directory")
    args = parser.parse_args()
    
    main(args.input, args.output, args.p2rank)
