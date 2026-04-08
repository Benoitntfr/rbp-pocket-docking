#!/usr/bin/env python3
"""
Étape 4: Extraire features des pockets pour clustering

Features extraites par pocket:
    - score: P2Rank confidence score
    - probability: P2Rank probability
    - sas_points: solvent accessible surface points
    - surf_atoms: surface atoms count
    - center_x/y/z: pocket center coordinates
    - residue_count: number of residues
    - aa_composition: % of each amino acid type (20 features)
    - hydrophobic_ratio: % hydrophobic residues
    - charged_ratio: % charged residues
    - polar_ratio: % polar residues

Output: CSV avec une ligne par pocket
"""
import csv
import argparse
from pathlib import Path
from collections import Counter

# Amino acid properties
HYDROPHOBIC = set('AILMFVPWG')
CHARGED = set('DEKRH')
POLAR = set('STYNQC')
ALL_AA = 'ACDEFGHIKLMNPQRSTVWY'

def parse_predictions_csv(csv_path: Path) -> list:
    """Parse P2Rank predictions CSV"""
    pockets = []
    with open(csv_path, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            pockets.append({
                'rank': int(row['rank']),
                'score': float(row['score']),
                'probability': float(row['probability']),
                'sas_points': int(row['sas_points']),
                'surf_atoms': int(row['surf_atoms']),
                'center_x': float(row['center_x']),
                'center_y': float(row['center_y']),
                'center_z': float(row['center_z']),
            })
    return pockets

def parse_residues_csv(csv_path: Path) -> dict:
    """Parse P2Rank residues CSV, group by pocket"""
    pockets_residues = {}
    with open(csv_path, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            pocket_rank = int(row.get('pocket', row.get('   pocket', '0')).strip())
            residue = row.get('residue_label', row.get('   residue_label', '')).strip()
            
            if pocket_rank not in pockets_residues:
                pockets_residues[pocket_rank] = []
            
            # Extract single-letter AA from residue label (e.g., "ALA_123" -> "A")
            if residue:
                aa_3letter = residue.split('_')[0] if '_' in residue else residue[:3]
                aa_1letter = AA_3TO1.get(aa_3letter.upper(), 'X')
                pockets_residues[pocket_rank].append(aa_1letter)
    
    return pockets_residues

# 3-letter to 1-letter amino acid mapping
AA_3TO1 = {
    'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',
    'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
    'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
    'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'
}

def compute_aa_features(residues: list) -> dict:
    """Compute amino acid composition features"""
    if not residues:
        return {f'aa_{aa}': 0.0 for aa in ALL_AA} | {
            'hydrophobic_ratio': 0.0,
            'charged_ratio': 0.0,
            'polar_ratio': 0.0,
            'residue_count': 0
        }
    
    counts = Counter(residues)
    total = len(residues)
    
    features = {f'aa_{aa}': counts.get(aa, 0) / total for aa in ALL_AA}
    features['hydrophobic_ratio'] = sum(counts.get(aa, 0) for aa in HYDROPHOBIC) / total
    features['charged_ratio'] = sum(counts.get(aa, 0) for aa in CHARGED) / total
    features['polar_ratio'] = sum(counts.get(aa, 0) for aa in POLAR) / total
    features['residue_count'] = total
    
    return features

def extract_features(pockets_dir: str, output_file: str, min_score: float = 0.0, top_n: int = None):
    """Extract features from all pockets"""
    pockets_path = Path(pockets_dir)
    
    all_features = []
    
    # Find all prediction files
    for protein_dir in pockets_path.iterdir():
        if not protein_dir.is_dir():
            continue
        
        protein_id = protein_dir.name
        predictions_file = list(protein_dir.glob("*_predictions.csv"))
        residues_file = list(protein_dir.glob("*_residues.csv"))
        
        if not predictions_file:
            print(f"[WARN] No predictions for {protein_id}")
            continue
        
        pockets = parse_predictions_csv(predictions_file[0])
        residues_by_pocket = parse_residues_csv(residues_file[0]) if residues_file else {}
        
        # Filter and limit pockets
        pockets = [p for p in pockets if p['score'] >= min_score]
        if top_n:
            pockets = pockets[:top_n]
        
        for pocket in pockets:
            rank = pocket['rank']
            residues = residues_by_pocket.get(rank, [])
            aa_features = compute_aa_features(residues)
            
            features = {
                'protein_id': protein_id,
                'pocket_rank': rank,
                **pocket,
                **aa_features
            }
            all_features.append(features)
    
    # Write to CSV
    if all_features:
        output_path = Path(output_file)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        fieldnames = list(all_features[0].keys())
        with open(output_path, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(all_features)
        
        print(f"Extracted {len(all_features)} pockets from {len(set(f['protein_id'] for f in all_features))} proteins")
        print(f"Output: {output_file}")
    else:
        print("No pockets found!")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, help="P2Rank output directory")
    parser.add_argument("--output", default="data/features/pocket_features.csv")
    parser.add_argument("--min-score", type=float, default=0.0, help="Minimum pocket score")
    parser.add_argument("--top-n", type=int, default=3, help="Top N pockets per protein")
    args = parser.parse_args()
    
    extract_features(args.input, args.output, args.min_score, args.top_n)
