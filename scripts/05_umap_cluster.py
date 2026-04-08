#!/usr/bin/env python3
"""
Étape 5: UMAP + HDBSCAN clustering des pockets RBP

Produit:
    - outputs/umap_pockets.html (interactif Plotly)
    - outputs/umap_pockets.png (statique)
    - outputs/cluster_stats.csv (stats par cluster)

Usage:
    python 05_umap_cluster.py --input data/features/pocket_features.csv --output outputs/

Dépendances:
    pip install umap-learn hdbscan plotly pandas scikit-learn matplotlib
"""
import argparse
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def main(features_file: str, output_dir: str, n_neighbors: int = 15, min_dist: float = 0.1, min_cluster_size: int = 5):
    try:
        import umap
        import hdbscan
        import plotly.express as px
        import plotly.graph_objects as go
    except ImportError as e:
        print(f"Missing dependency: {e}")
        print("Install: pip install umap-learn hdbscan plotly")
        return
    
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    # Load features
    df = pd.read_csv(features_file)
    print(f"Loaded {len(df)} pockets from {df['protein_id'].nunique()} proteins")
    
    # Select numeric features for UMAP
    feature_cols = [c for c in df.columns if c.startswith('aa_') or c in [
        'score', 'probability', 'sas_points', 'surf_atoms',
        'hydrophobic_ratio', 'charged_ratio', 'polar_ratio', 'residue_count'
    ]]
    
    X = df[feature_cols].values
    
    # Normalize
    from sklearn.preprocessing import StandardScaler
    X_scaled = StandardScaler().fit_transform(X)
    
    # UMAP
    print("Running UMAP...")
    reducer = umap.UMAP(
        n_neighbors=n_neighbors,
        min_dist=min_dist,
        n_components=2,
        random_state=42,
        metric='euclidean'
    )
    embedding = reducer.fit_transform(X_scaled)
    
    df['umap_x'] = embedding[:, 0]
    df['umap_y'] = embedding[:, 1]
    
    # HDBSCAN clustering
    print("Running HDBSCAN...")
    clusterer = hdbscan.HDBSCAN(
        min_cluster_size=min_cluster_size,
        min_samples=3,
        metric='euclidean'
    )
    df['cluster'] = clusterer.fit_predict(embedding)
    
    n_clusters = df['cluster'].nunique() - (1 if -1 in df['cluster'].values else 0)
    n_noise = (df['cluster'] == -1).sum()
    print(f"Found {n_clusters} clusters, {n_noise} noise points")
    
    # Interactive Plotly
    fig = px.scatter(
        df,
        x='umap_x',
        y='umap_y',
        color='cluster',
        hover_data=['protein_id', 'pocket_rank', 'score', 'residue_count', 
                   'hydrophobic_ratio', 'charged_ratio'],
        title=f'UMAP des pockets RBP ({n_clusters} clusters)',
        labels={'umap_x': 'UMAP 1', 'umap_y': 'UMAP 2'},
        color_continuous_scale='viridis'
    )
    fig.update_traces(marker=dict(size=8, opacity=0.7))
    fig.write_html(output_path / "umap_pockets.html")
    print(f"Saved: {output_path / 'umap_pockets.html'}")
    
    # Static matplotlib
    plt.figure(figsize=(12, 8))
    scatter = plt.scatter(
        df['umap_x'], df['umap_y'],
        c=df['cluster'],
        cmap='tab20',
        s=50,
        alpha=0.7
    )
    plt.colorbar(scatter, label='Cluster')
    plt.xlabel('UMAP 1')
    plt.ylabel('UMAP 2')
    plt.title(f'UMAP des pockets RBP - {len(df)} pockets, {n_clusters} clusters')
    plt.tight_layout()
    plt.savefig(output_path / "umap_pockets.png", dpi=150)
    print(f"Saved: {output_path / 'umap_pockets.png'}")
    
    # Cluster stats
    cluster_stats = df.groupby('cluster').agg({
        'protein_id': 'nunique',
        'score': 'mean',
        'residue_count': 'mean',
        'hydrophobic_ratio': 'mean',
        'charged_ratio': 'mean',
        'polar_ratio': 'mean',
        'sas_points': 'mean'
    }).round(3)
    cluster_stats.columns = ['n_proteins', 'avg_score', 'avg_residues', 
                            'avg_hydrophobic', 'avg_charged', 'avg_polar', 'avg_sas']
    cluster_stats.to_csv(output_path / "cluster_stats.csv")
    print(f"Saved: {output_path / 'cluster_stats.csv'}")
    
    # Print summary
    print("\n=== CLUSTER SUMMARY ===")
    print(cluster_stats.to_string())
    
    # Save full results
    df.to_csv(output_path / "pockets_with_clusters.csv", index=False)
    print(f"\nSaved: {output_path / 'pockets_with_clusters.csv'}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", default="outputs/")
    parser.add_argument("--n-neighbors", type=int, default=15)
    parser.add_argument("--min-dist", type=float, default=0.1)
    parser.add_argument("--min-cluster-size", type=int, default=5)
    args = parser.parse_args()
    
    main(args.input, args.output, args.n_neighbors, args.min_dist, args.min_cluster_size)
