# RBP Pocket Docking

Matching des pockets de binding des RBP (Receptor Binding Proteins) de phages avec les pockets de récepteurs bactériens pour prédire les interactions phage-bactérie.

## Thèse

L'état de l'art ([PhageHostLearn](https://github.com/dimiboeckaerts/PhageHostLearn), Nature Comms 2024) prédit l'interaction phage-bactérie au niveau protéine entière:

```
RBP entière → ESM-2 embedding → XGBoost → P(infection)
```

**Limite**: Le binding se fait sur quelques résidus (pocket), pas sur toute la protéine. Signal dilué.

**Notre approche**: Prédire au niveau pocket:

```
RBP → Structure 3D → Pocket RBP ─┐
                                 ├→ Matching → P(infection)
Récepteur → Structure 3D → Pocket récepteur ─┘
```

## Pipeline

| Étape | Script | Input | Output |
|-------|--------|-------|--------|
| 1 | `01_extract_fasta.py` | `RBPbase.csv` | `all_rbp.fasta` |
| 2 | `02_predict_structures.py` | FASTA | PDB (ESMFold) |
| 3 | `03_predict_pockets.py` | PDB | Pockets (P2Rank) |
| 4 | `04_extract_features.py` | Pockets | `pocket_features.csv` |
| 5 | `05_umap_cluster.py` | Features | UMAP + clusters |

## Installation

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

## Usage

```bash
python scripts/01_extract_fasta.py
python scripts/02_predict_structures.py  # ~2h (API calls)
python scripts/03_predict_pockets.py     # Requires P2Rank
python scripts/04_extract_features.py
python scripts/05_umap_cluster.py
```

## Data

- `RBPbase.csv`: 274 RBP annotées (source: PhageHostLearn)
- Structures 3D: prédites via [ESMFold API](https://esmatlas.com/)
- Pockets: détectés via [P2Rank](https://github.com/rdk/p2rank)

## Roadmap

- [x] Extraction FASTA
- [ ] Prédiction structures RBP
- [ ] Détection pockets RBP
- [ ] Clustering pockets RBP
- [ ] Structures récepteurs bactériens
- [ ] Matching pocket × pocket

## Références

- PhageHostLearn: [Nature Comms 2024](https://doi.org/10.1038/s41467-024-48675-6)
- ESMFold: [Lin et al. 2023](https://doi.org/10.1126/science.ade2574)
- P2Rank: [Krivák & Hoksza 2018](https://doi.org/10.1186/s13321-018-0285-8)