## Ligand Clustering for Bias Reduction and Robust ML Splits

This project provides a pipeline for clustering chemical ligands based on structural similarity. The motivation is to avoid bias in downstream analysis and ensure robust train/test splits when applying machine learning methods.

**Problem Statement**

In cheminformatics and molecular modeling, many ligands belong to the same chemical series (congeners). If members of the same cluster end up in both training and test datasets, this may lead to information leakage, model overfitting, and misleading performance metrics.

To mitigate this:

  * We cluster similar ligands using fingerprint-based similarity.

  * We then ensure that no cluster is split between train and test.

**Project Structure**
`````
├── ecfp_cluster_split.ipynb     # Main notebook: clustering & splitting
├── init_data.csv                # Input ligand data
├── src/
│   ├── clustering.py            # Clustering functions
│   └── visualization.py         # Plotting functions
├── pics/
│   └── clustering.pdf           # Cluster visualization (see below)
├── LICENSE
└── README.md
`````

**Notebook**

The full clustering workflow and group-aware data splitting are demonstrated in the Jupyter notebook: 📘 ecfp_cluster_split.ipynb

It includes:

    - Molecular fingerprint generation (ECFP)

    - Butina clustering algorithm, which is a form of hierarchical clustering specifically designed for chemical fingerprints using a distance (or similarity) threshold

    - Visualization of clustering results

    - An example of GroupKFold split that keeps clusters intact


**Visualization**

Cluster assignments are visualized in the file below: 📎 clustering.pdf

The PDF is large and detailed — open it separately to explore all clusters and ligand assignments.

**Key Features**

    - Group-aware cross-validation (GroupKFold)

    - Easy-to-extend modular structure (src/)

    - Clean visualization utilities

    - Reproducible and readable Jupyter pipeline

**Getting Started**

    Install required libraries (e.g. RDKit, scikit-learn, pandas, matplotlib) and run the notebook:

    jupyter notebook ecfp_cluster_split.ipynb

**License**

    This project is licensed under the MIT License.
