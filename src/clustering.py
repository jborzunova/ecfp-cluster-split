import rdkit
from rdkit.Chem import rdFingerprintGenerator
from rdkit import Chem
from rdkit import DataStructs
from rdkit.ML.Cluster import Butina
from itertools import combinations, product


def tanimoto_distance_matrix(fp_list):
    """Calculate distance matrix for fingerprint list"""
    dissimilarity_matrix = []
    # Notice how we are deliberately skipping the first and last items in the list
    # because we don't need to compare them against themselves
    for i in range(1, len(fp_list)):
        # Compare the current fingerprint against all the previous ones in the list
        similarities = DataStructs.BulkTanimotoSimilarity(fp_list[i], fp_list[:i])
        # Since we need a distance matrix, calculate 1-x for every element in similarity matrix
        dissimilarity_matrix.extend([1 - x for x in similarities])
    return dissimilarity_matrix


def cluster_fingerprints(fingerprints, cutoff=0.2):
    """Cluster fingerprints
    Parameters:
        fingerprints
        cutoff: threshold for the clustering
    """
    # Calculate Tanimoto distance matrix
    distance_matrix = tanimoto_distance_matrix(fingerprints)
    # Now cluster the data with the implemented Butina algorithm:
    clusters = Butina.ClusterData(distance_matrix, len(fingerprints), cutoff, isDistData=True)
    clusters = sorted(clusters, key=len, reverse=True)
    return clusters


def compute_avg_similarities(clusters, fingerprints):
    intra_sims = []
    inter_sims = []

    # Intra-cluster similarity
    for cluster in clusters:
        if len(cluster) < 2:
            continue  # Can't compare within a cluster of 1 item
        for i, j in combinations(cluster, 2):
            sim = DataStructs.TanimotoSimilarity(fingerprints[i], fingerprints[j])
            intra_sims.append(sim)

    # Inter-cluster similarity
    for c1_idx in range(len(clusters)):
        for c2_idx in range(c1_idx + 1, len(clusters)):
            cluster1 = clusters[c1_idx]
            cluster2 = clusters[c2_idx]
            for i, j in product(cluster1, cluster2):
                sim = DataStructs.TanimotoSimilarity(fingerprints[i], fingerprints[j])
                inter_sims.append(sim)

    avg_intra = sum(intra_sims) / len(intra_sims) if intra_sims else 0.0
    avg_inter = sum(inter_sims) / len(inter_sims) if inter_sims else 0.0

    return avg_intra, avg_inter
