import numpy as np
import pandas as pd
import pickle
import umap
import matplotlib.pyplot as plt
from pathlib import Path

import numpy as np
import pandas as pd
import pickle
import umap
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from sklearn.cluster import KMeans, AgglomerativeClustering
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import linkage, dendrogram


def load_reference_cache(species="human", cache_dir="numba_cache"):
    """Load matrices + metadata from cache"""
    cache_file = Path(cache_dir) / f"ref_matrices_{species}.npy"
    metadata_file = Path(cache_dir) / f"ref_metadata_{species}.pkl"
    
    if not cache_file.exists() or not metadata_file.exists():
        raise FileNotFoundError(f"No cache found for {species} in {cache_dir}")
    
    matrices = np.load(cache_file, mmap_mode="r")
    with open(metadata_file, "rb") as f:
        cache_data = pickle.load(f)
        metadata = cache_data["metadata"]
    
    return matrices, metadata


def reduce_dimensions(X, method="umap", n_components=2, **kwargs):
    """Dimensionality reduction helper"""
    if method == "umap":
        reducer = umap.UMAP(n_components=n_components, random_state=42, **kwargs)
    elif method == "tsne":
        reducer = TSNE(n_components=n_components, random_state=42, **kwargs)
    elif method == "pca":
        reducer = PCA(n_components=n_components, random_state=42)
    else:
        raise ValueError(f"Unknown method {method}")
    return reducer.fit_transform(X)


def plot_embedding(embedding, metadata, color_mode="family", highlight_hlas=None, title="Embedding", output_path=None):
    """Plot 2D embedding with flexible coloring"""
    if color_mode == "allele":
        labels = metadata["hla"]
    elif color_mode == "family":
        # Robust extract: A, B, C from "HLA-A*02:01"
        labels = metadata["hla"].str.extract(r'([ABC])', expand=False).fillna("Other")
    elif color_mode == "highlight":
        if highlight_hlas is None:
            raise ValueError("Must provide highlight_hlas when color_mode='highlight'")
        labels = metadata["hla"].apply(lambda h: h if h in highlight_hlas else "Other")
    else:
        raise ValueError(f"Unknown color_mode: {color_mode}")
    
    codes, uniques = pd.factorize(labels)

    plt.figure(figsize=(10, 8))
    scatter = plt.scatter(
        embedding[:, 0],
        embedding[:, 1],
        c=codes,
        cmap="Spectral",
        s=40,
        alpha=0.8
    )
    plt.title(title, fontsize=14)
    plt.xlabel("Dim 1")
    plt.ylabel("Dim 2")

    handles, _ = scatter.legend_elements()
    plt.legend(handles, uniques, title="HLA", bbox_to_anchor=(1.05, 1), loc="upper left")
    plt.tight_layout()
    if output_path:
        plt.savefig(output_path, dpi=300)
        print(f"Saved {output_path}")
    plt.close()


def plot_similarity_heatmap(X, metadata, metric="cosine", output_path="similarity_heatmap.png"):
    """Plot pairwise similarity heatmap of motifs"""
    dists = squareform(pdist(X, metric=metric))
    sim = 1 - dists  # similarity from distance
    
    plt.figure(figsize=(12, 10))
    sns.heatmap(sim, cmap="viridis", xticklabels=False, yticklabels=False)
    plt.title(f"HLA motif similarity heatmap ({metric})")
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    print(f"Saved {output_path}")
    plt.close()


def plot_dendrogram(X, metadata, method="average", metric="cosine", output_path="dendrogram.png"):
    """Hierarchical clustering dendrogram"""
    Z = linkage(X, method=method, metric=metric)
    plt.figure(figsize=(14, 6))
    dendrogram(Z, labels=metadata["hla"].values, leaf_rotation=90, leaf_font_size=6)
    plt.title(f"Hierarchical clustering ({method}, {metric})")
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    print(f"Saved {output_path}")
    plt.close()


def run_analysis(species="human", cache_dir="numba_cache"):
    """Run all analyses: UMAP, t-SNE, PCA, heatmap, dendrogram"""
    matrices, metadata = load_reference_cache(species, cache_dir)
    X = matrices.reshape(matrices.shape[0], -1)

    # Embeddings
    for method in ["umap", "tsne", "pca"]:
        embedding = reduce_dimensions(X, method=method, n_components=2)
        plot_embedding(
            embedding, metadata,
            color_mode="family",
            title=f"{method.upper()} projection of HLA motifs",
            output_path=f"{method}_projection.png"
        )

    # Similarity heatmap
    plot_similarity_heatmap(X, metadata, metric="cosine")

    # Dendrogram
    plot_dendrogram(X, metadata, method="average", metric="cosine")

def inspect_npk(npy_path="numba_cache/ref_matrices_human.npy",
                pkl_path="numba_cache/ref_metadata_human.pkl",
                show_example=True, max_hla=10):
    """
    Inspect contents of the reference .npy and .pkl files with a professional summary.

    Parameters
    ----------
    npy_path : str
        Path to the .npy file containing reference matrices.
    pkl_path : str
        Path to the .pkl file containing metadata and configuration.
    show_example : bool
        If True, prints the first matrix as an example.
    max_hla : int
        Maximum number of unique HLA alleles to display in summary.
    """

    # Load .npy
    matrices = np.load(npy_path, mmap_mode="r")
    n_matrices, max_positions, n_aas = matrices.shape
    mem_mb = matrices.nbytes / (1024 ** 2)

    print("=" * 60)
    print("**NPY SUMMARY**")
    print("-" * 60)
    print(f"Total matrices      : {n_matrices:,}")
    print(f"Peptide length (max): {max_positions}")
    print(f"Amino acids encoded : {n_aas}")
    print(f"Memory footprint    : {mem_mb:.2f} MB")
    print("=" * 60)

    if show_example:
        print("[NPY] Example matrix [0]:")
        print(matrices[0])

    # Load .pkl
    with open(pkl_path, "rb") as f:
        cache_data = pickle.load(f)

    metadata = cache_data["metadata"]

    print("\n" + "=" * 60)
    print("**PKL SUMMARY**")
    print("-" * 60)
    print(f"Species             : {cache_data['species']}")
    print(f"Max positions       : {cache_data['max_positions']}")
    print(f"Amino acids         : {', '.join(cache_data['amino_acids'])}")

    unique_hlas = metadata["hla"].unique()
    print(f"Unique HLA alleles  : {len(unique_hlas):,}")
    print(f"Example HLAs        : {', '.join(unique_hlas[:max_hla])}{' ...' if len(unique_hlas) > max_hla else ''}")

    print("\nMetadata head:")
    print(metadata.head())

    print("=" * 60)

    return matrices, cache_data


# if __name__ == "__main__":
#     inspect_npk()
#     run_analysis()
