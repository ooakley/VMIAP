"""
A script to cluster genomic intervals based on their genomic context.

Example:
python selection_clustering.py
"""
import pandas as pd
import umap
from sklearn.feature_selection import VarianceThreshold
from sklearn.cluster import KMeans


def main():
    # Importing datasets:
    total_df = pd.read_pickle("data/iap_variances.pkl")
    total_df["val_result"] = total_df["val_result"].replace("-1", "Untested")

    # Selecting high variance features:
    selector = VarianceThreshold(threshold=0.2)
    selector.fit(total_df.iloc[:, 8:])
    support = selector.get_support(indices=True)
    kbest_total = total_df.iloc[:, 8:][total_df.iloc[:, 8:].columns[support]]

    # UMAP fitting & embedding:
    reducer = umap.UMAP(n_neighbors=15, n_components=2, verbose=False, random_state=1)
    reducer.fit(kbest_total)
    embedding_total = reducer.transform(kbest_total)

    # k-means clustering on UMAP embeddings:
    kmeans = KMeans(n_clusters=5, random_state=0).fit(embedding_total)
    total_category_labels = pd.Series(kmeans.labels_)
    total_category_labels.index = total_df["element_id"].astype(int).to_list()

    # Assigning cluster values to interval dataframe:
    sorted_total_df = total_df.iloc[:, 0:7].copy(deep=True)
    sorted_total_df["cluster_assignments"] = total_category_labels

    label_order = ["A", "B", "C", "D", "E"]

    for i in range(5):
        sorted_total_df.loc[:, "cluster_assignments"] = \
            sorted_total_df.loc[:, "cluster_assignments"].replace(i, label_order[i])

    # Saving dataframe:
    sorted_total_df.to_pickle("data/iap_clustered.pkl")


if __name__ == "__main__":
    main()
