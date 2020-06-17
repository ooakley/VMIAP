"""
A script to cluster genomic intervals based on their genomic context.

Example:
python selection_clustering.py
"""
import pandas as pd
from sklearn.feature_selection import VarianceThreshold
from sklearn.preprocessing import MaxAbsScaler
from sklearn.cluster import KMeans


def main():
    # Importing the relevant dataframes:
    activity_var_df = pd.read_pickle("data/activity_variances.pkl")
    ctcf_var_df = pd.read_pickle("data/iap_variances.pkl")
    tissue_var_df = pd.read_pickle("data/activity_tissue_variances.pkl")
    concat_list = [activity_var_df.sort_index(),
                   ctcf_var_df.iloc[:, 7:].sort_index(),
                   tissue_var_df.iloc[:, 7:].sort_index()]
    total_df = pd.concat(concat_list, axis=1)

    # Selecting high variance features:
    selector = VarianceThreshold(threshold=0.005)
    selector.fit(total_df.iloc[:, 7:])
    support = selector.get_support(indices=True)
    kbest_total = total_df.iloc[:, 7:][total_df.iloc[:, 7:].columns[support]]

    # Normalising variances & counts:
    transformer = MaxAbsScaler().fit(kbest_total)
    kbest_total_abs = transformer.transform(kbest_total)

    # k-means clustering on filtered dataset:
    kmeans = KMeans(n_clusters=5, random_state=0).fit(kbest_total_abs)
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
