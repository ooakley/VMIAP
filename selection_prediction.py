"""
A script to select features from the dataset useful in predicting variable methylation status
and to output a dataset containing these predictions.

Example:
python selection_prediction.py
"""
import seaborn as sns
import pybedtools
import pandas as pd
import numpy as np
import umap
from sklearn.feature_selection import VarianceThreshold, SelectKBest, chi2, mutual_info_classif
from sklearn.semi_supervised import LabelSpreading


def mi_kbest_selector(data, labels, k=30):
    selector = SelectKBest(mutual_info_classif, k)
    selector.fit(data, labels)
    return data[data.columns[
        selector.get_support(indices=True)]], selector.get_support(indices=True)


def main():
    "Import data, perform feature selection & imputation of methylation state."

    # Importing data:
    total_df = pd.read_pickle("data/iap_counts.pkl")
    total_df["val_result"] = total_df["val_result"].replace("-1", "Untested")

    # Generating integer encodings of target values:
    total_df["integer_encodings"] = total_df["val_result"]
    total_df.loc[:, "integer_encodings"] = \
        total_df.loc[:, "integer_encodings"].replace("Untested", -1)
    total_df.loc[:, "integer_encodings"] = \
        total_df.loc[:, "integer_encodings"].replace("True ME", 1)
    total_df.loc[:, "integer_encodings"] = \
        total_df.loc[:, "integer_encodings"].replace("False-positive", 2)
    total_df.loc[:, "integer_encodings"] = \
        total_df.loc[:, "integer_encodings"].replace("Tissue-specific", 3)
    cols = []
    cols = total_df.columns.tolist()
    cols = cols[:7] + [cols[-1]] + cols[7:-1]
    total_df = total_df[cols]

    # Extracting dataframe of all validated IAPs:
    val_df = total_df[total_df.val_result != "Untested"]

    # Selecting most relevant features based on mutual information:
    k_best_df, support = mi_kbest_selector(val_df.iloc[:, 9:], val_df["integer_encodings"])
    kbest_total = total_df.iloc[:, 9:][total_df.iloc[:, 9:].columns[support]]

    # Imputing labels for the entire dataset:
    print(total_df.loc[:, "val_result"].value_counts())
    label_prop_model = LabelSpreading()
    label_prop_model.fit(kbest_total, total_df.loc[:, "integer_encodings"])
