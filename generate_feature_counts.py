"""
A script to evaluate feature counts in and around a set of genomic intervals.

Example:
python generate_feature_counts.py
"""
import os
import pandas as pd
import pybedtools
import pybedtools.featurefuncs as featurefuncs


def generate_counts(epigenome_bed, intervals_bed, window_size):
    # Defining names & centre of each IAP
    names = ["chrom", "start", "end", "element_id", "length", "strand", "hits"]
    interval_center = intervals_bed.each(
        featurefuncs.center, width=10).saveas()

    # Finding hits
    center_slop = interval_center.slop(b=window_size, genome="mm10")
    counts = center_slop.intersect(epigenome_bed, c=True)
    hits_df = counts.to_dataframe(names=names)
    hits_df.index = hits_df["element_id"].astype(int).to_list()

    return hits_df["hits"]


def main():
    # Loading datasets:
    dataset_list = os.listdir("data/biomart_data/")
    total_df = pd.read_pickle("data/labelled_iaps.pkl")
    interval_bed = pybedtools.BedTool(
        "data/clean_beds/mm10.IAP.mended.extent.bed")

    i = 0
    for filename in dataset_list:
        # Generating progress indicator:
        if i % 20 == 0:
            print(str(i) + " datasets processed.")
        i += 1

        # Loading epigenome data:
        dataset = pd.read_pickle("data/biomart_data/" + filename)
        dataset["chromosome_name"] = [
            "chr" + str(x) for x in list(dataset["chromosome_name"])]
        dataset_bed = pybedtools.BedTool.from_dataframe(dataset)

        # Generating counts:
        feature_id = filename[:-4]
        total_df[feature_id + "_counts_100"] = generate_counts(
            dataset_bed, interval_bed, 100)
        total_df[feature_id + "_counts_1k"] = generate_counts(
            dataset_bed, interval_bed, 1000)
        total_df[feature_id + "_counts_5k"] = generate_counts(
            dataset_bed, interval_bed, 5000)
        total_df[feature_id + "_counts_10k"] = generate_counts(
            dataset_bed, interval_bed, 10000)
        total_df[feature_id + "_counts_50k"] = generate_counts(
            dataset_bed, interval_bed, 50000)

    total_df.to_pickle("data/iap_counts.pkl")


if __name__ == "__main__":
    main()
