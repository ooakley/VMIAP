"""
A script to evaluate feature distances from a set of genomic intervals.

Example:
python generate_feature_distances.py
"""
import os
import pandas as pd
import pybedtools
import pybedtools.featurefuncs as featurefuncs


def generate_distances(epigenome_bed, intervals_bed):
    # Defining names & centre of each IAP
    names = ["chrom", "start", "end", "element_id", "length", "strand", "1",
             "2", "3", "distance"]
    interval_center = intervals_bed.each(
        featurefuncs.center, width=10).saveas().sort()

    # Finding hits
    distances = interval_center.closest(epigenome_bed.sort(),
                                        d=True, t="first")
    distances_df = distances.to_dataframe(names=names)
    distances_df.index = distances_df["element_id"].astype(int).to_list()
    return distances_df["distance"]


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

        # Excluding enhancer data:
        if filename[:8] == "Enhancer":
            print("skipped")
            continue

        # Loading epigenome data:
        dataset = pd.read_pickle("data/biomart_data/" + filename)
        dataset["chromosome_name"] = [
            "chr" + str(x) for x in list(dataset["chromosome_name"])]
        dataset_bed = pybedtools.BedTool.from_dataframe(dataset)

        # Generating counts:
        feature_id = filename[:-4]
        total_df[feature_id + "_distance1"] = generate_distances(
            dataset_bed, interval_bed)
        # total_df[feature_id + "_distance2"] = generate_distances(
        #    dataset_bed, interval_bed, 2)

    total_df.to_pickle("data/iap_distances.pkl")


if __name__ == "__main__":
    main()
