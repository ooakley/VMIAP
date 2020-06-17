"""
A script to evaluate feature variance from in the surrounding genomic context of a set of genomic
intervals.

Example:
python generate_feature_dev_variance.py --shuffle False
"""
import argparse
import os
import pandas as pd
import pybedtools
import pybedtools.featurefuncs as featurefuncs
import numpy as np

EPIGENOME_GROUPS = ["heart"]
FEATURES = ["CTCF_Binding_Site", "Enhancer", "Promoter"]


def _parse_args():
    parser = argparse.ArgumentParser(
        description='Evaluate variances of counts around defined genomic intervals')
    parser.add_argument(
        '--shuffle', type=bool, default=False,
        help='Whether or not to shuffle the intervals randomly.')
    parser.add_argument(
        '--window', type=int, default=10000,
        help='The length of window around genomic intervals to be queried.')
    return parser.parse_args()


def generate_counts(epigenome_bed, intervals_bed, window_size):
    # Defining names & centre of each IAP:
    names = ["chrom", "start", "end", "element_id", "length", "strand", "hits"]
    interval_center = intervals_bed.each(
        featurefuncs.center, width=10).saveas()

    # Finding hits:
    center_slop = interval_center.slop(b=window_size, genome="mm10")
    counts = center_slop.intersect(epigenome_bed, c=True)
    hits_df = counts.to_dataframe(names=names)
    hits_df.index = hits_df["element_id"].astype(int).to_list()

    return hits_df["hits"]


def main():
    # Parsing command line arguments:
    args = _parse_args()

    # Loading datasets:
    dataset_list = os.listdir("data/biomart_data/")
    total_df = pd.read_pickle("data/labelled_iaps.pkl")
    interval_bed = pybedtools.BedTool(
        "data/clean_beds/mm10.IAP.mended.extent.bed")

    # Shuffling dataset if necessary:
    if args.shuffle:
        interval_bed = interval_bed.shuffle(genome='mm10', chrom=True, seed=0).saveas()
        total_df = interval_bed.to_dataframe()
        total_df = total_df.rename(
            columns={"name": "element_id", "score": "length"})
        total_df.index = total_df["element_id"].astype(int).to_list()

    # Iterating through groups of features & epigenomes:
    i = 0  # progress indicator
    max_epi_no = 20  # maximum number of epigenomes queried
    for genomic_feature in FEATURES:
        # Filtering entire biomart download set for specific feature:
        gen_str = len(genomic_feature)
        feature_filtered = [
            x for x in dataset_list if x[0:gen_str] == genomic_feature]

        for epigenome_group in EPIGENOME_GROUPS:
            # Generating progress indicator:
            if i % 10 == 0:
                print(str(i) + " epigenomes processed.")
            i += 1

            # Filtering filename list:
            tot_str = gen_str + len(epigenome_group)
            epigenome_feature_filtered = [
                x for x in feature_filtered if x[gen_str+1:tot_str+1] == epigenome_group]

            # Initialising count array:
            count_array = np.zeros((total_df.shape[0], 1))

            # Generating & allocating counts:
            j = 0
            for epigenome_name in epigenome_feature_filtered:
                j += 1
                if j > max_epi_no:
                    print("skip")
                    continue
                dataset = pd.read_pickle("data/biomart_data/" + epigenome_name)
                dataset["chromosome_name"] = [
                    "chr" + str(x) for x in list(dataset["chromosome_name"])]
                dataset_bed = pybedtools.BedTool.from_dataframe(dataset)

                counts_instance = generate_counts(
                    dataset_bed, interval_bed, args.window).to_numpy()
                counts_instance = np.expand_dims(counts_instance, 1)
                count_array = np.concatenate((count_array, counts_instance), axis=1)

            # Computing variance over count array:
            var_array = np.var(count_array[:, 1:], axis=1)
            mean_array = np.mean(count_array[:, 1:], axis=1)

            # Formatting:
            var_series = pd.Series(var_array)
            var_series.index = total_df["element_id"].astype(int).to_list()
            mean_series = pd.Series(mean_array)
            mean_series.index = total_df["element_id"].astype(int).to_list()

            # Assignment:
            column_id = genomic_feature + "_" + epigenome_group
            total_df[column_id + "_variance"] = var_series
            print(column_id + "_variance:")
            print(var_series)
            total_df[column_id + "_mean"] = mean_series

    # Saving dataframe:
    if args.shuffle:
        total_df.to_pickle("data/rand_feature_dev_variance.pkl")
    else:
        total_df.to_pickle("data/iap_feature_dev_variance.pkl")


if __name__ == "__main__":
    main()
