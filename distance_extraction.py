"""
A script to extract the distance to the nearest genomic element of interest
(such as a regulatory element) from a set of regions (such as a bed file of
IAPs).

Example:
python distance_extraction.py mm10.IAP.mended.extent enhancers
"""
import pybedtools
import argparse


def _parse_args():
    parser = argparse.ArgumentParser(description='Extract distances to genomic\
                                     features for given set of intervals.')
    parser.add_argument('intervals', metavar='N', type=str,
                        help='The filename of the intervals .bed to be\
                        queried.')
    parser.add_argument('features', metavar='N', type=str,
                        help='The filenames of the features .bed to be\
                        queried.')
    return parser.parse_args()


def find_distances(interval_bed, feature_bed, feature_id):
    # Finding closest tss sequences & calculating their distances:
    dist_bed = interval_bed.sort().closest(
        feature_bed.sort(), D="b", t="first")

    # Generating list of column names for dataframe:
    interval_names = list(interval_bed.to_dataframe().columns)
    feature_names = list(feature_bed.to_dataframe().columns)
    feature_names = [e + "1" for e in feature_names]
    col_names = []
    col_names.extend(interval_names)
    col_names.extend(feature_names)
    col_names.append(feature_id)

    # Generating clean dataframe:
    dist_df = dist_bed.to_dataframe(names=col_names)
    dist_df = dist_df.drop(columns=feature_names)
    return dist_df


def main():
    # Parsing arguments:
    args = _parse_args()
    intervals = args.intervals
    features = args.features

    # Loading beds & finding distances:
    path = "data/clean_beds/"
    interval_bed = pybedtools.BedTool(path + intervals + ".bed")
    feature_bed = pybedtools.BedTool(path + features + ".bed")
    dist_df = find_distances(interval_bed, feature_bed)

    # Saving results:
    dist_df.to_pickle("data/extracted_dataframes/"+features+"_distance.pkl")


if __name__ == "__main__":
    main()
