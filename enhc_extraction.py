"""
A script to extract distances of genomic elements from the nearest regulatory
elements.

Example:
python reg_extraction.py
"""
import pybedtools


def enhancer_extraction(gen_bed, enhc_bed):
    # Extracting distances from enhancers:
    enhc_dist_bed = gen_bed.sort().closest(enhc_bed.sort(), D="b", t="first")

    # Generating list of column names for dataframe:
    gen_names = list(gen_bed.to_dataframe().columns)
    enhc_names = list(enhc_bed.to_dataframe().columns)
    enhc_names = [e + "1" for e in enhc_names]
    col_names = []
    col_names.extend(gen_names)
    col_names.extend(enhc_names)
    col_names.append("enhc_dist")

    # Generating clean dataframe:
    enhc_dist_df = enhc_dist_bed.to_dataframe(names=col_names)
    enhc_dist_df = enhc_dist_df.drop(columns=enhc_names)
    return enhc_dist_df


def main():
    iap_bed = pybedtools.BedTool('data/clean_beds/mm10.IAP.mended.extent.bed')
    enhc_bed = pybedtools.BedTool('data/clean_beds/enhancers.bed')
    enhc_dist_df = enhancer_extraction(iap_bed, enhc_bed)
    enhc_dist_df.to_pickle("data/extracted_dataframes/enhc_dist.pkl")


if __name__ == "__main__":
    main()
