"""
A script to integrate the validated IAP labels with the entire IAP dataset.

Example:
python assemble_dataset.py
"""
import pybedtools
import pandas as pd


def check_intersection(x, val_df):
    if x["name"] in val_df["index"]:
        x["val_result"] = val_df.loc[x["index"], "val_result"]
        return x
    else:
        x["val_result"] = "untested"
        return x


def main():
    # Importing validated IAPs and total IAP dataset:
    total_iap = pybedtools.BedTool(
        "data/clean_beds/mm10.IAP.mended.extent.bed")
    validate_iap = pybedtools.BedTool(
        "data/feature_beds/IAP_validation.July2019.stranded.bed")

    # Generating dataframes:
    total_df = total_iap.to_dataframe()
    val_df = validate_iap.to_dataframe()

    # Adding element ids to validated IAP dataset:
    column_names = ["chrom1", "start1", "end1", "strand1", "name", "un1",
                    "un2", "val_result", "chrom2", "start3", "end3",
                    "element_id", "length1", "strand2", "length2"]
    val_id_df = validate_iap.intersect(total_iap, wo=True).to_dataframe(
                   names=column_names)
    val_id_df = val_id_df.drop(columns=["chrom1", "start1", "end1", "strand1",
                                        "name", "un1", "un2", "chrom2", "start3",
                                        "end3", "length1", "strand2", "length2"])
    val_id_df.loc[:, "element_id"] = val_id_df.loc[:, "element_id"].apply(pd.to_numeric)

if __name__ == "__main__":
    main()
