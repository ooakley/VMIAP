"""
A script to integrate the validated IAP labels with the entire IAP dataset.

Example:
python assemble_dataset.py
"""
import pybedtools


def main():
    # Importing validated IAPs and total IAP dataset:
    total_iap = pybedtools.BedTool(
        "data/clean_beds/mm10.IAP.mended.extent.bed")
    validate_iap = pybedtools.BedTool(
        "data/feature_beds/IAP_validation.July2019.stranded.bed")

    # Generating dataframes:
    total_df = total_iap.to_dataframe()
    val_df = validate_iap.to_dataframe()

    # Tidying dataframes:
    val_df = val_df.rename(columns={"name": "strand", "strand": "score",
                                    "score": "gene", "itemRgb": "val_result"})
    total_df = total_df.rename(columns={"name": "element_id",
                                        "score": "length"})

    # Generating column headers:
    total_names = list(total_df.columns)
    val_names = list(val_df.columns)
    val_names = [e + "1" for e in val_names]
    col_names = []
    col_names.extend(val_names)
    col_names.extend(total_names)
    col_names.extend(["length1"])

    # Assigning element ids to validated IAPs:
    val_id_df = validate_iap.intersect(
        total_iap, wo=True).to_dataframe(names=col_names)
    val_id_df = val_id_df.drop(columns=val_names[:-1])

    # Making the indices identical to element IDs:
    total_df.index = total_df["element_id"].astype(int).to_list()
    val_id_df.index = val_id_df["element_id"].astype(int).to_list()

    # Assigning validation status to all autosomal IAPs:
    total_df["val_result"] = ""
    for element_id, row in total_df.iterrows():
        if element_id in val_id_df["element_id"]:
            total_df.loc[element_id, "val_result"] = val_id_df.loc[
                element_id, "val_result1"]
        else:
            total_df.loc[element_id, "val_result"] = "untested"

    # Printing results & saving to file:
    print(total_df["val_result"].value_counts())
    total_df.to_pickle("data/labelled_iaps.pkl")


if __name__ == "__main__":
    main()
