"""
A script to reliably access and download the relevant regulatory feature
files from biomart.

Example:
python ensembl_downloads.py
"""
import pybiomart
import os
import pickle


def download_regulatory_features(bm_data, feat, attr_list, epi_list, chr_list):
    # Accessing, downloading & saving files
    for epigenome in epi_list:
        filename = (feat.replace(" ", "_") + "_" +
                    epigenome.replace(" ", "_") + ".pkl")
        filepath = "data/biomart_data/" + filename
        if os.path.isfile(filepath):
            print(filename + " already downloaded.")
            continue
        print("Downloading: " + filepath)
        query_filters = {"regulatory_feature_type_name": feat,
                         "epigenome_name": epigenome,
                         "chromosome_name": chr_list}
        feature_sites = bm_data.query(attributes=attr_list,
                                      filters=query_filters,
                                      use_attr_names=True)
        feature_sites.to_pickle(filepath)


def main():
    # Connecting to the relevant biomart server
    bm_server = pybiomart.Server(host='http://www.ensembl.org')
    mart = bm_server["ENSEMBL_MART_FUNCGEN"]
    bm_data = mart["mmusculus_regulatory_feature"]

    # Defining necessary query lists
    chromosome_list = [str(x) for x in list(range(20)[1:])]
    feature_list = ["CTCF Binding Site", "Enhancer", "Promoter",
                    "TF Binding Site"]
    attribute_list = ["chromosome_name", "chromosome_start", "chromosome_end"]
    with open('data/epigenome_list.data', 'rb') as filehandle:
        epigenome_list = pickle.load(filehandle)

    # Downloading necessary features:
    for feature in feature_list:
        download_regulatory_features(bm_data, feature, attribute_list,
                                     epigenome_list, chromosome_list)


if __name__ == "__main__":
    main()
