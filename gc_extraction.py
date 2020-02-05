"""
A script to extract nearest TSS from the flanking sequences of given IAPs

Example:
python gc_extraction.py
"""
import pybedtools
import statistics as stat
import math
from Bio import SeqIO
from Bio.SeqUtils import GC


def gcflank_extraction(input_bed, fa_file):
    gcflank_seq = input_bed.flank(genome='mm10', b=2000).sequence(fi=fa_file)
    gcflank_ext = [GC(rec.seq) for rec in SeqIO.parse(gcflank_seq.seqfn,
                                                      "fasta")]
    gcflank_list = []
    for i in range(math.floor(len(gcflank_ext)/2)):
        j = i*2
        gcflank_list.append(stat.mean(gcflank_ext[j:j+1]))

    print(len(gcflank_list))
    return gcflank_list


def main():
    mm_fasta = pybedtools.BedTool('data/GRCm38.p6.genome.fa')
    iap = pybedtools.BedTool('data/clean_beds/mm10.IAP.mended.extent.bed')
    gcflank_iap = gcflank_extraction(iap, mm_fasta)
    iap_df = iap.to_dataframe()
    iap_df["flank_gc"] = gcflank_iap
    iap_df.to_pickle("data/extracted_dataframes/flank_gc.pkl")


if __name__ == "__main__":
    main()
