"""
A script that cleans .bed files, removing any non-consensus intervals.

Example:
python clean_bedfile.py mm10.IAP.mended.extent
"""

import pybedtools
import argparse

def _parse_args():
    parser = argparse.ArgumentParser(description=
        'Clean .bed files for the mouse genome.')
    parser.add_argument('filename', metavar='N', type=str,
                    help='The filename of the .bed to be cleaned.')
    return parser.parse_args()

def main():
    args = _parse_args()
    filename = args.filename
    dirty_iap = pybedtools.BedTool('data/ME_beds/' + filename + '.bed')
    dirty_iap = dirty_iap.sort(chrThenSizeA=True)
    chr_list = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8',
                'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15',
                'chr16', 'chr17', 'chr18', 'chr19', 'chrX', 'chrY']
    iap = dirty_iap.filter(lambda x: any(x.chrom == chr for chr in chr_list))
    iap.saveas('data/clean.' + filename + '.bed')

if __name__ == '__main__':
    main()
