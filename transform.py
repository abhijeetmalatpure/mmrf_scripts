import os, pandas as pd
import requests
from json import encoder, decoder
from ratelimit import limits
from subprocess import Popen, PIPE

ONE_SEC=1

# Usage example append_seg("/N/absolute_file_path/file_name.seg", "dir.CNA_Exome")
def append_seg(f):
    tsv = pd.read_csv(f, delimiter='\t', header=0)
    if os.path.isfile("/N/slate/abhmalat/MMRF_CoMMpass_IA16a/segment_combined.seg"):
        tsv.to_csv("/N/slate/abhmalat/MMRF_CoMMpass_IA16a/segment_combined.seg", mode='a+', header=False, index=False, sep='\t')
    else:
        tsv.to_csv("/N/slate/abhmalat/MMRF_CoMMpass_IA16a/segment_combined.seg", mode='a+', header=True, index=False, sep='\t')
    return


@limits(calls=5, period=ONE_SEC)
def call_api(ids):
    hgnc_url = "http://rest.genenames.org/fetch/ensembl_gene_id/"
    response = []
    for id in ids:
        response.append(str(Popen(['curl','-H','Accept:application/json',
                          hgnc_url + id], stdout=PIPE, stderr=PIPE,
                                     shell=False).communicate()[0]))
    return response


def get_hugo_symbol(f):
    tsv = pd.read_csv(open(f), delimiter='\t', header=0, nrows=5)
    tsv["Hugo_Symbol"] = call_api(tsv["Gene"])


def main():

    mmrf_root = "/N/slate/abhmalat/MMRF_CoMMpass_IA16a"
    fl = os.path.join(mmrf_root, 'copy_number_estimates', 'MMRF_CoMMpass_IA16a_CNA_Exome_PerGene_LargestSegment.txt')

    get_hugo_symbol(fl)

    # for root, dirs, files in os.walk(mmrf_root, topdown=True):
    #     for file in files:
    #         if not root == mmrf_root:
    #             # Call append_seg to combine all .seg files into segment_combined.seg
    #             if file.endswith(".seg"):
    #                 append_seg(file)



if __name__ == '__main__':
    main()