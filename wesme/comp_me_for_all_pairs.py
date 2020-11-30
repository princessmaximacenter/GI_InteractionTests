"""
compute WeSME pvalue for all pairs
each mutation type separately

>> python comp_me_for_all_pairs.py BRCA smut 0.1

"""
import argparse
import os
import logging
import pandas
logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger()

# added data_dir and results_dir by JD
from config import data_dir, wrs_dir, results_dir
from mut_ex import ws_mut_ex, misc_mut_ex
from data_proc import io_utils

# arguments
parser = argparse.ArgumentParser()
parser.add_argument("can", help="cancer type", type=str)
parser.add_argument("mtype", help="mutation type (smut, amp, del, cnv or alt)", type=str)
parser.add_argument("pth", help="pvalue threshold to write in the output file", type=float)
parser.add_argument("-m", "--mut_input", type=str, help="mutation file name")
parser.add_argument("-s", "--sampling_dir", type=str, help="sampling directory")
parser.add_argument("-o", "--output", type=str, help="output file name")

args = parser.parse_args()
cantype, mtype, pth = args.can, args.mtype, args.pth
logging.debug("comp ME for %s, %s" %(cantype, mtype))

# INPUT_FILES/DIRECTORIES
# JD: changed "data/" in data_dir
if args.mut_input is not None:
    mut_file = args.mut_input
else:
    mut_dir = data_dir + cantype + "/"
    mut_file = mut_dir + "_".join([cantype, mtype]) + "_list.txt"

if args.sampling_dir is not None:
    cur_wrs_dir = args.sampling_dir
else:
    cur_wrs_dir = wrs_dir + cantype + "/" + mtype + "/"

# OUTPUT_FILES
if args.output is not None:
    pv_file = args.output
else:
    pv_file = results_dir + cantype + "_" + mtype + "_me_pvs_"+str(pth)+".txt"

# JD: create results_dir if not exists
try:
    os.makedirs(os.path.dirname(pv_file))
except OSError:
    if not os.path.isdir(os.path.dirname(pv_file)):
        raise

# read mutation profile
alt_list_dic, samples = io_utils.read_mut_list(mut_file, samples=True)
gene_ks = dict([(gene, len(alt_list_dic[gene])) for gene in alt_list_dic])  # gene -> cover size dic
nsamples = len(samples)

# read weighted sampling files for all k
# use the weighted sampling files in cur_wrs_dir
kfiles = os.listdir(cur_wrs_dir)
# store weighted sampling for each k: k --> list of arrays
ws_k_cover_dic = {}
for kf in kfiles:
    # read random file for k
    k = int(kf.split(".")[0].split("_")[1])
    ws_k_cover_dic[k], samples = io_utils.read_mut_list(cur_wrs_dir + kf, genes=False)

# measure ME/CO
h_pv_dic = {}  #added by JD: hypergeometric
ws_pv_dic = {}  # WeSME
jaccard_dic = {}  # jaccard index

# compute ME for all pairs
genes = alt_list_dic.keys()
ws_ex_cover_sizes_dic = {}
for i in range(len(genes)):
    gene1 = genes[i]
    if i%1000 == 0:
        logging.info("%d " % i)
    for j in range(i+1, len(genes)):
        gene2 = genes[j]
        gene_pairs = (gene1, gene2)
        # compute pvalue using weighted sampling
        covers = [alt_list_dic[gene1], alt_list_dic[gene2]]
        ws_pv, ws_ex_cover_sizes_dic = ws_mut_ex.compute_me_co_pv_ws\
            (covers, ws_k_cover_dic, max_pair_num=10**5, me_co="me", ws_ex_cover_sizes_dic=ws_ex_cover_sizes_dic)
        if ws_pv > pth:  # if pvalue is not significant enough
            continue
        ws_pv_dic[gene_pairs] = ws_pv # WeSME pvalue
        # compute pvalue using hypergeometric test
        h_pv_dic[gene_pairs] = misc_mut_ex.compute_me_pv_hypergeom(alt_list_dic[gene1], alt_list_dic[gene2], nsamples)
        # compute jaccard index
        jaccard_dic[gene_pairs] = misc_mut_ex.compute_jaccard(alt_list_dic[gene1], alt_list_dic[gene2])


# write the results
rows = []
for gp in ws_pv_dic:
    rows.append((gp[0], gp[1], jaccard_dic[gp], h_pv_dic[gp], ws_pv_dic[gp]))

# write the results
all_pvs = pandas.DataFrame(data=rows, columns=["gene1", "gene2", "jaccard index", "pv (fisher)", "pv (ws)"])
all_pvs.to_csv(pv_file, sep="\t", index=False)
