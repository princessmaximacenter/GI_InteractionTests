"""
compute WeSME pvalue for all pairs
each mutation type separately

>> python comp_me_for_all_pairs_PAN.py BRCA smut 0.1

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

# JD
from operator import add

# arguments
parser = argparse.ArgumentParser()
parser.add_argument("nws", help="number of weighted samplings", type=str)
parser.add_argument("mtype", help="mutation type (smut, amp, del, cnv or alt)",
                    type=str)
parser.add_argument("pth", help="pvalue threshold to write in the output file",
                    type=float)
parser.add_argument("proj", help="project name (dkfz/stjude", type=str)
parser.add_argument("come", help="interaction type: co or me",
                    type=str)
parser.add_argument("-m", "--mut_input", type=str, help="mutation file name")
parser.add_argument("-s", "--sampling_dir", type=str, help="sampling directory")
parser.add_argument("-o", "--output", type=str, help="output file name")

args = parser.parse_args()
nws, mtype, pth, proj, come = args.nws, args.mtype, args.pth, args.proj, args.come
logging.debug("comp %s for PAN, %s" %(come, mtype))

# INPUT_FILES/DIRECTORIES
# JD: changed "data/" in data_dir
# JD: changed input directories into ones without cancer types
mut_file = None
if args.mut_input is not None:
    mut_file = args.mut_input
else:
    mut_dir = data_dir

if args.sampling_dir is not None:
    cur_wrs_dir = args.sampling_dir
else:
    cur_wrs_dir = wrs_dir

# OUTPUT_FILES
if args.output is not None:
    pv_file = args.output
else:
    pv_file = results_dir + "PAN" + "_" + mtype + "_" + come + "_pvs_" + str(pth) + ".txt"

# JD: create results_dir if not exists
try:
    os.makedirs(os.path.dirname(pv_file))
except OSError:
    if not os.path.isdir(os.path.dirname(pv_file)):
        raise

# JD: read cancertypes from file
cantypes = [line.rstrip('\n') for line in open(proj + '/cantypes.txt')]

alt_list_dic = {}
alt_samples = {}
samples = {}
gene_ks = {}
nsamples = {}
kfiles = {}
ws_k_cover_dic = {}
ws_ex_cover_sizes_dic = {}
ex_cover_sizes_dic = {}
ex_cover_size = {}
co_cover_size = {}
gene_pairs_ct = {}

ws_pv_dic = {}  # WeSME
jaccard_dic = {}  # jaccard index

# set nr of perm
init_pair_num=10**3
min_rank=100
max_pair_num=10**6

for ct in cantypes:
    # we will use this dic later
    ws_ex_cover_sizes_dic[ct] = {}

    # set ct specific directories and files
    # if file is given, replace PAN for cancertype
    if mut_file is not None:
        mut_file_ct = mut_file.replace("PAN",ct)
    else:
        mut_dir_ct = mut_dir + ct + "/"
        mut_file_ct = mut_dir_ct + "_".join([ct, mtype]) + "_list.txt"

    cur_wrs_dir_ct = cur_wrs_dir +  "/" + ct + "/" + mtype + "/" + str(nws) + "/"

    # read mutation profile
    if not os.path.exists(mut_file_ct):
        break
    alt_list_dic[ct], alt_samples[ct] = \
    io_utils.read_mut_list(mut_file_ct, samples=True)
    gene_ks[ct] = dict([(gene, len(alt_list_dic[ct][gene])) \
           for gene in alt_list_dic[ct]])  # gene -> cover size dic
    nsamples[ct] = len(alt_samples[ct])

    # read weighted sampling files for all k
    # use the weighted sampling files in cur_wrs_dir/ct
    # Changed: do this for each cancer type, make dictonary of dictonaries
    kfiles[ct] = os.listdir(cur_wrs_dir_ct)
    # store weighted sampling for each k: k --> list of arrays

    ws_k_cover_dic[ct] = {}
    for kf in kfiles[ct]:
        # read random file for k
        k = int(kf.split(".")[0].split("_")[1])
        ws_k_cover_dic[ct][k], samples[ct] = \
        io_utils.read_mut_list(cur_wrs_dir_ct + kf, genes=False)

    # collect gene pairs
    genes = alt_list_dic[ct].keys()
    for i in range(len(genes)):
       gene1 = genes[i]
       for j in range(i+1, len(genes)):
           gene2 = genes[j]
           gene_pairs = tuple(sorted([gene1, gene2]))
           covers = [alt_list_dic[ct][gene1], alt_list_dic[ct][gene2]]
           # for each gene pair get ex_cover
           set1 = set(covers[0])
           set2 = set(covers[1])
           #cover_sizes = ws_mut_ex.compute_ex_cover(covers)
           cover_sizes = set1.symmetric_difference(set2)
           cover_sizes_co = set1.intersection(set2)
           ex_cover_size[gene_pairs] = ex_cover_size.get(gene_pairs, 0) + \
           len(cover_sizes)
           co_cover_size[gene_pairs] = co_cover_size.get(gene_pairs, 0) + \
           len(cover_sizes_co)
           ex_cover_sizes_dic.setdefault(gene_pairs, {}).update({ct: [len(c) for c in covers]})


# only gene pairs for which at least 2 samples have only one of these genes
# mutated are included
# only use gene pairs with me count > 2
if come=="me":
    ex_cover_size_flt = {k:v for (k,v) in ex_cover_size.iteritems() if v > 2}
else:
    ex_cover_size_flt = {k:v for (k,v) in co_cover_size.iteritems() if v > 0}

# for each gene pair:
for gp in ex_cover_size_flt:

    ws_ex_cover_sizes_sum=[]

    # for each ct with pair:
    for ct in ex_cover_sizes_dic[gp]:
        # get cover sizes from perm sets
        # re-use if already collected for other gene pair
        # ws_ex_cover_sizes_sum is len(excover) for all permutations        #
        ws_ex_cover_sizes_tmp, ws_ex_cover_sizes_dic[ct] = \
        ws_mut_ex.compute_ws_cover_sizes(ex_cover_sizes_dic[gp][ct],
                                         ws_k_cover_dic[ct], init_pair_num,
                                         ws_ex_cover_sizes_dic[ct],
                                         sort = False)
        # add length(cover sizes) to length(cover sizes) of previous
        # if to slow, look at using numpy
        if len(ws_ex_cover_sizes_sum)>0:
            ws_ex_cover_sizes_sum = map(add, ws_ex_cover_sizes_sum,
                                 ws_ex_cover_sizes_tmp[0:init_pair_num])
        else:
            ws_ex_cover_sizes_sum = ws_ex_cover_sizes_tmp[0:init_pair_num]

    while True:  # dynamic sampling
        # compute pvalue
        ws_rank = ws_mut_ex.compute_rank(ex_cover_size[gp],
                                         ws_ex_cover_sizes_sum, me_co=come,
                                         ordered=False)
        pair_num = len(ws_ex_cover_sizes_sum)
        if ws_rank < min_rank and pair_num < max_pair_num:  # if not precise enough
            ws_ex_cover_sizes_sum=[]
            for ct in ex_cover_sizes_dic[gp]:
                ws_ex_cover_sizes_tmp, ws_ex_cover_sizes_dic[ct] = \
                    ws_mut_ex.compute_ws_cover_sizes(ex_cover_sizes_dic[gp][ct],
                                             ws_k_cover_dic[ct], pair_num*10,
                                             ws_ex_cover_sizes_dic[ct],
                                             sort = False)
                # add length(cover sizes) to length(cover sizes) of previous
                # if to slow, look at using numpy
                if len(ws_ex_cover_sizes_sum)>0:
                    ws_ex_cover_sizes_sum = map(add, ws_ex_cover_sizes_sum,
                                 ws_ex_cover_sizes_tmp[0:pair_num*10])
                else:
                    ws_ex_cover_sizes_sum = ws_ex_cover_sizes_tmp[0:pair_num*10]

        else:
            ws_pv = ws_rank/float(pair_num)
            break

    if ws_pv > pth:  # if pvalue is not significant enough
        continue
    ws_pv_dic[gp] = ws_pv # WeSME pvalue

    # calc pv
    # if pv < threshold: redo with larger nr of perm

# write the results
rows = []
for gp in ws_pv_dic:
    rows.append((gp[0], gp[1], -1, ws_pv_dic[gp]))
    
# write the PAN mut file
# TODO: remove double samples
# if mut file doesn't exist, create one
if not os.path.exists(mut_file):
    if not os.path.exists(os.path.dirname(mut_file)):
        os.makedirs(os.path.dirname(mut_file))
            
    list_dic_PAN = {}
    samples_PAN = []
    list_dic_ct = {}
    startix = 0
    for ct in cantypes:
        len_ct=len(alt_samples[ct])
        samples_PAN.extend(alt_samples[ct])
        # change indices in alt_list_dic: increase with startix
        for g in alt_list_dic[ct]:            
            list_dic_ct[g]=[ix+startix for ix in alt_list_dic[ct][g]] 
            if g not in list_dic_PAN:
                list_dic_PAN[g] = list_dic_ct[g]
            else:
                list_dic_PAN[g].extend(list_dic_ct[g])    
        startix += len_ct
        
    io_utils.write_alt_list(list_dic_PAN, filename = mut_file , 
                                samples = samples_PAN)

# write the results
all_pvs = pandas.DataFrame(data=rows,
                           columns=["gene1", "gene2", "jaccard index",
                                    "pv (ws)"])
all_pvs.to_csv(pv_file, sep="\t", index=False)
