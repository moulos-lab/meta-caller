import pandas as pd
import os, sys
from scipy import stats
import numpy as np
from statsmodels.stats.multitest import fdrcorrection
cwd = os.getcwd()

def weighted(args):
    ref_set_file = f'{cwd}/{args.name}_callers.txt'
    if not os.path.isfile(ref_set_file):
        sys.exit()
    else:
        ref_set = pd.read_table(ref_set_file , dtype = {'id' : object}, na_values = ['99'])

    weights = {'macs2'     : 0.1873048639,
               'Q'         : 0.2988098977,
               'homer'     : 0.1656223856,
               'PeakRanger': 0.1505676135,
               'spp'       : 0.1976952393
               }

    ref_set = ref_set.fillna(1.0)
    ref_set['macs2'] = ref_set['macs2'] ** weights['macs2']
    ref_set['spp'] = ref_set['spp'] ** weights['spp']
    ref_set['homer'] = ref_set['homer'] ** weights['homer']
    ref_set['Q'] = ref_set['Q'] ** weights['Q']
    ref_set['PeakRanger'] = ref_set['PeakRanger'] ** weights['PeakRanger']

    ref_set['meta-caller']  = 1.0
    ref_set['combine2'] = ref_set[['macs2', 'spp', 'homer', 'Q', 'PeakRanger', 'meta-caller']].cumprod(axis=1)['meta-caller']
    ref_set['meta-caller'] = ref_set['combine2'].replace({1 : 99.0})
    ref_set.drop(columns = ['combine2'], inplace = True)

    combine = ref_set[['id', 'meta-caller']].copy()
    #fdr correction
    pvals = pd.Series(combine['meta-caller'].values.tolist())
    fdr = fdrcorrection(pvals)
    combine['fdr'] = fdr[1].tolist()
    #keep peaks with p-value <= than the one given
    combine_p = combine[combine['meta-caller']<= args.p].copy()
    combine_p.sort_values(['meta-caller'], ascending = True, inplace = True)
    #report the start&end positions of these peaks
    peaks_bed = pd.read_table(f'{cwd}/{args.name}.bed', dtype = {'chr': object})
    common = pd.merge(combine_p,peaks_bed, on=['id'])
    common = common[['chr', 'start', 'end','meta-caller', 'name', 'scores', 'fdr']]
    filename = f"{cwd}/{args.name}_metacaller.tsv"
    common.to_csv(filename, sep = "\t", index = False)
    #report all peaks if keep = 1 or 2
    if args.keep != 3 :
        filename2 = f'{cwd}/{args.name}_all_metacaller.tsv'
        common = pd.merge(combine, peaks_bed, on=['id'])
        common = common[['chr', 'start', 'end', 'meta-caller', 'name', 'scores', 'fdr']]
        common.to_csv(filename2, sep = "\t", index = False)

def simes(args):
    ref_set = pd.read_table(f'{cwd}/{args.name}_callers.txt', dtype = {'id' : 'object'})
    ref_simes = ref_set[['macs2', 'PeakRanger', 'Q', 'spp', 'homer']].copy()

    def simes_r(row):
        row.sort_values(inplace = True)
        for j in range(1,5):
            row[j] = row[j] / (j + 1)
        return(np.nanmin(row.values))

    ref_set['simes'] = ref_simes.apply(simes_r, axis = 1)
    ref_set = ref_set[['id', 'simes']]
    #fdr correction
    pvals = pd.Series(ref_set['simes'].values.tolist())
    fdr = fdrcorrection(pvals)
    ref_set['fdr'] = fdr[1].tolist()
    #keep peaks with p-value <= than the one given
    ref_set_p = ref_set[ref_set['simes'] <= args.p].copy()
    ref_set_p.sort_values(['simes'], ascending = True, inplace = True)
    #report the start&end positions of these peaks
    peaks_bed = pd.read_table(f'{cwd}/{args.name}.bed', dtype = {'chr': object})
    common = pd.merge(ref_set_p, peaks_bed, on=['id'])
    common = common[['chr', 'start', 'end', 'simes', 'name', 'scores', 'fdr']]
    filename = f"{cwd}/{args.name}_simes.tsv"
    common.to_csv(filename, sep = "\t", index = False, header = True)
    if args.keep != 3:
        filename2 = f'{cwd}/{args.name}_all_simes.tsv'
        common = pd.merge(ref_set, peaks_bed, on=['id'])
        common = common[['chr', 'start', 'end', 'simes', 'name', 'scores', 'fdr']]
        common.to_csv(filename2, sep = "\t", index = False, header = True)

def fishers(args):
    ref_set = pd.read_table(f'{cwd}/{args.name}_callers.txt', dtype = {'id' : 'object'})
    callers = ref_set[['macs2', 'homer', 'Q', 'spp', 'PeakRanger']]
    fishers = pd.Series(callers.values.tolist())
    ref_set['pvalues'] = fishers
    ref_set['fishers'] = [stats.combine_pvalues(x)[1] for x in ref_set['pvalues']]
    ref_set = ref_set[['id', 'fishers']]
    #fdr correction
    pvals = pd.Series(ref_set['fishers'].values.tolist())
    fdr = fdrcorrection(pvals)
    ref_set['fdr'] = fdr[1].tolist()
    #keep peaks with p-value <= than the one given
    ref_set_p = ref_set[ref_set['fishers'] <= args.p].copy()
    ref_set_p.sort_values(['fishers'], ascending = True, inplace = True)
    #report the start&end positions of these peaks
    peaks_bed = pd.read_table(f'{cwd}/{args.name}.bed', dtype = {'chr': object})
    common = pd.merge(ref_set_p, peaks_bed, on=['id'])
    common = common[['chr', 'start', 'end', 'fishers', 'name', 'scores', 'fdr']]
    filename = f"{cwd}/{args.name}_fishers.tsv"
    common.to_csv(filename, sep = "\t", index = False, header = True)
    if args.keep != 3:
        filename2 = f'{cwd}/{args.name}_all_fishers.tsv'
        common = pd.merge(ref_set, peaks_bed, on=['id'])
        common = common[['chr', 'start', 'end', 'fishers', 'name', 'scores', 'fdr']]
        common.to_csv(filename2, sep ="\t", index = False, header = True)
