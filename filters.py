mergeBedimport pandas as pd
import os
import subprocess
import tempfile
import shutil

cwd = os.getcwd()

def filter_Q(summitsf,args):
    summits = pd.read_table(summitsf, usecols = ['Chromosome', 'pos', 'pos+1', 'p-value'])
    #adjust start and end positions
    ext_len = int(args.length / 2)
    summits['pos+1'] = summits['pos'] + ext_len
    summits["pos"] = summits["pos"] - ext_len
    # add 'id' column for the final name
    summits['id'] = summits.index
    summits['name'] = "Q_" + summits['id'].apply(str)
    # convert-log10(pvalue)
    summits['p-value'] = 10**(-summits['p-value'])
    #sort by chrom and start position
    summits.sort_values(['Chromosome', 'pos'], ascending = [True, True], inplace = True)
    # write adjsted to 200 length peaks
    summits = summits[['Chromosome', 'pos', 'pos+1', 'name', 'p-value']]
    with tempfile.NamedTemporaryFile() as adjusted_peaks:
        summits.to_csv(adjusted_peaks.name, sep = "\t", index = False, header = False)
        merge(adjusted_peaks.name, 'Q')

def filter_macs2(peaks,summits,args):
    peaks = pd.read_table(peaks, names=['chr', 'start', 'end', 'name', 'score', 'strand', 'fold_enr', 'pvalue', 'qvalue', 'summit'], dtype={'chr' : object})
    summits = pd.read_table(summits, names=['chr', 'sum', 'sum+1', 'name', 'score'], dtype={'chr' : object})
    #add 'length' column
    peaks['length'] = peaks['end'] - peaks['start']
    #keep peaks >= 50 & <=1000
    filtered_peaks = peaks[(peaks['length'] >= args.mnl) & (peaks['length'] <= args.mxl)]
    #keep summits based on peaks filtered by length
    common = pd.merge(summits, filtered_peaks, on=['name'])
    filtered_summits = summits[summits['name'].isin(common['name'])].copy()
    ext_len = int(args.length / 2)
    filtered_summits['sum+1'] = filtered_summits['sum'] + ext_len
    filtered_summits["sum"] = filtered_summits["sum"] - ext_len
    #convert -log10(pvalue)
    filtered_summits['pvalue'] = 10**(-peaks['pvalue'])
    #keep only the 5 columns
    filtered_summits = filtered_summits[['chr', 'sum', 'sum+1', 'name', 'pvalue']]
    #sort by chrom and start position
    filtered_summits.sort_values(['chr', 'sum'], ascending = [True, True], inplace = True)
    # write adjusted peaks
    with tempfile.NamedTemporaryFile() as adjusted_peaks:
        filtered_summits.to_csv(adjusted_peaks.name, sep = "\t", index = False, header = False)
        merge(adjusted_peaks.name, 'macs2')

def filter_homer(peaks,args):
    #tmp file
    only_peaks = tempfile.NamedTemporaryFile()
    #remove comment lines
    subprocess.run(['grep', '^[^#]', peaks], stdout = only_peaks)
    peaks = pd.read_table(only_peaks.name, names=['id', 'chr', 'start', 'end', 'strand', 'nomtagcount', 'focusratio', 'score', 'Totaltags', 'ControlTags', 'FoldChVsCon', 'pvalue_VS_control', 'clonalfoldch', 'pvalue_vs_local'], dtype = {'chr': object})
    #find summit start&end position
    peaks['length'] = peaks['end'] - peaks['start']
    peaks['start'] = peaks['start'].add(peaks['length'].floordiv(2))
    peaks['end'] = peaks['start'] + 1
    final = peaks[['chr']].copy()
    ext_len = int(args.length / 2)
    final['start'] = peaks['start'].sub(ext_len)
    final['end'] = peaks['start'].add(ext_len)
    final['name'] = "homer_" + peaks['id']
    final['pvalue'] = peaks['pvalue_VS_control']
    #sort by chrom and start position
    final.sort_values(['chr', 'start'], ascending = [True, True], inplace = True)
    # write 'adjusted to 200 length' peaks
    with tempfile.NamedTemporaryFile() as adjusted_peaks:
        final.to_csv(adjusted_peaks.name, sep = "\t", index = False, header = False)
        merge(adjusted_peaks.name, 'homer')

def filter_PeakRanger(peaks,details,args):
    #tmp file for peaks
    temp_peaks = tempfile.NamedTemporaryFile()
    #remove comment lines
    subprocess.run(['grep', '^[^#]', peaks], stdout = temp_peaks)
    peaks = pd.read_table(temp_peaks.name, names=['chr', 'start', 'end', 'region_id', 'score', 'strand'], dtype={'chr' : object})
    #tmp file for details
    temp_details = tempfile.NamedTemporaryFile()
    subprocess.run(['grep', '^[^#]', details], stdout = temp_details)
    details = pd.read_table(temp_details.name, names=['chr', 'start', 'end', 'nearby_genes', 'region_id', 'summit','pvalue', 'fdr', 'strand', 'treads', 'creads'], dtype={'chr' : object})
    # add 'length'column
    peaks['length'] = peaks['end'] - peaks['start']
    #add column 'id' for the final name of peak
    details['id'] = details.index
    filtered_peaks = peaks[(peaks['length'] >= args.mnl) & (peaks['length'] <= args.mxl)]
    common = pd.merge(details, filtered_peaks, on=['region_id'])
    filtered_details = details[details['region_id'].isin(common['region_id'])].copy()
    #adjust start and end positions
    ext_len = int(args.length / 2)
    filtered_details['end'] = filtered_details['summit'] + ext_len
    filtered_details["start"] = filtered_details["summit"] - ext_len
    filtered_details['name'] = "PeakRanger_" + filtered_details['id'].apply(str)
    #keep 5 columns and sort by score column
    filtered_details = filtered_details[['chr', 'start', 'end', 'name', 'pvalue']]
    #sort by chrom and start position
    filtered_details.sort_values(['chr', 'start'], ascending = [True, True], inplace = True)
    # write adjsted to 200 length peaks
    with tempfile.NamedTemporaryFile() as adjusted_peaks:
        filtered_details.to_csv(adjusted_peaks.name, sep = "\t", index = False, header = False)
        merge(adjusted_peaks.name, 'PeakRanger')

def filter_spp(summitsf,args):
    summits = pd.read_table(summitsf, usecols = ['chr', 'pos', 'FDR'], dtype={'pos':int})
    #adjust start and end positions
    ext_len = int(args.length /2)
    summits['pos+1'] = summits['pos'] + ext_len
    summits['pos'] = summits['pos'] - ext_len
    summits['id'] = summits.index
    summits['name'] = "spp_" + summits['id'].apply(str)
    #keep only the 5 columns
    summits = summits[['chr', 'pos', 'pos+1', 'name', 'FDR']]
    #sort by chrom and start position
    summits.sort_values(['chr', 'pos'], ascending = [True, True], inplace = True)
    # write adjsted to 200 length peaks
    with tempfile.NamedTemporaryFile() as adjusted_peaks:
        summits.to_csv(adjusted_peaks.name, sep = "\t", index = False, header = False)
        merge(adjusted_peaks.name, 'spp')

def merge(peaks, caller):
    with tempfile.NamedTemporaryFile(delete = False) as merged:
        subprocess.run(['mergeBed', '-i', peaks, '-c', '4,5','-o','collapse,collapse'], stderr = subprocess.STDOUT ,stdout = merged)
        re_adjust_peaks(merged.name, caller)

def re_adjust_peaks(peaks,caller):
    merged_peaks = pd.read_table(peaks, names = ['chr', 'start', 'end', 'name', 'pvalue'], dtype = {'chr' : object})
    #keep only the peaks that are merged
    more_than_1_peaks = merged_peaks[merged_peaks['name'].str.contains(",")].copy()
    if more_than_1_peaks.shape[0] == 0:
        shutil.copy(peaks, f'{cwd}/adjusted_merged_{caller}.bed')
        os.remove(peaks)
    else:
        only_1_peak = merged_peaks[~merged_peaks['name'].str.contains(",")]
        #add a 'length' column so we can 'trim' accordingly start&end positions
        more_than_1_peaks['length'] = more_than_1_peaks['end'] - more_than_1_peaks['start']
        #calculate the bases to trim/add
        more_than_1_peaks['ext_len'] = (more_than_1_peaks['length'] - 200).floordiv(2)
        more_than_1_peaks['start'] = more_than_1_peaks['start'] + more_than_1_peaks['ext_len'] + (more_than_1_peaks['length']-200).mod(2)
        more_than_1_peaks['end'] = more_than_1_peaks['end'] - more_than_1_peaks['ext_len']
        # drop additional columns from DataFrame
        more_than_1_peaks.drop(columns = ['ext_len', 'length'], inplace = True)
        def split_pvalue_name(df,caller):
            caller_num_letters = len(caller)
            #pvalue part
            df2 = df[['pvalue']].copy()
            df2['ext_pvalue'] = df2['pvalue'].str.split(",")
            df2['min'] = [min(x) for x in df2['ext_pvalue']]
            df['pvalue'] = df2['min']
            #name part
            df3 = df[['name']].copy()
            df3['rm_caller'] = df3['name'].str.replace(caller,"")
            df3['ext_name'] = df3['rm_caller'].str.split(",")
            df3['new_name'] = ["".join(x) for x in df3['ext_name']]
            df3['final_name'] = caller + df3['new_name'].apply(str)
            df['name'] = df3['final_name']
            return df
        #for the final dataframe
        final_set = pd.concat([only_1_peak, split_pvalue_name(more_than_1_peaks, caller)])
        final_set.sort_index(inplace = True)
        filename = f"{cwd}/adjusted_merged_{caller}.bed"
        final_set.to_csv(filename, sep="\t", index = False, header = False)
        os.remove(peaks)
