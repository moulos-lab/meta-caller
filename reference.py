import pandas as pd
import os, sys
import subprocess
import tempfile
from multiprocessing import Pool

cwd = os.getcwd()

def check_adjusted_res(name, **callers):
    with tempfile.NamedTemporaryFile(mode = 'w', delete = False) as all_bed:
        for caller in callers.keys():
            f1 = f'{cwd}/adjusted_merged_{caller}.bed'
            if os.path.isfile(f1):
                with open(f1) as infile:
                    for line in infile:
                        all_bed.write(line)
            else:
                sys.exit(f'File from {caller} missing.')

    with tempfile.NamedTemporaryFile(mode = 'w', delete = False) as sorted_bed:
        sort_pr = subprocess.Popen(['sort', '-k1,1', '-k2,2n', all_bed.name], stdout = sorted_bed, stderr = subprocess.DEVNULL)
        sort_pr.wait()

    with tempfile.NamedTemporaryFile(mode = 'w', delete = False) as merged:
        merge_pr = subprocess.Popen(['mergeBed', '-i', sorted_bed.name ,'-c', '4,4,5', '-o', 'count,collapse,collapse'], stdout = merged)
        merge_pr.wait()

    os.remove(all_bed.name)
    # to dataframe
    ref_set = pd.read_table(merged.name, names=['chr', 'start', 'end', 'count', 'name', 'scores'], dtype = {'chr' : object})
    ref_set['i'] = ref_set.index
    ref_set['id'] = "Ref_" + ref_set['i'].apply(str)
    ref_set.drop(columns=['i'],inplace = True)
    filename = f'{cwd}/{name}.bed'
    ref_set.to_csv(filename, sep = "\t", index = False)
    for caller in callers.keys():
        os.remove(f'{cwd}/adjusted_merged_{caller}.bed')
    os.remove(merged.name)
    os.remove(sorted_bed.name)
    return filename

def calculate(func, args):
    result = func(*args)
    return result

def get_pvalues(filename, name, caller):
    df = pd.read_table(filename, dtype = {'chr' : object})
    df_with_caller = df[df['name'].str.contains(caller)].copy()
    df_without_caller = df[~df['name'].str.contains(caller)].copy()
    #for peaks only in reference_set we set it to:99
    df_without_caller[f'{caller}'] = 99
    #for peaks in both sets we set it to:pvalue
    df_split_name = df_with_caller['name'].str.split(",", expand = True)
    df_split_pvalue = df_with_caller['scores'].str.split(",", expand = True)
    num_of_col = len(list(df_split_name.columns))
    num_of_rows = df_with_caller.shape[0]
    df_pvalues = pd.DataFrame(index=df_with_caller.index, columns = ['name', 'pvalue'])
    for i in range(num_of_col):
        df_con = pd.concat([df_split_name.loc[:,i], df_split_pvalue.loc[:,i]], axis = 1, keys =['name', 'pvalue'])
        df_con.dropna(inplace = True)
        df_found = df_con[df_con['name'].str.contains(caller)]
        indexes = df_found.index.tolist()
        for ind in indexes:
            if df_pvalues.loc[ind].isnull().all():
                df_pvalues.loc[ind] = df_found.loc[ind,['name', 'pvalue']]
            else:
                #keep the smaller pvalue for each caller
                if df_pvalues.loc[ind,'pvalue'] > df_found.loc[ind,'pvalue']:
                    df_pvalues.loc[ind] = df_found.loc[ind,['name', 'pvalue']]
    df_pvalues['pvalue'] = df_pvalues['pvalue'].astype(float)
    df_with_caller[f'{caller}'] = df_pvalues['pvalue'].copy()
    df_with_caller = df_with_caller[['id', caller]]
    df_without_caller = df_without_caller[['id', caller]]
    final_df = pd.concat([df_with_caller, df_without_caller])
    final_df.sort_index(inplace = True)
    return final_df

def pool_ref_matrix(filename,cores,name,**callers):
    df = pd.read_table(filename, usecols = ['id'])
    with Pool(processes=cores) as pool:
        TASKS = [(get_pvalues, (filename, name ,caller)) for caller in callers.keys()]

        results = [pool.apply_async(calculate, t) for t in TASKS]
        for caller,r in zip(callers.keys(),results):
            df = pd.concat([df, r.get()[caller]], axis = 1)
    matrix_file = f"{cwd}/{name}_callers.txt"
    df.to_csv(matrix_file, sep = "\t", index = False)
