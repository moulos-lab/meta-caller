#!/usr/bin/env python3

import sys
import os
import shutil
import subprocess
from multiprocessing import Pool
import argparse
from filters import filter_Q, filter_macs2, filter_homer, filter_PeakRanger, filter_spp
from reference import check_adjusted_res, pool_ref_matrix
from combine import weighted, simes, fishers
cwd = os.getcwd()

def get_arguments():
    parser = argparse.ArgumentParser(description='Combined p-value based on 5 peak callers p-values', prog='meta-caller', usage = '%(prog)s [options]')

    #file arguments
    input_f = parser.add_argument_group('files (REQUIRED)')
    input_f.add_argument('-t', '--treatment', required = True, help = 'Treatment file, must be BAM format', type = argparse.FileType('r'))
    input_f.add_argument('-c', '--control', required = True, help = 'Control file, must be BAM format', type = argparse.FileType('r'))
    #optional arguments
    parser.add_argument('--name', type = str, help = 'A name for the project', default = 'NA')
    parser.add_argument('--cores', type = int, default = 1, help = 'Number of cores to be used')
    #callers
    parser.add_argument('--macs2_path', type = str, help = 'Path to the macs2 executable')
    parser.add_argument('--peakranger_path', type = str, help = 'Path to peakranger executable')
    parser.add_argument('--Q_path', type = str, help = 'Path to Q executable')
    parser.add_argument('--homer_path', type = str, help = 'Path to homer executable' )
    #filtering arguments
    filter_o = parser.add_argument_group('filters')
    filter_o.add_argument('--mnl', type = int, help = 'Minimum length accepted for a peak (>=)', default = 50)
    filter_o.add_argument('--mxl', type = int, help = 'Maximum length accepted for a peak (<=)', default = 1000)
    filter_o.add_argument('-l', '--length', type = int, help = 'Extension length, default = 200', default = 200)
    #output arguments
    output = parser.add_argument_group('output-options')
    output.add_argument('-p', type = float, help = 'P-value cut off for the peaks to be reported', default = 0.05)
    output.add_argument('-s', help = 'If simes method calculation will be reported', action = 'store_true')
    output.add_argument('-f', help = 'If fisher\'s method calculation will be reported', action = 'store_true')
    output.add_argument('--keep', type = int, help = 'Set --keep 1, if you wish to keep each peak caller\'s results, set --keep 2 if you want a report with all the peaks. By default --keep is set to 3 and only peaks with p-value < than the cut-off will be reported.', default = 3, choices = [1,2,3])

    if len(sys.argv) == 1:
        parser.print_help(file=sys.stderr)
        sys.exit(1)

    args = parser.parse_args()
    if args.mnl < 5:
        sys.exit('Error: Minimum length allowed is 5')
    if args.mxl > 10000:
        sys.exit('Error: Maximum length allowed is 10.000')
    if args.mnl > args.mxl:
        sys.exit('Error: Minimum length cannot be bigger than maximum')
    if args.length%2!=0 or args.length < 0:
        sys.exit('Error: Length must be a positive even number')
    return args

def check_bam_files():
    print('----------------------')
    print('Checking files:\n')
    treat = args.treatment.name
    if not os.path.isfile(treat) or not treat.endswith(".bam"):
        sys.exit('Treatment file does not exist or does not have the correct .bam extension, please check again.')
    else:
        print(f'\t Treatment file: {treat}')
    control = args.control.name
    if not os.path.isfile(control) or not control.endswith(".bam"):
        sys.exit('Control file does not exist or does not have the correct .bam extension, please check again')
    else:
        print(f'\t Control file: {control}')

def get_tool_path(args_path):

    if not shutil.which(args_path):
        sys.exit(f'Error:{args_path} not the correct path')
    else:
        print(f'\t {args_path.split("/")[-1]}\'s path is correct')

def check_dependencies():
    callers = {}
    print('----------------------')
    print('Checking Dependencies:\n')
    print('----------------------')
    print('Checking peak-callers:\n')

    if args.macs2_path == None :
        try:
            macs2_version = subprocess.check_output(['macs2', '--version'], encoding = 'UTF-8')
            print(f'\t {macs2_version.rstrip()} is installed')
            callers['macs2'] = 'macs2'
        except:
            sys.exit('Error: macs2 is either not installed or not added to the path, use \'--macs2_path to specify the path to executable')
    else:
        callers['macs2'] = args.macs2_path
        get_tool_path(args.macs2_path)

    if args.Q_path == None :
        try:
            Q_version = subprocess.check_output(['Q', '--version'], encoding = 'UTF-8')
            print(f'\t {Q_version.rstrip()} is installed')
            callers['Q'] = 'Q'
        except:
            sys.exit('Error: Q is either not installed or not added to the path, use \'--Q_path\' to specify the path to the executable')
    else:
        callers['Q'] = args.Q_path
        get_tool_path(args.Q_path)

    if args.peakranger_path == None:
        try:
            PeakRanger_version = subprocess.check_output(['peakranger','ranger', '--version'], encoding = 'UTF-8')
            print(f'\t PeakRanger version {PeakRanger_version.strip()} is installed')
            callers['PeakRanger'] = 'peakranger'
        except:
            sys.exit('Error: PeakRanger is either not installed or not added to the path, use \'--peakranger_path\' to specify the path to the executable')
    else:
        callers['PeakRanger'] = args.peakranger_path
        get_tool_path(args.peakranger_path)

    if args.homer_path == None:
        try:
            subprocess.run(['findPeaks'], stdout = subprocess.DEVNULL, stderr = subprocess.DEVNULL)
            print('\t homer is installed')
            callers['homer'] = 'homer'
        except:
            sys.exit('Error: Homer is either not installed or not added to the path, use \'--homer_path\' to specify the path to the executable')
    else:
        callers['homer'] = args.homer_path
        get_tool_path(args.homer_path)

    try:
        R_version = subprocess.run(['R', '--version'], stdout = subprocess.DEVNULL)
        if R_version.returncode == 0:
            try:
                spp_version = subprocess.check_output(['Rscript', f'{cwd}/check_spp.r'], encoding = 'UTF-8')
                print(f'\t spp {spp_version.split("‘")[1].split("’")[0]} version is installed')
                callers['spp'] = 'spp'
            except:
                sys.exit('Error: R package spp is not installed, try \'githubinstall("spp", ref = "88bbd37")\'')
    except:
        sys.exit('Error: R is not installed')

    print('----------------------')
    print('Checking other tools:\n')
    try:
        samtools_version = subprocess.check_output(['samtools', '--version'], encoding = 'UTF-8')
        print(f'\t {samtools_version.splitlines()[0]} version is installed')
    except:
        sys.exit('Error: samtools is not installed, please install it and run the program again')

    try:
        bedtools_version = subprocess.check_output(['bedtools', '--version'], encoding ='UTF-8')
        print(f'\t {bedtools_version.rstrip()} is intalled')
    except:
        sys.exit('Error: bedtools is not installed, please install it and run the program again')

    return callers

def macs2_run(treat, control, name):
    macs2_log = open(f'{cwd}/macs2/{name}.log', 'w')
    macs2_pr = subprocess.run([callers['macs2'], 'callpeak', '-t', treat, '-c', control, '-n', 'macs2', '-p', '0.9', '--outdir', f'{cwd}/macs2/'],  stdout = macs2_log ,stderr = macs2_log)
    macs2_log.close()
    print('\t macs2 completed')

def Q_run(treat, control, name, ncores):
    Q_log = open(f'{cwd}/Q/{name}.log', 'w')
    Q_pr = subprocess.run([callers['Q'], '-t', treat, '-c', control, '-n', '10000000', '-o', f'{cwd}/Q/{name}','-p', f'{ncores}', '-v'], stdout = Q_log, stderr = Q_log)
    Q_log.close()
    print('\t Q completed')

def PeakRanger_run(treat, control, name, ncores):
    PeakRanger_log = open(f'{cwd}/PeakRanger/{name}.log', 'w')
    Peak_pr = subprocess.run([callers['PeakRanger'], 'ranger', '-d', treat, '-c', control, '--format', 'bam', '-l', '200', '-o', f'{cwd}/PeakRanger/{name}', '-p', '0.9', '-t', f'{ncores}', '--verbose'], stdout = PeakRanger_log, stderr = PeakRanger_log)
    PeakRanger_log.close()
    print('\t PeakRanger completed')

def homer_run(treat, control, name):
    homer_log = open(f'{cwd}/homer/{name}.log', 'w')
    homer_mktag = "makeTagDirectory".join(str(callers['homer']).rsplit("homer", 1))
    homer_findPeaks = "findPeaks".join(str(callers['homer']).rsplit("homer", 1))
    samtools_sam_t = subprocess.run(['samtools', 'view', '-h', treat, '-o', 'treat.sam'])
    samtools_sam_p = subprocess.run(['samtools', 'view', '-h', control, '-o', 'control.sam'])
    homer_t_mkt = subprocess.run([homer_mktag, f'{cwd}/homer/treatment', 'treat.sam', '-format', 'sam'], stderr = subprocess.DEVNULL)
    homer_c_mkt = subprocess.run([homer_mktag, f'{cwd}/homer/control', 'control.sam', '-format', 'sam'], stderr = subprocess.DEVNULL)
    os.remove('treat.sam')
    os.remove('control.sam')
    homer_t_fp = subprocess.run([homer_findPeaks, f'{cwd}/homer/treatment/', '-style', 'factor', '-o', 'auto', '-i', f'{cwd}/homer/control/', '-poisson', '0.9', '-F', '2', '-P', '1', '-L', '1', '-LP', '1', '-C', '0'], stdout = homer_log, stderr = homer_log)
    homer_log.close()
    print('\t homer completed')

def spp_run(treat, control, name, ncores):
    spp_log = open(f'{cwd}/spp/{name}.log', 'w')
    spp_pr = subprocess.run(['Rscript', f'{cwd}/spp.r', treat, control, name, f'{ncores}'], stdout = spp_log, stderr = spp_log)
    spp_log.close()
    print('\t spp completed')

def calculate(func, args):
    result = func(*args)
    return result

def pool_run_peak_callers():
    treat = args.treatment.name
    control = args.control.name
    name = args.name
    mnl = args.mnl
    mxl = args.mxl
    nc = args.cores
    make_dirs()
    print('----------------------')
    print('Peak Calling:\n')
    with Pool(processes = nc) as pool:
        if nc <= 5:
            ncores = 1
            spp_cores = ncores
        else:
            nc = nc - 5
            if nc%3 == 0:
                ncores = int(nc/3) + 1
                spp_cores = ncores
            else:
                mod_cores = nc%3
                ncores = int(nc/3) + 1
                spp_cores = ncores + mod_cores


        TASKS = [(homer_run, (treat, control, name))] + \
                [(spp_run, (treat, control, name, spp_cores))] + \
                [(macs2_run, (treat, control, name))] + \
                [(Q_run, (treat, control, name, ncores))] + \
                [(PeakRanger_run, (treat, control, name, ncores))]

        results = [pool.apply_async(calculate, t) for t in TASKS]
        for r in results:
            r.get()

    args.treatment.close()
    args.control.close()
    filter_peaks()

def filter_peaks():
    print('----------------------')
    print('Filtering Peaks:\n')
    filter_macs2(f'{cwd}/macs2/macs2_peaks.narrowPeak', f'{cwd}/macs2/macs2_summits.bed',args)
    print('\t macs2 completed')
    filter_Q(f'{cwd}/Q/{args.name}-Q-summit-info.tab',args)
    print('\t Q completed')
    filter_PeakRanger(f'{cwd}/PeakRanger/{args.name}_region.bed', f'{cwd}/PeakRanger/{args.name}_details',args)
    print('\t PeakRanger completed')
    filter_homer(f'{cwd}/homer/treatment/peaks.txt', args)
    print('\t homer completed')
    filter_spp(f'{cwd}/spp/{args.name}.binding.positions.txt', args)
    print('\t spp completed')
    if args.keep != 1:
        rm_dirs()

def make_dirs():
    for caller in callers.keys():
        path = os.path.join(cwd, caller)
        if not os.path.isdir(path):
            os.makedirs(path)

def rm_dirs():
    for caller in callers.keys():
        path = os.path.join(cwd, caller)
        shutil.rmtree(path, ignore_errors = False)

if __name__ == '__main__':
    args = get_arguments()
    callers = check_dependencies()
    check_bam_files()
    pool_run_peak_callers()
    f = check_adjusted_res(args.name, **callers)
    pool_ref_matrix(f, args.cores, args.name,**callers)
    weighted(args)
    if args.s == True:
       simes(args)
    if args.f == True :
       fishers(args)
    os.remove(f'{cwd}/{args.name}.bed')
    os.remove(f'{cwd}/{args.name}_callers.txt')