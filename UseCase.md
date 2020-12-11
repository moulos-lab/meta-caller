# Use case



1. **Pre-processing**

The files must be .bam format. Using the metacaller with real data sets, we saw that chromosomes other the 1...22 and X,Y could cause spp to crash, so they should be removed beforehand. Also the chromosomes should be in the format of 1,2,...,Y and not CHR1,CHR2,...,CHRY or chr1,chr2,...,chrY. Finally the bigger file must be downsampled since some of the callers don't do that before the peak calling (homer, PeakRanger, spp) but even macs2 that by default downsamples the bigger one, it does that by linear scaling, which gives different results (only Q's peaks are the same).

2. **Options**

The help menu of metacaller is the following:

```
usage: meta-caller [options]

Combined p-value based on 5 peak callers p-values

optional arguments:
  -h, --help            show this help message and exit
  --name NAME           A name for the project
  --macs2_path MACS2_PATH
                        Path to the macs2 executable
  --peakranger_path PEAKRANGER_PATH
                        Path to peakranger executable
  --Q_path Q_PATH       Path to Q executable
  --homer_path HOMER_PATH
                        Path to homer executable

files (REQUIRED):
  -t TREATMENT, --treatment TREATMENT
                        Treatment file, must be BAM format
  -c CONTROL, --control CONTROL
                        Control file, must be BAM format

filters:
  --mnl MNL             Minimum length accepted for a peak (>=)
  --mxl MXL             Maximum length accepted for a peak (<=)
  -l LENGTH, --length LENGTH
                        Extension length, default = 200

output-options:
  -p P                  P-value cut off for the peaks to be reported
  -s                    If simes method calculation will be reported
  -f                    If fisher's method calculation will be reported
  --keep {1,2,3}        Set --keep 1, if you wish to keep each peak caller's
                        results, set --keep 2 if you want a report with all
                        the peaks. By default --keep is set to 3 and only
                        peaks with p-value < than the cut-off will be
                        reported.
                        
```

Only the treatment (-t) and control (-c) files are required, all the other options are optional. If we don't provide a name for the project by default will be 'NA', but be careful if we run the meta-caller again with different arguments, and again we don't assign a name to the project the initial files will be overwritten. It is advised to use the --name argument. 

### simple use:

``` ./metacaller.py -t treat.bam -c control.bam ```

### advanced use:

- path arguments

meta-caller gives us the option, if we have the different peak callers installed but not added to the path, to give their respective paths as arguments 

e.g ``` ./metacaller.py -t treat.bam -c control.bam --name exampleuse --Q_path ~/Desktop/tools/Q/bin/Q```


- filter arguments

by default meta-caller keeps only peaks >=50 & <=1000bp during the filtering. This process is applied only on macs2, PeakRanger and spp** since homer's and Q's peaks have fixed length. These lengths can change, using --mnl and --mxl. Finally each peak by default is extended to 200bp (or trimmed), this can also change using the argument --l and a positive even integer. 

- output options

since we run the peak callers with a very loose p-value (0.9), a lot of peaks and noise are going to be reported. By default in the final report the cut-off p-value is set to 0.05, but it can be adjusted using the -p argument. Also the user can choose if he wants a report with the combined p-value using the fisher's method (-f) or/and simes method(-s). Finally there are 3 options about the files to be saved, by default (--keep 3) only a file with the following columns

|chr|start|end|metacaller's p-value|name|scores|fdr|
|---|-----|---|--------------------|----|------|---|

and if -f and -s are choosen a similar file for each method. The --keep 2 option, also saves a file with all the peaks regardless of their significance (p-value > cut-off). And the --keep 1 option, keeps along with all the previous and the initial files from each peak caller. 

 ```./metacaller.py -t treat.bam -c control.bam --name advanceduse -p 0.001 -s -f --mnl 100 --mxl 2000 --keep 1```

