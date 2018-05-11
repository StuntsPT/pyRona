# Usage

This section describes how to use *pyRona*. All options as of version 0.3.3 are described here.

## Parameters

* Getting help
    - Simply calling the program, or issuing the single argument `-h` will print all the available options on the console

* Positional arguments:
    - The first argument should be the name of the program whose results you wish to analyse. *pyRona* can take one of two options here:
        - `baypass`:Process Baypass results
        - `lfmm`:Process LFMM results

    - Depending on which you choose, some of the other options will be different.

* I/O options
    - *pyRona* **requires** several files as input, that should be specified in the following options:
        - Common options to both *Baypass* and *LFMM*:
            - `-pc`: File containing current environmental variables, formatted as per the upstream program input (*Baypass* or *LFMM*)
            - `-fc`: File containing future environmental variables, formatted as per the upstream program input (*Baypass* or *LFMM*)
            - `-out`: File where the plot will be saved. This option is extension aware, and entering the extension as "PDF", "SVG" or "PNG", will make *pyRona* save the plot in the respective format
            - `-map`: File where a RONA map will be saved. This options is also extension aware. In order to draw the map, `cartopy` must be installed on your system (see [Installation](install.md)).
        - *Baypass* specific options:
            - `-pop`: File with names of populations (one population per line)
            - `-beta`: *BayPass* output "summary_betai"
            - `-pij`: *BayPass* output "PIJ"
        - *LFMM* specific options:
            - `-assoc`: File containing a matrix of *p*-values as outputted by *LFMM*
        - `-geno`: The input file used for *LFMM*. Should be in the *LFMM* format

* Parameters
    - *pyRona* requires some parameters to be set in order to perform the analysis. These are:
        - Common options to both *Baypass* and *LFMM*:
            - `-covars`: [**optional**] Number of covars to calculate the RONA for (default: 3)
            - `-remove-outliers`: [**optional**] Pass this option if you want to *pyRona* to calculate the regression trend-line *after* removing any coordinates considered outliers
            - `-immutables`: Covariates to skip from the environmental variables file. Useful to skip variables that are the same in te present covars and future covars file, such as latitude (default: 1 2 3)
            - `-ronatype`: Defines the RONA is calculated. `absdiff` performs calculations as described in Rellstab et al. 2016 - using the absolute value of the differences, `diff` uses the differences without modulus, and `dist` accounts simply for the distance between the future condition and the trendline (default: `absdiff`)
        - *Baypass* specific options:
            - `-bf`: Bayes Factor threshold. This is the value above which associations are considered significant
        - *LFMM* specific options:
            - `-P`: *p*-value threshold. This is the value below which associations are considered significant

* Other options (**all optional**)
    - *pyRona* allows for setting some further miscellaneous options:
        - `draw-ind-plots`: Passing this options with a path to a directory will cause *pyRona* to draw each individual "covar-marker" plot. If it is omitted, no individual plots are drawn.
        - `-no-weighted-means`: Set this option to "1" or "True" if you wish to use *means* instead of *weighted means* for the RONA calculation


## Example runs:

```
pyRona baypass -pc Popfiles/ENVFILE -fc Popfiles/ENVFILE_rpc26 -beta Analyses/Baypass/mcmc_aux/Qsuber_GBS_07_loki_mcmc_aux_summary_betai.out -pij Analyses/Baypass/mcmc_aux/Qsuber_GBS_07_loki_mcmc_aux_summary_pij.out -pop Popfiles/popnames_single_GEO.txt -bf 15 -outliers 0 -out ~/Desktop/rpc26.pdf -covars 4
```

This command will execute *pyRona* set to analyse *Baypass* output using the input files specified by `-pc`, `-fc`, `-beta`, `-pij` and `-pop`, using the value `15` as the BF threshold, skipping outlier removal and calculating the RONA value for the most frequent 4 associated environmental variables. The plot will be saved as a PDF under `~/Desktop/rpc26.pdf`.


```
pyRona lfmm -pc ~/aa/present_covars.txt -fc ~/aa/RCP26_covars.txt -out ~/Desktop/LFMM_rpc26.pdf -P 0.01 -assoc ~/aa/scaled_lfmm_results.csv -geno ~/aa/qsuber.lfmm
```

This command will execute *pyRona* set to analyse *LFMM* output using the input files specified by `-pc`, `-fc`, `-assoc` and `-geno`, using the value `0.01` as the *p*-value threshold, skipping outlier removal (default) and calculating the RONA value for the most frequent 3 associated environmental variables (default). The plot will be saved as a PDF under `~/Desktop/LFMM_rpc26.pdf`.
