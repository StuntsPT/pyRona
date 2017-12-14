# Usage

This section describes how to use *pyRona*. All options as of version 0.2.0 are described here.

## Parameters

* Getting help
    - Simply calling the program, or issuing the single argument `-h` will print all the available options on the console

* I/O options
    - *pyRona* **requires** several files as input, that should be specified inn the following options:
        - `-pc`: File containing current environmental variables (formatted as per *BayPass* input)
        - `-fc`: File containing future environmental variables (formatted the same as the present conditions file)
        - `-pop`: File with names of populations (one population per line)
        - `-beta`: *BayPass* output "summary_betai"
        - `-pij`: *BayPass* output "PIJ"
    - Furthermore, *pyRona* also requires an output to be specified, for saving the RONA plot.
        - `-out`: File where the plot will be saved. This option is extension aware, and entering the extension as "PDF", "SVG" or "PNG", will make *pyRona* save the plot in the respective format.
* Parameters
    - *pyRona* requires some parameters to be set in order to perform the analysis. These are:
        - `-bf`: Bayes Factor threshold. This is the value above which associations are considered significant.
        - `-covars`: [**optional**] Number of covars to calculate the RONA for (default: 3)
        - `-outliers`: [**optional**] Number of outliers to remove - "0" skips outlier removal, "1" removes a maximum of 1 outlier (if there is one), "2" removes any number of markers considered outliers (default: 2)
        - `-immutables`: Number of covariates to skip from the environmental variables file. Usefull to skip variables that are the same in te present covars and future covars file, such as latitude (default: 3).
        - `-ronatype`: Defines the RONA is calculated. `absdiff` performs calculations as described in Rellstab et al. 2016 - using the absolute value of the differences, `diff` uses the differences without modulus, and `dist` accounts simply for the distance between the future condition and the trendline (default: `absdiff`).
* Other options (**all optional**)
    - *pyRona* allows for setting some further miscellaneous options:
        - `-no-plots`: Do not draw the individual regression plots.
        - `-no-weighted-means`: Use this option if you wish to use *means* instead of *weighted means* for the RONA calculation.


## Example run:

```
pyRona -pc Popfiles/ENVFILE -fc Popfiles/ENVFILE_rpc26 -be Analyses/Baypass/mcmc_aux/Qsuber_GBS_07_loki_mcmc_aux_summary_betai.out -pij Analyses/Baypass/mcmc_aux/Qsuber_GBS_07_loki_mcmc_aux_summary_pij.out -pop Popfiles/popnames_single_GEO.txt -bf 15 -outliers 0 -out ~/Desktop/rpc26.pdf -covars 4
```

This command will execute *pyRona* using the input files specified by `-pc`, `-fc`, `-be`, `-pij` and `-pop`, using the value `15` as the BF threshold, skipping outlier removal and calculating the RONA value for the most frequent 4 associated environmental variables. The plot will be saved as a PDF under `~/Desktop/rpc26.pdf`.
