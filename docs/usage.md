# Usage
This section describes how to use *pyRona*. It is not yet complete, and for now it is just a shallow copy of the help text.

```
How to get online help:
  -h, --help            show this help message and exit

Input/Output options:
  -pc PRESENT_COVARS_FILE
                        File with Present environmental data.
  -fc FUTURE_COVARS_FILE
                        File with Future environmental data.
  -pop POPNAMES_FILE    File with population names.
  -beta BAYPASS_SUMMARY_BETAI_FILE
                        Baypass summary betai file.
  -pij BAYPASS_PIJ_FILE
                        Baypass pij file.
  -out OUTFILE          Path to where RONA plot should be saved. Supports PDF, SVG and PNG extensions.

Program execution options:
  -bf BAYES_FACTOR      Bayes factor treshold for considering associations.
  -covars NUM_COVARS    Number of covars to calculate the RONA for.
  -outliers {0,1,2}     Number of outliers to remove. 0 does no outier removal, 1 removes **at most** 1 outlier and 2 removes **any** number of outliers that match the distance criteria.
  -immutables IMMUTABLES [IMMUTABLES ...]
                        List of immutable covariates. These are not even parsed from the betai file. By default the first 3 covars are skipped. You can enter any other values here.
  -ronatype {diff,absdiff,dist}
                        Type of RONA to calculate. Default is absolute difference as in Rellstab et al. 2016. Other options are 'difference' (not abs) and 'distance' (future vs. trendline).

Miscellaneous options:
  -no-plots             Pass this option if you don't want individual regression plots to be drawn.
  -no-weighted-means    Pass this option if you don't want to useweighted means for RONA calculations.
```

Example run:

```
Coming soon
```
