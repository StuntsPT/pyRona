# Wrapping LFMM

*pyRona* can take the output of *LFMM* as input to calculate the RONA values. Due to this, *pyRona* includes an `R` script named `LFMM_workflow.R` that can be used to automate and wrap *LFMM* analyses.
Not that it is highly recommended that you give *LFMM*'s [manual](https://bcm-uga.github.io/lfmm/reference/index.html), [tutorial](https://bcm-uga.github.io/lfmm/articles/lfmm), and [paper](http://dx.doi.org/10.1093%2Fmolbev%2Fmst063) a through read before using this script.


## LFMM_workflow.R and LFMM2_workflow.R

This script will automate the workflow for the [LFMM](https://bcm-uga.github.io/lfmm/index.html) and [LFMM2](https://github.com/bcm-uga/LEA)
software by Olivier François and Kevin Caye, which is described in [this paper](http://dx.doi.org/10.1093%2Fmolbev%2Fmst063) (LFMM) and [this paper](https://doi.org/10.1093/molbev/msz008) (LFMM2).
It does **no error handling** of any kind, nor any logging. It just automates
the procedures outlined in the manual with some degrees of freedom.
Please be careful when using it. It may kill your kittens and/or burn your
house down, but worst of all, it will tend to make you lazy regarding the inner
workings of *LFMM*.
The script takes no arguments, but all the variables you should need to edit
are presented at the start of the script, coupled with a short description.


### Variables

In the beginning of the script there are several line with empty variables. You should fill in the correct values for your case in order to use it.
Each option is pretty much self documented with both an explanation of what is expected and an example.


### Functions

The rest of the script is comprised of functions to wrap the LFMM functionality. You should not have to change anything here to get a complete run as all the parameters that are likely be changed can be from the "Variables" section.


### Running the script

Simply enter values for the variables at the beginning of the script and run it in R. It is recommended that you take a close look at the "PCA variance explained" plot to choose the best "K" to use, which will likely mean performing multiple runs.
