# Wrapping BayPass

*pyRona* can take the output of *BayPass* as input to calculate the RONA values. Due to this, *pyRona* includes an `R` script named `Baypass_workflow.R` that can be used to automate and wrap *BayPass* analyses.
Not that it is highly recommended that you read both *BayPass*'s [manual](http://www1.montpellier.inra.fr/CBGP/software/baypass/files/BayPass_manual_2.1.pdf) (and [paper](http://www.genetics.org/content/201/4/1555)) a through read before using this script.


## Baypass_workflow.R

This script will automate the workflow for the awesome [BayPass](http://www1.montpellier.inra.fr/CBGP/software/baypass/)
software by M. Gautier, which is described in [this paper](http://www.genetics.org/content/early/2015/10/20/genetics.115.181453).
It does **no error handling** of any kind, nor any logging. It just automates
the procedures outlined in the manual with some degrees of freedom.
Please be careful when using it. It may kill your kittens and/or burn your
house down, but worst of all, it will tend to make you lazy regarding the inner
workings of *BayPass*.
The script takes no arguments, but all the variables you should need to edit
are presented at the start of the script, coupled with a short description.

### Variables

In the beginning of the script there are several line with empty variables. You should fill in the correct values for your case in order to use the script. Each option is pretty much self documented with both an explanation of what is expected and an example.
