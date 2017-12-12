# pyRona
A python implementation of "Risk of non Adaptedness" method (with a bit of R too!)

[![Build Status](https://travis-ci.org/StuntsPT/pyRona.svg?branch=master)](https://travis-ci.org/StuntsPT/pyRona)

## Baypass_workflow.R

This script will automate the workflow for
the awesome [BayPass](http://www1.montpellier.inra.fr/CBGP/software/baypass/)
software by M. Gautier, which is described in
[this paper](http://www.genetics.org/content/early/2015/10/20/genetics.115.181453).
It does **no error handling** of any kind, nor any logging. It just automates
the procedures outlined in the manual with some degrees of freedom.
Please be careful when using it. It may kill your kittens and/or burn your
house down, but worst of all, it will tend to make you lazy regarding the inner
workings of BayPass. Please give BayPass's manual (and paper) a through read
before using this script.
The script takes no arguments, but all the variables you should need to edit
are presented at the start of the script, coupled with a short description.


## pyRona.py

This script is not yet properly documented since the implementation is still subject to change.
For now, you can use the program with `-h` parameter for a list of options.
Sorry about that.

## License
*pyRona* is licensed under the [GNU GPLv3](https://www.gnu.org/licenses/gpl-3.0-standalone.html).
