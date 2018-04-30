# *pyRona* changelog

## Changes made in version 0.3.4

* *Immutables* option is no longer on by default
* Improved help text of the *immutables* option
* Improved handling of LFMM population data (as opposed to individual data)
* Tests altered to conform to new options

---

## Changes made in version 0.3.3

* Better handling of individual plots
* Better PEP8 conformance

---

## Changes made in version 0.3.2

* Corrects a bug in when using outliers
* Better PEP8 conformance
* Implements automated test coverage measurements

---

## Changes made in version 0.3.1

* Corrects a bug in the installation of the R script `LFMM_workflow.R`
* Better PEP8 conformance

---

## Changes made in version 0.3.0

* `pyRona` can now use LFMM results as an input to calculate RONA
    * The test suite was also increased to match the changes

---

## Changes made in version 0.2.0

* Major code re-organization
* Can now be installed via `setup.py`
* [Pypi submission](https://pypi.python.org/pypi/pyRona/)
* Adds unit tests
* Travis-CI integration
* Proper documentation is started

---

## Changes made in version 0.1.3

* Autosave plots instead of interactive (which allows for a lot better automation)
* Bugfix with `arg.immutables` argument.

---

## Changes made in version 0.1.2

* Better PEP8 conformance
* Weighted R² means in the summary (when requested)
* Better label handling in the plots
* Option to choose which covariates are immutable and skip them

---

## Changes made in version 0.1.1

* Added a new line to the summary (Number of SNPs associated with each covariate)
* Can now parse multiple outputs from *BayPass* and not only for the MCMC_core model.

---

## Changes up to version 0.1.0

First fully working version.
There is still a lot to do to turn the program into a true "hit", but the basics are in.
