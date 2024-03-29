# *pyRona* changelog

## Changes made in version 0.4.4

### Bug fixes
* Fixed a bug with individual plots when using LFMM (which would cause the wrong covariate to be shown in the plot). Another big thank you to Gabriele Nocchi for finding it!

### New feature
* You can now provide a file with covariate names (one per line), and all plots will use this information instead of simply numbering the covariates.
* If you do not provide this file, by default your covariates will be numbered, **starting with 1** (previous to version 0.4.4 they were 0 indexed).

---

## Changes made in version 0.4.3

### Bug fixes
* Fixed a bug with delimiters on LFMM input files (files using tab as a delimiter would cause python tracebacks). A big thank you to Gabriele Nocchi for finding it!

### Documentation Updates
* Corrected a switched "Latitude" and "Longitude" information in the manual. Once more, big thank you Gabriele Nocchi for spotting it and reporting!

---

## Changes made in version 0.4.2

### Updates
* Minor improvements to LFMM2_Workflow.R

---

## Changes made in version 0.4.1

### Updates
* Corrects a bug introduced by changes in numpy 1.19.0

---

## Changes made in version 0.4.0

### Updates
* Adds a new workflow script to be able to use lfmm2

---

## Changes made in version 0.3.8

### Updates
* Changes the Bioconductor install method to the current one

---

## Changes made in version 0.3.7

### New features
* Adds infrastructure to use Gitlab's CI server parallel to Travis CI.

---

## Changes made in version 0.3.6

* Fixed a broken help link

---

## Changes made in version 0.3.5

* Adds a new feature: RONA maps
    * To use it, just add the option `-map` with a path to where you want to save the map, such as `~/my_RONA_map.png`
    * In order to use this feature, the dependency [cartopy](http://scitools.org.uk/cartopy/) is required.
    * The [manual](http://pyrona.readthedocs.io/en/latest/install/) now contains instructions on installing `cartopy` depending on your OS
        * Instructions for OSX are still missing. Help, anyone?
    * For now, the options assumes that the first variable in the environmental file is "Latitude" and the second one is "Longitude"

---

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
