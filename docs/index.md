# *pyRona*

## Description

A python implementation of "Risk of non Adaptedness" method (with a bit of R too!)
*pyRona* calculates the "Risk of non Adaptedness" (RONA) of "demes" based on allelic frequencies, current environmental variables and projected future environmental variables. It is the first public implementation of the method developed in [Rellstab et al. 2016](doi.wiley.com/10.1111/mec.13889). You can find more detailed information in the [method description](description.md) section.


## Requirements

 - Python 3; *pyRona* has not been tested with python 2
 - numpy
 - matplotlib


## Where to get it

* Source code - [pyRona on github](https://github.com/StuntsPT/pyRona)
* Source distribution on pypi - [Sturcture_threader on Pypi](https://pypi.python.org/pypi/pyRona/)
    * You can easily install *pyRona* by issuing the command `pip3 install pyRona`


## Contents

* [Installation & dependencies](install.md)
* [Method description](description.md)
* [Usage](usage.md)
  * [Wrapping BayPass](baypass.md)
  * [Wrapping LFMM](lfmm.md)
* [Output](output.md)
* [Citation](citation.md)
* [Future Plans](future.md)
* [FAQ](faq.md)


## Bug reporting

Found a bug or would like a feature added? Or maybe drop some feedback?
Just [open a new issue](https://github.com/StuntsPT/pyRona/issues/new).


## License

*pyRona* is licensed under the [GNU GPLv3](https://www.gnu.org/licenses/gpl-3.0-standalone.html).
