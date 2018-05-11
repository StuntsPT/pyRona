# Installation

Due to the way different Operating Systems handle dependencies specific instructions for the preferred installation method are provided for each of the most used OS's "GNU/Linux", "MacOS" and "Windows".

## Preferred method, by platform

### GNU/Linux

1. Install python 3. Although python 3 is already installed by default in most modern Linux distributions, sometimes it may not be available (you can type `python3 --version` from a terminal to see if python 3 is installed. If it is, you will have something similar to `Python 3.6.1` printed on your terminal, if it isn't you will either see a "command not found" error, or a helpful message on how to install python 3). In "Debian based" distributions (such as Ubuntu) you can install python 3 by opening a terminal an running the command `sudo apt-get install python3`. In other Linux distributions you can similarly use your package manger to install it. If you do not have administration privileges in your environment, please ask your sysadmin to install python 3 for you. This is the only step that has a hard requirement on administrative privileges.
2. Install `pip`. `pip` is a [package manager for python](https://en.wikipedia.org/wiki/Pip_(package_manager)). If `pip` is not already installed in  your system, you can follow the official instructions on how to get it [here](https://pip.pypa.io/en/stable/installing/). **Make sure you run get-pip.py using python 3 in order to be able to use *pyRona*.** Like this: `python3 get-pip.py`
3. Install *pyRona*. Now that you have python 3 and `pip` installed, installing *pyRona* is just one command away: `pip3 install pyRona --user`. The `--user` option installs the software to a local directory, ensuring you do not need administration privileges to perform the installation.
4. Using *pyRona*. Running the command from step 3 will install the program to `~/.local/bin`. You can either run it by calling it directly `~/.local/bin/pyRona` or by adding the location `~/.local/bin` to your shell `$PATH` ([here is a good guide on how to do it](https://unix.stackexchange.com/questions/26047/how-to-correctly-add-a-path-to-path)) and just calling `pyRona`.


### MacOS

1. Install python 3. MacOS comes with python 2.7 installed by default, (Mac OSX versions before "Sierra" does not have any version of python installed) but in order to run *pyRona* you will need python 3.4 or above. You can [follow this comprehensive guide to do it](http://docs.python-guide.org/en/latest/starting/install3/osx/).
2. Install `pip`. `pip` is a [package manager for python](https://en.wikipedia.org/wiki/Pip_(package_manager)). You can use the guide from step 1 to install it on your system.
3. Install *pyRona*. Now that you have python 3 and `pip` installed, installing *pyRona* is just one terminal command away: `pip3 install pyRona --user`. The `--user` option installs the software to a local directory, ensuring you do not need administration privileges to perform the installation.
4. Using *pyRona*. Running the command from step 3 will install the program to `~/Library/Python/[PYTHON_VERSION]/bin` (where `PYTHON_VERSION` is the version of python you used, eg. `3.6`). You can either run it by calling it directly `~/Library/Python/[PYTHON_VERSION]/bin/pyRona` or by adding the location `~/Library/Python/[PYTHON_VERSION]/bin` to your shell `$PATH` ([here is a good guide on how to do it](https://unix.stackexchange.com/questions/26047/how-to-correctly-add-a-path-to-path)) and just calling `pyRona`.


### Windows

1. Install python 3. No version of Windows comes with python installed by default, but you can install it from [here](https://www.python.org/downloads/). If you need help installing python 3 for windows, here is [the official guide](https://docs.python.org/3/using/windows.html).
2. Install `pip`.  `pip` is a [package manager for python](https://en.wikipedia.org/wiki/Pip_(package_manager)). Using the instructions from the guide from step 1 will install `pip` for you.
3. Install *pyRona*. Now that you have python 3 and `pip` installed, installing *pyRona* is just one terminal command away: `C:\Python3.6\python.exe -m pip install pyRona`. Don't forget to change the path "python3.6" to whatever version of python 3 you have installed.
4. Using *pyRona*. Running the command from step 3 will install the program to `C:\Python3.6\Scripts`. You can run it by calling it directly `C:\Python3.6\Scripts\pyRona.exe`.


## Important note

Note that *pyRona* will also automatically install two *R* scripts called `Baypass_workflow.R` and `LFMM_workflow.R` under `~/.local/bin` that can be used to automate the usage of the upstream software *BayPass* and *LFMM*, whose output can be used as input for *pyRona*. For more information on these scripts, please see the [baypass](baypass.md) and [LFMM](lfmm.md) sections.


## Alternative methods (AKA 'expert mode')
You can also run *pyRona* by simply cloning the repository (or
downloading one of the releases), and placing the contents of the directory
"pyRona", on any location on your `$PATH` environment variable var.

Another alternative, is the `setup.py` method. This can be used either by running `python3 setup.py install` (or even
better, `pip3 install .`) from the distribution's root directory (where
`setup.py` is located).


## Installing cartopy

*pyRona* requires the [cartopy](http://scitools.org.uk/cartopy/) library to draw RONA maps. It is optional and the module will only be required if you specifically ask *pyRona* to draw a map.
Below are instructions on how to get it to run on three different Operating Systems. After running them, *pyRona* should be able to draw RONA maps using the `-map` argument.

### GNU/Linux

Under GNU/Linux OSes, `cartopy` requires the `proj4` and `geos` libraries. On Ubuntu and Debian based OSes, you can install it using the following command:

```
sudo apt-get install libproj-dev libgeos-dev
```

Then you can proceed to install `cartopy` and its dependencies `scipy` and `cython`:

```
pip install cartopy scipy cython --user
```


### Windows

Under Windows OSes, you can install the unofficial `cartopy` binaries, which can be found [here](https://www.lfd.uci.edu/~gohlke/pythonlibs/#cartopy).
You should still install `scipy` and `cython` via `pip`:

```
C:\Python3.6\python.exe -m pip install scipy cython
```

### MacOS

Installing `cartopy` under OSX does not seem to be an easy task. I really would appreciate some help here. The `pip` parts should be similar to those of GNU/Linux, but getting `proj4` and `geos` libraries installed might prove trickier.
