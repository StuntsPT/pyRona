#!/usr/bin/python3

# Copyright 2017-2020 Francisco Pina Martins <f.pinamartins@gmail.com>
# This file is part of pyRona.
# pyRona is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# pyRona is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with pyRona. If not, see <http://www.gnu.org/licenses/>.


import sys
from setuptools import setup


class NotSupportedException(BaseException):
    """
    Dummy class to work as an exception in case someone attempts to install
    pyRona using pip2.
    """
    pass


if sys.version_info.major < 3:
    raise NotSupportedException("Only Python 3.x Supported")


# Set some variables (PKGBUILD inspired)
VERSION = "0.4.2"
URL = "https://gitlab.com/StuntsPT/pyRona"


setup(
    name="pyRona",
    version=VERSION,
    packages=["pyRona",
              "pyRona.plotters"],
    install_requires=["numpy",
                      "matplotlib"],
    description=('A python implementation of "Risk of non Adaptedness"'
                 'method (with a bit of R too!)'),
    url=URL,
    download_url="{0}/archive/v{1}.tar.gz".format(URL, VERSION),
    author="Francisco Pina-Martins",
    author_email="f.pinamartins@gmail.com",
    license="GPL3",
    classifiers=["Intended Audience :: Science/Research",
                 "License :: OSI Approved :: GNU General Public License v3 ("
                 "GPLv3)",
                 "Natural Language :: English",
                 "Operating System :: POSIX :: Linux",
                 "Topic :: Scientific/Engineering :: Bio-Informatics",
                 "Programming Language :: Python :: 3 :: Only",
                 "Programming Language :: Python :: 3.4",
                 "Programming Language :: Python :: 3.5",
                 "Programming Language :: Python :: 3.6",
                 "Programming Language :: Python :: 3.7",
                 "Programming Language :: Python :: 3.8"],
    data_files=[("bin", ["pyRona/R/Baypass_workflow.R",
                         "pyRona/R/LFMM_workflow.R"])],
    entry_points={
        "console_scripts": [
            "pyRona = pyRona.pyRona:main",
        ]
    },
)
