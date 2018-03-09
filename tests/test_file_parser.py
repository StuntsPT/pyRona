#!/usr/bin/python3
# Copyright 2016-2018 Francisco Pina Martins <f.pinamartins@gmail.com>
# and Joao Baptista <baptista.joao33@gmail.com>
# This file is part of pyRona.
# pyRona is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# pyRona is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with pyRona.  If not, see <http://www.gnu.org/licenses/>.


import pickle
import pyRona.file_parser as fp

import pytest
assert pytest  # Hacky solution to stop the linters from complaining


def test_parse_baypass_envfile():
    """
    Test the function parse_baypass_envfile of file_parser.py.
    """
    test_parsed = fp.parse_baypass_envfile("../tests/data/ENVFILE")
    test_parsed_list = [list(x) for x in test_parsed]

    with open("../tests/data/jar/file_parser.parse_baypass_envfile.pickle",
              "rb") as fle:
        control_parsed = pickle.load(fle)

    control_parsed_list = [list(x) for x in control_parsed]

    assert test_parsed_list == control_parsed_list


def test_popnames_parser():
    """
    Test the function popnames_parser of file_parser.py.
    """
    test_popname = fp.popnames_parser("../tests/data/popnames_single_GEO.txt")

    with open("../tests/data/jar/"
              "file_parser.popnames_parser.pickle", "rb") as fle:
        control_popname = pickle.load(fle)

    assert test_popname == control_popname


def test_baypass_s_b_p():
    """
    Test the function baypass_summary_betai_parser of file_parser.py.
    """
    test_betai = fp.baypass_summary_betai_parser(
        "../tests/data/Qsuber_GBS_mcmc_aux_summary_betai.out", 20,
        ["1", "2", "3"])

    with open("../tests/data/jar/"
              "file_parser.baypass_summary_betai_parser.pickle", "rb") as fle:
        control_betai = pickle.load(fle)

    assert test_betai == control_betai


def test_baypass_pij_parser():
    """
    Test the function baypass_pij_parser of file_parser.py.
    """
    test_pij = fp.baypass_pij_parser(
        "../tests/data/Qsuber_GBS_mcmc_aux_summary_pij.out",
        fp.baypass_summary_betai_parser(
            "../tests/data/Qsuber_GBS_mcmc_aux_summary_betai.out",
            20, ["1", "2", "3"]))

    with open("../tests/data/jar/"
              "file_parser.baypass_pij_parser.pickle", "rb") as fle:
        control_pij = pickle.load(fle)

    assert str(test_pij) == str(control_pij)


def test_lfmm_res_parser():
    """
    Test the function lfmm_results_parser of file_parser.py.
    """
    lfmm_results = fp.lfmm_results_parser(
        "../tests/data/Qsuber_lfmm_results.csv",
        0.01,
        ["1", "2", "3"])
    with open("../tests/data/jar/file_parser.lfmm_results_parser.pickle",
              "rb") as fle:
        control_pvalues = pickle.load(fle)

    assert str(lfmm_results) == str(control_pvalues)
