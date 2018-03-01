#!/usr/bin/python3
# Copyright 2016-2018 Francisco Pina Martins <f.pinamartins@gmail.com>
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


from collections import defaultdict
import numpy as np


# Common functions to Baypass & LFMM:
def parse_envfile(envfile_filename):
    """
    Parses an ENVFILE (baypass) input file and returns a list of np arrays
    (one per line)
    """
    envfile = open(envfile_filename, 'r')
    covariates = []
    for lines in envfile:
        covariates.append(np.array([float(x) for x in lines.split()]))

    envfile.close()

    return covariates


def popnames_parser(popnames_file):
    """
    Parses a file with population names and returns a list with these names.
    The order is the same as in the file.
    """
    popnames = []

    popfile = open(popnames_file, 'r')
    for lines in popfile:
        popnames.append(lines.strip())

    popfile.close()

    return popnames


# Baypass excusive functions
def baypass_summary_betai_parser(summary_filename, bf_treshold, immutables):
    """
    Parses a baypass summary file to extract any significant associations
    between a marker and a covariate.
    Returns a list of associations tuples:
    [(marker, covariate), (marker, covariate)...]
    """
    summary = open(summary_filename, 'r')
    header = summary.readline()  # Skip header and get BF column.
    for identifier in ("eBPmc", "BF(dB)"):
        try:
            bf_col = header.strip().split().index(identifier)
        except ValueError:
            pass
    associations = []
    for lines in summary:
        splitline = lines.strip().split()
        covariate = splitline[0]
        # Remove "immutable" covariates (Lat, Long, Alt)
        if covariate not in immutables:
            marker_id = splitline[1]
            bf_value = float(splitline[bf_col])
            if bf_value >= float(bf_treshold):
                associations.append((marker_id, covariate))

    summary.close()

    return associations


def baypass_pij_parser(pij_filename, associations):
    """
    Parses a baypass pij file to extract standardized allelic frequencies.
    Returns a dict with markers as keys and a np.array of allelic frequencies
    as values:
    {marker:np.array([freq_pop1, freq_pop2, ...])}
    """
    marker_list = [x for x, y in associations]
    frequencies = defaultdict(list)

    pij = open(pij_filename, 'r')
    for lines in pij:
        splitline = lines.strip().split()
        if splitline[1] in marker_list:
            frequencies[splitline[1]].append(float(splitline[4]))

    pij.close()

    # Make sure we have np.arrays for the freqencies
    for key in frequencies.keys():
        frequencies[key] = np.array(frequencies[key])

    return frequencies


# LFMM exclusive functions
def parse_allele_freqs(allele_freqs_file, associations):
    """
    Parses a file that contains allelic frequencies (in this case, the LFMM
    input file). Uses the associations list to filter which frequencies get
    returned.
    Returns a dict with markers as keys and a np.array of allelic frequencies
    as values:
    {marker:np.array([freq_pop1, freq_pop2, ...])}
    Only markers with significant associations are placed in the dict.
    """
    marker_list = [x for x, y in associations]
    freqs_array = np.genfromtxt(allele_freqs_file, delimiter="\t")
    marker = 0
    frequencies = {}
    for snp in freqs_array.T:
        marker += 1
        if str(marker) in marker_list:
            frequencies[str(marker)] = snp

    return frequencies


def lfmm_results_parser(lfmm_results_filename, assoc_threshold, immutables):
    """
    Parses a lfmm results file to extract any significant associations
    between a marker and a covariate.
    Returns a list of associations tuples:
    [(marker, covariate), (marker, covariate)...]
    """
    associations = []
    snp = 0
    imut_indeces = [int(x) for x in immutables]
    results = open(lfmm_results_filename, "r")
    for lines in results:
        snp += 1
        lines = [float(x) for x in lines.split(",")]
        try:
            snp_assocs = [(str(snp), i) for i in line_index
                          if lines[i] < assoc_threshold]
        except NameError:
            line_index = list(range(len(lines)))
            for index in sorted(imut_indeces, reverse=True):
                del line_index[index]
            snp_assocs = [(str(snp), i) for i in line_index
                          if lines[i] < assoc_threshold]

        associations += snp_assocs

    results.close()

    print(len(associations))
    return associations


if __name__ == "__main__":
    lfmm_results_parser("/home/francisco/pvalues.csv", 0.05, [1, 2, 3])
