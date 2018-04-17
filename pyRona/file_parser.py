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
from collections import OrderedDict
import numpy as np


# Baypass exclusive functions
def parse_baypass_envfile(envfile_filename):
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
    The order must be the same as in the baypass input files (geno and env).
    """
    popnames = []

    popfile = open(popnames_file, 'r')
    for lines in popfile:
        popnames.append(lines.strip())

    popfile.close()

    return popnames


def baypass_summary_betai_parser(summary_filename, bf_treshold, immutables):
    """
    Parses a baypass summary file to extract any significant associations
    between a marker and a covariate.
    Returns a list of associations tuples:
    [('marker', 'covariate'), ('marker', 'covariate')...]
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
def lfmm_results_parser(lfmm_results_filename, assoc_threshold, immutables):
    """
    Parses a lfmm results file to extract any significant associations
    between a marker and a covariate.
    Returns a list of associations tuples:
    [('marker', 'covariate'), ('marker', 'covariate')...]
    """
    associations = []
    snp = 0
    imut_indeces = [int(x) for x in immutables]
    results = open(lfmm_results_filename, "r")
    for lines in results:
        snp += 1
        lines = [float(x) for x in lines.split(",")]
        try:
            snp_assocs = [(str(snp), str(i)) for i in line_index
                          if lines[i] < assoc_threshold]
        except NameError:
            line_index = list(range(len(lines)))
            for index in sorted(imut_indeces, reverse=True):
                del line_index[index]
            snp_assocs = [(str(snp), str(i)) for i in line_index
                          if lines[i] < assoc_threshold]

        associations += snp_assocs

    results.close()

    return associations


def lfmm_to_pop_allele_freqs(lfmm_filename, env_filename, associations,
                             popnames=False):
    """
    Parses a LFMM input file and an evironemtal variables file to group
    the allele frequencies into "populations".
    The enviromental file should have the population names in the first
    column.
    Returns a dict with markers as keys and a np.array of allelic frequencies
    as values:
    {marker:np.array([freq_pop1, freq_pop2, ...])}
    """
    def _process_alleles(allele_matrix, indices):
        """
        Turns allele counts int allele frequencies.
        """
        frequencies = []
        allele_matrix = list(allele_matrix)

        absolutes = []
        startindex = 0
        for i in indices:
            alleles = allele_matrix[startindex:i + 1]
            absolutes.append(alleles)
            al_count = len(alleles)
            total = (al_count - alleles.count(9)) * 2
            try:
                freq = sum([0 if x == 9 else x for x in alleles]) / total
            except ZeroDivisionError:
                freq = np.nan

            frequencies.append(freq)
            startindex = i + 1

        return frequencies

    env_vars = open(env_filename, "r")
    pops = []
    indices = []
    for lines in env_vars:
        popname = lines.split()[0]
        try:
            if popname != pops[-1]:
                indices.append(len(pops) - 1)
        except IndexError:
            pass
        pops.append(popname)
    indices.append(len(pops) - 1)

    env_vars.close()

    collapsed_pops = list(OrderedDict.fromkeys(pops))

    id_freqs = {}
    lfmm = np.genfromtxt(lfmm_filename, delimiter=" ", dtype=float)
    snp_num = 0

    markers, _ = zip(*associations)
    markers = set(markers)

    for snp in lfmm.T:
        snp_num += 1
        if str(snp_num) in markers:
            snp_data = _process_alleles(snp, indices)
            if snp_data is None:
                id_freqs[str(snp_num)] = None
            else:
                id_freqs[str(snp_num)] = np.array(snp_data)

    if popnames:
        return collapsed_pops, id_freqs
    else:
        return id_freqs


def parse_lfmm_envfile(envfile_filename):
    """
    Parses an ENVFILE (lfmm) input file and returns a list of np arrays
    (one per line). Covariate values are averaged per "population".
    """
    envfile = open(envfile_filename, 'r')

    env_data = {}
    for lines in envfile:
        lines = lines.strip().split()
        if lines[0] not in env_data:
            env_data[lines[0]] = np.array([float(x) for x in lines[1:]])
        else:
            env_data[lines[0]] = np.vstack((env_data[lines[0]],
                                            np.array([float(x) for x in
                                                      lines[1:]])))

    envfile.close()

    covariates = []
    for pop in env_data.values():
        covariates.append(np.array([np.average(x) for x in pop.T]))
    covariates = np.array(covariates)
    return list(covariates.T)
