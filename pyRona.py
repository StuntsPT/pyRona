#!/usr/bin/python3
# Copyright 2016 Francisco Pina Martins <f.pinamartins@gmail.com>
# This file is part of misc_plotters.
# misc_plotters is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# misc_plotters is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with misc_plotters.  If not, see <http://www.gnu.org/licenses/>.

# Usage: python3

from collections import defaultdict

import numpy as np
import matplotlib.pyplot as plt

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


def baypass_summary_beta2_parser(summary_filename, bf_treshold):
    """
    Parses a baypass summary file to extract any significant associations
    between a marker and a covariate.
    Returns a list of associations tuples:
    [(marker, covariate), (marker, covariate)...]
    """
    summary = open(summary_filename, 'r')
    summary.readline() # Skip header
    associations = []
    for lines in summary:
        splitline = lines.strip().split()
        marker_id = splitline[1]
        covariate = splitline[0]
        bf_value = float(splitline[4])
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
            frequencies[splitline[1]].append(float(splitline[0]))

    pij.close()

    # Make sure we have np.arrays for the freqencies
    for key in frequencies.keys():
        frequencies[key] = np.array(frequencies[key])

    return frequencies


def calculate_rona(marker, covar, present_covar, future_covar, allele_freqs,
                   popnames):
    """
    Calculates the "Risk of non adaptation" (RONA) of each popuation for a
    given association.
    """
    # TODO: Improve this function, get rid of redundant code and make the code more readable!

    # Calculate trendline:
    fit = np.polyfit(present_covar, allele_freqs, 1)
    fit_fn = np.poly1d(fit)

    print("Marker: %s; covar %s" % (marker, covar))
    for pres, fut, freq, pops in zip(present_covar, future_covar, allele_freqs,
                                     popnames):
        pres_trendline_value = fit_fn(pres)
        fut_trendline_value = fit_fn(fut)

        pres_dist = freq - pres_trendline_value
        fut_dist = freq - fut_trendline_value
        diff_dist = abs(pres_dist) - abs(fut_dist)
        rel_dist = diff_dist / max(allele_freqs)

        print("%s: %s" % (pops.strip(), rel_dist))


def plot_associations(marker, covar, present_covar, future_covar, allele_freqs):
    """
    Plots the relevant associtions, along with a trendline
    """
    all_covars = np.append(present_covar, future_covar)

    # Calculate trendline:
    fit = np.polyfit(present_covar, allele_freqs, 1)
    fit_fn = np.poly1d(fit)


    # set-up the plot
    plt.xlabel("Covariate %s" % covar)
    plt.ylabel('Marker %s standardized allele freqs.' % marker)
    plt.title('Linear regression plot')

    plt.plot(present_covar, allele_freqs, 'bo', label='Present observations')
    plt.plot(future_covar, allele_freqs, 'go', label='Future predictions')
    plt.plot(all_covars, fit_fn(all_covars), 'r--', label='Regression line')

    plt.xlim(min(all_covars) - np.average(present_covar) * 0.1,
             max(all_covars) + np.average(present_covar) * 0.1)
    plt.ylim(0, max(allele_freqs) + 2)

    plt.show()


if __name__ == "__main__":
    from sys import argv
    # TODO: Implement argparse!
    present_covariates = parse_envfile(argv[1]) # ENVFILE
    future_covariates = parse_envfile(argv[2]) # ENVFILE
    assocs = baypass_summary_beta2_parser(argv[3], argv[4]) # summary; BF
    al_freqs = baypass_pij_parser(argv[5], assocs) # summary_pij
    popnames = open(argv[6], 'r').readlines()
    # Get first 3 assocs, for testing
    for assoc in assocs[:3]:
        marker, covar = assoc
        calculate_rona(marker, covar, present_covariates[int(covar) + 1],
                       future_covariates[int(covar) + 1], al_freqs[marker],
                       popnames)
        plot_associations(marker, covar, present_covariates[int(covar) + 1],
                          future_covariates[int(covar) + 1], al_freqs[marker])
