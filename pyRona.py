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

import numpy as np
import matplotlib.pyplot as plt

from collections import defaultdict


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
        marker = splitline[1]
        covariate = splitline[0]
        bf_value = float(splitline[4])
        if bf_value >= float(bf_treshold):
            associations.append((marker, covariate))

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


def example():
    y = np.array([0.84097025,0.75629977,0.81300453,0.92482685,0.97089265,0.95719490,1.05375884,0.92784199,0.90493792,0.87466980,1.06920686,0.95997053,0.99620446,0.93575349,0.93254295,0.88172193]) # Marker 420

    y = np.array([6.233, 10.000, 10.033, 1.333, 11.300, 7.000, 4.867, 9.633, 7.400, 10.967, 10.433, 6.600, 6.233, 8.667, 11.633, 8.667]) # Marker 43

    x = np.array([6.233, 10.000, 10.033, 1.333, 11.300, 7.000, 4.867, 9.633, 7.400, 10.967, 10.433, 6.600, 6.233, 8.667, 11.633, 8.667]) # Covariate 8

    x = np.array([69, 175, 121, 60, 134, 110, 93, 133, 71, 144, 112, 111, 77, 121, 116, 126]) # Covariate 10

    fit = np.polyfit(x, y, 1)
    fit_fn = np.poly1d(fit)

    print(fit)

    # Confidence interval
    coord_y = [np.min(fit_fn(x)), np.max(fit_fn(x))]
    coord_x = [np.min(x), np.max(x)]

    predict_y = fit[0] * x + fit[1]
    print(predict_y)

    error_y = y - predict_y

    predict_x = np.arange(np.min(x),np.max(x)+1,1)

    mean_x = np.mean(x)
    n = len(x)
    t = 2.31 # Two tailed 95%
    sum_sq_err = np.sum(np.power(error_y, 2))

    confs = t * np.sqrt((sum_sq_err/(n - 2)) * (1.0 / n + (np.power((predict_x - mean_x), 2) / ((np.sum(np.power(x, 2))) - n * (np.power(mean_x, 2))))))

    # now predict y based on test x-values
    predict_y = fit[0] * predict_x + fit[1]

    # get lower and upper confidence limits based on predicted y and confidence intervals
    lower = predict_y - abs(confs)
    upper = predict_y + abs(confs)

    # /Confidence interval

    # set-up the plot
    #plt.axes().set_aspect('equal')
    plt.xlabel('Covariate 10')
    plt.ylabel('Marker 43 standard freqs')
    plt.title('Linear regression and confidence limits')

    plt.plot(x, y, 'bo', label='Observations')
    plt.plot(coord_x, coord_y, 'r-', label='Regression line')

    plt.plot(predict_x, lower, 'b--',label='Lower confidence limit (95%)')
    plt.plot(predict_x, upper, 'b--',label='Upper confidence limit (95%)')


    plt.xlim(0, 200)
    plt.ylim(0, 15)

    plt.show()

if __name__ == "__main__":
    from sys import argv
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
