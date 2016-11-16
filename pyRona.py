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


from collections import defaultdict
from sys import argv

import argparse as ap
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


def calculate_rona(marker_name, covar_name, present_covar, future_covar,
                   allele_freqs, popnames_file, draw_plot=True):
    """
    Calculates the "Risk of non adaptation" (RONA) of each popuation for a
    given association.
    Also plots the associations if requested.
    """
    # Get population names
    popnames = open(popnames_file, 'r').readlines()

    # Remove outliers
    outlier_pos = md_remove_outliers(present_covar, allele_freqs)
    present_covar = np.delete(present_covar, outlier_pos)
    future_covar = np.delete(future_covar, outlier_pos)
    allele_freqs = np.delete(allele_freqs, outlier_pos)
    popnames = np.delete(popnames, outlier_pos)

    # Calculate trendline:
    fit = np.polyfit(present_covar, allele_freqs, 1)
    fit_fn = np.poly1d(fit)

    print("Marker: %s; covar %s" % (marker_name, covar_name))
    for pres, fut, freq, pops in zip(present_covar, future_covar, allele_freqs,
                                     popnames):
        pres_trendline_value = fit_fn(pres)
        fut_trendline_value = fit_fn(fut)

        pres_distance = freq - pres_trendline_value
        fut_distance = freq - fut_trendline_value
        diff_distance = abs(pres_distance) - abs(fut_distance)
        rel_distance = diff_distance / max(allele_freqs)

        print("%s: %s" % (pops.strip(), rel_distance))

    if draw_plot is True:
        all_covars = np.append(present_covar, future_covar)

        # Set-up the plot
        plt.xlabel("Covariate %s" % covar_name)
        plt.ylabel('Marker %s standardized allele freqs.' % marker_name)
        plt.title('Linear regression plot')

        plt.plot(present_covar, allele_freqs, 'bo')
        plt.plot(future_covar, allele_freqs, 'go')
        plt.plot(all_covars, fit_fn(all_covars), 'r--')

        plt.xlim(min(all_covars) - np.average(present_covar) * 0.1,
                 max(all_covars) + np.average(present_covar) * 0.1)
        plt.ylim(0, max(allele_freqs) + 2)

        # Annotation
        for label, x, y in zip(popnames, present_covar, allele_freqs):
            plt.annotate(label.strip(), xy=(x, y), xytext=(-15, 15),
                         textcoords='offset points', ha='right',
                         va='bottom', bbox=dict(boxstyle='round,pad=0.1',
                                                fc='yellow',
                                                alpha=0.5),
                         arrowprops=dict(arrowstyle='->',
                                         connectionstyle='arc3,rad=0'))

        plt.show()


def mahalanobis_dist_calculator(x_coords, y_coords):
    """
    Calculates Mahalanobis Distance.
    Takes 2 np.array([]) as input which are used to calculate the MD distance.
    Returns a list with the MD distances between every xy point.
    http://kldavenport.com/mahalanobis-distance-and-outliers/
    """
    covariance_xy = np.cov(x_coords, y_coords, rowvar=0)
    inv_covariance_xy = np.linalg.inv(covariance_xy)
    xy_mean = np.mean(x_coords), np.mean(y_coords)
    x_diff = np.array([x_i - xy_mean[0] for x_i in x_coords])
    y_diff = np.array([y_i - xy_mean[1] for y_i in y_coords])
    diff_xy = np.transpose([x_diff, y_diff])

    mh_dist = []
    for i in range(len(diff_xy)):
        mh_dist.append(np.sqrt(np.dot(np.dot(np.transpose(diff_xy[i]),
                                             inv_covariance_xy), diff_xy[i])))

    return mh_dist


def md_remove_outliers(x_coords, y_coords):
    """
    Removes outliers based on Mahalanobis Distance.
    Takes 2 np.array([]) as input which are used to calculate the MD distance.
    Returns an np.array([]) with the indices of the removed outliers.
    http://kldavenport.com/mahalanobis-distance-and-outliers/
    """
    mahalanobis_dist = mahalanobis_dist_calculator(x_coords, y_coords)
    threshold = np.mean(mahalanobis_dist) * 1.5 # adjust 1.5 accordingly
    n_x, n_y, outliers = [], [], []
    for i in range(len(mahalanobis_dist)):
        if mahalanobis_dist[i] <= threshold:
            n_x.append(x_coords[i])
            n_y.append(y_coords[i])
        else:
            outliers.append(i) # position of removed pair

    return np.array(outliers)


def argument_parser(args):
    """
    Parses arguments and returns them in a neat variable.
    """

    # Argument list
    parser = ap.ArgumentParser(description="A program to calculate "
                                           "the 'Risk of Non Adaptation' "
                                           "(RONA), based on BayPass "
                                           "output.",
                               prog="pyRona",
                               formatter_class=ap.RawTextHelpFormatter)

    io_opts = parser.add_argument_group("Input/Output options")
    parameters = parser.add_argument_group("Program execution options")
    misc_opts = parser.add_argument_group("Miscellaneous options")

    parameters.add_argument("-bf", dest="bayes_factor", type=float,
                            default=20, required=True,
                            help="Bayes factor treshold for considering "
                                 "associations.")

    io_opts.add_argument("-pc", dest="present_covar_file", type=str,
                         required=True, help="File with Present environmental "
                                             "data.")

    io_opts.add_argument("-fc", dest="future_covar_file", type=str,
                         required=True, help="File with Future environmental "
                                             "data.")

    io_opts.add_argument("-pop", dest="popnames_file", type=str,
                         required=True, help="File with population names.")

    io_opts.add_argument("-be", dest="baypass_summary_beta2_file", type=str,
                         required=True, help="Baypass summary beta2 file.")

    io_opts.add_argument("-pij", dest="baypass_pij_file", type=str,
                         required=True, help="Baypass pij file.")

    misc_opts.add_argument("-plot", dest="plots", type=bool, required=False,
                           help="Set this option to 'False' if you don't want "
                                "plots to be drawn.",
                           default=True)

    misc_opts.add_argument("-purge-outliers", dest="outlers", type=bool,
                           help="Set this option to 'False' if you don't want "
                                "to purge outliers from the data.",
                           required=False, default=True)

    arguments = parser.parse_args(args)


    return arguments


def main(params):
    """
    Main function. Takes all the inputs as arguments and runs the remaining
    functions of the program.
    """
    arg = argument_parser(params)
    present_covariates = parse_envfile(arg.present_covars_file)
    future_covariates = parse_envfile(arg.future_covars_file)
    assocs = baypass_summary_beta2_parser(arg.baypass_summary_beta2_file,
                                          arg.bayes_factor)
    al_freqs = baypass_pij_parser(arg.baypass_pij_file, assocs)
    # Get first 3 assocs, for testing
    for assoc in assocs[:3]:
        marker, covar = assoc
        calculate_rona(marker, covar, present_covariates[int(covar) + 1],
                       future_covariates[int(covar) + 1], al_freqs[marker],
                       arg.popnames_file)


if __name__ == "__main__":
    main(argv[1:])
