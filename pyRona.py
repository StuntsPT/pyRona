#!/usr/bin/python3
# Copyright 2016 Francisco Pina Martins <f.pinamartins@gmail.com>
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
from sys import argv

import argparse as ap
import md_outiler_remover as mor
import numpy as np
import general_plotter as gp


class RonaClass:
    """
    Stores the RONA values for each covar
    """
    def __init__(self, covar):
        self.name = covar
        self.pop_names = []
        self.pop_ronas = defaultdict(list)
        self.corr_coef = {}
        self.avg_ronas = []
        self.stderr_ronas = []

    def basic_stats(self):
        """
        Gets the average RONA and stdev per population for each associated
        covariate. Stores the values in variables inside the class instance.
        """
        list_of_marker_values = []
        if len(self.pop_ronas) > 1:
            for marker_value in self.pop_ronas.values():
                list_of_marker_values.append(marker_value)

            list_of_marker_values = np.array(list_of_marker_values, dtype=float)
            for i in np.nditer(list_of_marker_values, flags=["external_loop"],
                               order="F"):
                self.avg_ronas += [np.average(i)]
                self.stderr_ronas += [np.std(i)/np.sqrt(len(i))]
        else:
            self.avg_ronas = [x for x in self.pop_ronas.values()][0]
            self.stderr_ronas = [0.0] * len(list(self.pop_ronas.values())[0])

    def count_markers(self):
        """
        Counts the number of markers in the instance.
        """
        return len(self.pop_ronas)


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
            frequencies[splitline[1]].append(float(splitline[4]))

    pij.close()

    # Make sure we have np.arrays for the freqencies
    for key in frequencies.keys():
        frequencies[key] = np.array(frequencies[key])

    return frequencies


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


def calculate_rona(marker_name, rona, present_covar, future_covar,
                   allele_freqs, plot, outliers):
    """
    Calculates the "Risk of non adaptation" (RONA) of each popuation for a
    given association.
    Also plots the associations if requested.
    """

    ronas = []

    # Remove outliers
    if outliers != 0:
        outlier_pos = mor.md_remove_outliers(present_covar, allele_freqs,
                                             outliers)
        present_covar = np.delete(present_covar, outlier_pos)
        future_covar = np.delete(future_covar, outlier_pos)
        allele_freqs = np.delete(allele_freqs, outlier_pos)
        rona.pop_names = np.delete(rona.pop_names, outlier_pos)

    # Calculate trendline:
    fit = np.polyfit(present_covar, allele_freqs, 1)
    fit_fn = np.poly1d(fit)

    # Get R²:
    rona.corr_coef[marker_name] = np.corrcoef(present_covar,
                                              allele_freqs)[1, 0] ** 2

    for pres, fut, freq in zip(present_covar, future_covar, allele_freqs):

        pres_trendline_value = fit_fn(pres)
        fut_trendline_value = fit_fn(fut)

        pres_distance = freq - pres_trendline_value
        fut_distance = freq - fut_trendline_value
        distance_diff = abs(pres_distance) - abs(fut_distance)
        amplitude = max(allele_freqs) - min(allele_freqs)
        rel_distance = distance_diff / amplitude

        rona.pop_ronas[marker_name] += [rel_distance]


    if plot is True:
        gp.draw_individual_plots(present_covar, future_covar, rona, marker_name,
                                 allele_freqs, fit_fn)

    header = "%s\t%s" % (marker_name, rona.corr_coef[marker_name])
    return [header] + ronas


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

    parameters.add_argument("-outliers", dest="outliers", type=int, default=2,
                            required=False, choices=[0, 1, 2],
                            help="Number of outliers to remove. 0 does no "
                                 "outier removal, 1 removes **at most** 1 "
                                 "outlier and 2 removes **any** number of "
                                 "outliers that match the distance criteria.")

    io_opts.add_argument("-pc", dest="present_covars_file", type=str,
                         required=True, help="File with Present environmental "
                                             "data.")

    io_opts.add_argument("-fc", dest="future_covars_file", type=str,
                         required=True, help="File with Future environmental "
                                             "data.")

    io_opts.add_argument("-pop", dest="popnames_file", type=str,
                         required=True, help="File with population names.")

    io_opts.add_argument("-beta", dest="baypass_summary_beta2_file", type=str,
                         required=True, help="Baypass summary beta2 file.")

    io_opts.add_argument("-pij", dest="baypass_pij_file", type=str,
                         required=True, help="Baypass pij file.")

    misc_opts.add_argument("-no-plots", dest="plots", action='store_false',
                           help="Pass this option if you don't want "
                                "plots to be drawn.",
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

    ronas = {}
    for assoc in assocs:
        marker, covar = assoc

        # Instanciate class
        if covar not in ronas:
            rona = RonaClass(covar)
            rona.pop_names = popnames_parser(arg.popnames_file)
        else:
            rona = ronas[covar]

        calculate_rona(marker, rona, present_covariates[int(covar) - 1],
                       future_covariates[int(covar) - 1], al_freqs[marker],
                       arg.plots, arg.outliers)

        ronas[covar] = rona

    for k, rona in ronas.items():
        rona.basic_stats()
        if rona.name == "8":
        # print(k)
        # print(rona.avg_ronas)
        # print(rona.stderr_ronas)
            gp.draw_rona_plot(rona, 1)

if __name__ == "__main__":
    main(argv[1:])
