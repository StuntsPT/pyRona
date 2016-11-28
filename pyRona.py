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
import file_parser as fp


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

    def basic_stats(self, use_weights):
        """
        Gets the average RONA and stdev per population for each associated
        covariate. Stores the values in variables inside the class instance.
        """
        if len(self.pop_ronas) > 1:
            # Sort markers:
            markers = sorted([x for x in self.corr_coef.keys()])

            list_of_marker_values = np.array([self.pop_ronas[x] for x in
                                              markers],
                                             dtype=float)
            corr_weights = np.array([self.corr_coef[x] for x in markers],
                                    dtype=float)

            for i in np.nditer(list_of_marker_values, flags=["external_loop"],
                               order="F"):
                if use_weights is True:
                    self.avg_ronas += [np.average(i, weights=corr_weights)]
                else:
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
        # TODO: Select plot type
        pres_trendline_value = fit_fn(pres)
        fut_trendline_value = fit_fn(fut)

        pres_distance = freq - pres_trendline_value
        fut_distance = freq - fut_trendline_value
        distance_diff = abs(pres_distance) - abs(fut_distance)

        # distance_diff = abs(fut_distance)

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

    misc_opts.add_argument("-no-weighted-means", dest="use_weights",
                           action='store_false',
                           help="Pass this option if you don't want to use"
                                "weighted means for RONA calculations.",
                           required=False, default=True)

    arguments = parser.parse_args(args)

    return arguments


def main(params):
    """
    Main function. Takes all the inputs as arguments and runs the remaining
    functions of the program.
    """
    arg = argument_parser(params)
    present_covariates = fp.parse_envfile(arg.present_covars_file)
    future_covariates = fp.parse_envfile(arg.future_covars_file)
    assocs = fp.baypass_summary_beta2_parser(arg.baypass_summary_beta2_file,
                                          arg.bayes_factor)
    al_freqs = fp.baypass_pij_parser(arg.baypass_pij_file, assocs)

    ronas = {}
    for assoc in assocs:
        marker, covar = assoc

        # Instanciate class
        if covar not in ronas:
            rona = RonaClass(covar)
            rona.pop_names = fp.popnames_parser(arg.popnames_file)
        else:
            rona = ronas[covar]

        calculate_rona(marker, rona, present_covariates[int(covar) - 1],
                       future_covariates[int(covar) - 1], al_freqs[marker],
                       arg.plots, arg.outliers)

        ronas[covar] = rona

    # Delete immutable covariates:
    del ronas["1"]
    del ronas["2"]
    del ronas["3"]

    sortable_representation = {}
    for k, rona in ronas.items():
        rona.basic_stats(arg.use_weights)
        sortable_representation[k] = len(rona.pop_ronas)

    top_represented = sorted(sortable_representation,
                             key=sortable_representation.get, reverse=True)[:3]
    top_ronas = [ronas[x] for x in top_represented]

    gp.draw_rona_plot(top_ronas)


if __name__ == "__main__":
    main(argv[1:])
