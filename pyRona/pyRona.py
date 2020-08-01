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
from sys import argv

import numpy as np
try:
    import md_outlier_remover as mor
    import plotters.general_plotter as gp
    import file_parser as fp
    from argparser import argument_parser
except ImportError:
    import pyRona.md_outlier_remover as mor
    import pyRona.plotters.general_plotter as gp
    import pyRona.file_parser as fp
    from pyRona.argparser import argument_parser


class RonaClass:
    """
    Stores the RONA values for each covar
    """
    POP_NAMES = []

    def __init__(self, covar):
        self.name = covar
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
            markers = sorted(list(self.corr_coef.keys()))

            list_of_marker_values = np.array([self.pop_ronas[x] for x in
                                              markers],
                                             dtype=float)
            corr_weights = np.array([self.corr_coef[x] for x in markers],
                                    dtype=float)

            for i in list_of_marker_values.T:
                not_nans = ~np.isnan(i)
                if use_weights is True:
                    if True in not_nans:
                        self.avg_ronas += [np.average(i[not_nans],
                                                      weights=corr_weights
                                                      [not_nans])]
                    else:
                        self.avg_ronas += [np.nan]
                else:
                    self.avg_ronas += [np.average(i[not_nans])]

                self.stderr_ronas += [np.std(i[not_nans])
                                      / np.sqrt(len(i[not_nans]))]
        else:
            self.avg_ronas = [x for x in self.pop_ronas.values()][0]
            self.stderr_ronas = [0.0] * len(list(self.pop_ronas.values())[0])

    def count_markers(self):
        """
        Counts the number of markers in the instance.
        """
        return len(self.pop_ronas)


def calculate_rona(marker_name, rona, present_covar, future_covar,
                   allele_freqs, plot, outliers, rtype):
    """
    Calculates the "Risk of non adaptation" (RONA) of each popuation for a
    given association.
    Also plots the associations if requested.
    """
    # Remove outliers
    if outliers is True:
        outlier_pos = mor.md_remove_outliers(present_covar, allele_freqs)
        outlier_pos = np.array(outlier_pos, dtype='int')
        for i in outlier_pos:
            present_covar[i] = np.nan
            future_covar[i] = np.nan
            allele_freqs[i] = np.nan

        rona.pop_names = np.delete(RonaClass.POP_NAMES, outlier_pos)
    else:
        rona.pop_names = RonaClass.POP_NAMES

    # Calculate trendline:
    not_nan = ~np.isnan(present_covar)
    fit = np.polyfit(present_covar[not_nan], allele_freqs[not_nan], 1)
    fit_fn = np.poly1d(fit)

    # Get R²:
    rona.corr_coef[marker_name] = np.corrcoef(present_covar[not_nan],
                                              allele_freqs[not_nan])[1, 0] ** 2

    for pres, fut, freq in zip(present_covar, future_covar, allele_freqs):

        pres_distance = freq - fit_fn(pres)
        fut_distance = freq - fit_fn(fut)
        distance_diff = abs(pres_distance) - abs(fut_distance)

        amplitude = max(allele_freqs) - min(allele_freqs)

        if rtype == "diff":
            rel_distance = distance_diff / amplitude

        elif rtype == "absdiff":
            rel_distance = abs(distance_diff) / amplitude

        elif rtype == "dist":
            rel_distance = abs(fut_distance)

        rona.pop_ronas[marker_name] += [rel_distance]

    if plot is not None:
        gp.draw_individual_plots(present_covar, future_covar, rona,
                                 marker_name, allele_freqs, fit_fn, plot)


def results_summary(ronas, use_weights):
    """
    This function outputs a summary of the RONAS for each population and
    covariate.
    """
    pop_names = ronas[0].pop_names
    for i, j in enumerate(pop_names):
        if i == 0:
            print("Covar\t%s" % "\t".join([x.name for x in ronas]))
            print("#SNPs\t%s" % "\t".join([str(x.count_markers()) for x in
                                           ronas]))
        print("%s\t%s" % (j, "\t".join([str(x.avg_ronas[i]) for x in ronas])))

    print("Min R^2\t%s" %
          "\t".join([str(np.nanmin(list(x.corr_coef.values()))) for x in
                     ronas]))

    print("Max R^2\t%s" %
          "\t".join([str(np.nanmax(list(x.corr_coef.values()))) for x in
                     ronas]))

    # if use_weights is True:
    #     means = [str(np.average(list(x.corr_coef.values()),
    #                             weights=list(x.corr_coef.values()))) for x in
    #              ronas]
    # else:
    means = [str(np.nanmean(list(x.corr_coef.values()))) for x in ronas]
    print("Average R^2\t%s" % "\t".join(means))


def ronas_filterer(ronas, use_weights, num_covars):
    """
    Filters RONAS to remove immutable covars, and return only the top "n" most
    represented covariables.
    """
    sortable_representation = {}
    for k, rona in ronas.items():
        rona.basic_stats(use_weights)
        sortable_representation[k] = len(rona.pop_ronas)

    top_represented = sorted(sortable_representation,
                             key=sortable_representation.get,
                             reverse=True)[:num_covars]
    top_ronas = [ronas[x] for x in top_represented]

    return top_ronas


def main():
    """
    Main function. Takes all the inputs as arguments and runs the remaining
    functions of the program.
    """
    if len(argv) < 2:
        arg_list = ["-h"]
    else:
        arg_list = argv[1:]

    arg = argument_parser(arg_list)

    if arg.upstream == "baypass":
        present_covariates = fp.parse_baypass_envfile(arg.present_covars_file)
        future_covariates = fp.parse_baypass_envfile(arg.future_covars_file)
        RonaClass.POP_NAMES = fp.popnames_parser(arg.popnames_file)
        assocs = fp.baypass_summary_betai_parser(
            arg.baypass_summary_betai_file,
            arg.bayes_factor, arg.immutables)
        al_freqs = fp.baypass_pij_parser(arg.baypass_pij_file, assocs)
    elif arg.upstream == "lfmm":
        present_covariates = fp.parse_lfmm_envfile(arg.present_covars_file)
        future_covariates = fp.parse_lfmm_envfile(arg.future_covars_file)
        assocs = fp.lfmm_results_parser(arg.lfmm_assoc_file,
                                        arg.p_thres,
                                        arg.immutables)
        RonaClass.POP_NAMES, al_freqs = fp.lfmm_to_pop_allele_freqs(
            arg.allele_freqs_file,
            arg.present_covars_file,
            assocs,
            popnames=True)

    ronas = {}
    for assoc in assocs:
        marker, covar = assoc

        # Instanciate class
        if covar not in ronas:
            rona = RonaClass(covar)
        else:
            rona = ronas[covar]

        calculate_rona(marker, rona, present_covariates[int(covar) - 1],
                       future_covariates[int(covar) - 1],
                       al_freqs[marker],
                       arg.plots, arg.outliers, arg.rtype)

        ronas[covar] = rona

    ronas = ronas_filterer(ronas, arg.use_weights, arg.num_covars)

    results_summary(ronas, arg.use_weights)
    gp.draw_rona_plot(ronas, arg.outfile)

    if arg.map_filename is not None:
        # The map plotting module is only imported if a map plot is requested.
        # This is to be able to keep 'cartopy' as an optional dependency.
        try:
            import plotters.map_plotter as mapper
        except ImportError:
            import pyRona.plotters.map_plotter as mapper
        mapper.map_plotter(ronas, present_covariates[1], present_covariates[0],
                           arg.map_filename)


if __name__ == "__main__":
    main()
