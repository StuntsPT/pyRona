#!/usr/bin/python3
# Copyright 2017 Francisco Pina Martins <f.pinamartins@gmail.com>
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


import argparse as ap


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

    parameters.add_argument("-covars", dest="num_covars", type=int,
                            default=3, required=False,
                            help="Number of covars to calculate the RONA for.")

    parameters.add_argument("-outliers", dest="outliers", type=int, default=2,
                            required=False, choices=[0, 1, 2],
                            help="Number of outliers to remove. 0 does no "
                                 "outier removal, 1 removes **at most** 1 "
                                 "outlier and 2 removes **any** number of "
                                 "outliers that match the distance criteria.")

    parameters.add_argument("-immutables", dest="immutables",
                            default=["1", "2", "3"], required=False, nargs="+",
                            help="List of immutable covariates. These are "
                                 "not even parsed from the betai file. By "
                                 "default the first 3 covars are skipped. "
                                 "You can enter any other values here.")

    parameters.add_argument("-ronatype", dest="rtype", type=str,
                            default="absdiff", required=False,
                            choices=["diff", "absdiff", "dist"],
                            help="Type of RONA to calculate. Default is "
                                 "absolute difference as in Rellstab et al. "
                                 "2016. Other options are 'difference' (not "
                                 "abs) and 'distance' (future vs. trendline).")

    io_opts.add_argument("-pc", dest="present_covars_file", type=str,
                         required=True, help="File with Present environmental "
                                             "data.")

    io_opts.add_argument("-fc", dest="future_covars_file", type=str,
                         required=True, help="File with Future environmental "
                                             "data.")

    io_opts.add_argument("-pop", dest="popnames_file", type=str,
                         required=True, help="File with population names.")

    io_opts.add_argument("-beta", dest="baypass_summary_betai_file", type=str,
                         required=True, help="Baypass summary betai file.")

    io_opts.add_argument("-pij", dest="baypass_pij_file", type=str,
                         required=True, help="Baypass pij file.")

    io_opts.add_argument("-out", dest="outfile", type=str,
                         required=True, help="Path to where RONA plot should "
                                             "be saved. Supports PDF, SVG and "
                                             "PNG extensions.")

    misc_opts.add_argument("-no-plots", dest="plots", action='store_false',
                           help="Pass this option if you don't want "
                                "individual regression plots to be drawn.",
                           required=False, default=True)

    misc_opts.add_argument("-no-weighted-means", dest="use_weights",
                           action='store_false',
                           help="Pass this option if you don't want to use"
                                "weighted means for RONA calculations.",
                           required=False, default=True)

    arguments = parser.parse_args(args)

    return arguments
