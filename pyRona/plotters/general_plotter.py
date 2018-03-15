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

from os import path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


COLOR_LIST = ['red', 'green', 'blue', 'cyan', 'magenta', 'black', 'orange',
              'yellow']


def draw_individual_plots(present_covar, future_covar, rona, marker_name,
                          allele_freqs, fit_fn, outpath):
    """
    Draws an individual marker vs covar plot.
    """
    all_covars = np.append(present_covar, future_covar)

    fig, axes = plt.subplots()
    # Set-up the plot
    axes.set_xlabel("Covariate %s" % rona.name)
    axes.set_ylabel('Marker %s standardized allele freqs.' % marker_name)
    axes.set_title('Linear regression plot')

    axes.plot(present_covar, allele_freqs, 'bo')
    axes.plot(future_covar, allele_freqs, 'go')
    axes.plot(all_covars, fit_fn(all_covars), 'r--')

    axes.set_xlim(min(all_covars) - np.nanmean(present_covar) * 0.1,
                  max(all_covars) + np.nanmean(present_covar) * 0.1)
    axes.set_ylim(min(allele_freqs) - min(allele_freqs) * 0.1,
                  max(allele_freqs) + max(allele_freqs) * 0.1)

    # Annotation
    for label, x, y in zip(rona.pop_names, present_covar, allele_freqs):
        axes.annotate(label.strip(), xy=(x, y), xytext=(-9, 9),
                      textcoords='offset points', ha='right',
                      va='bottom', bbox=dict(boxstyle='round,pad=0.1',
                                             fc='yellow',
                                             alpha=0.3),
                      arrowprops=dict(arrowstyle='->',
                                      connectionstyle='arc3,rad=0'))

    # Set outpath
    filename = "Cov{}_Mrk{}.pdf".format(rona.name, marker_name)
    if path.isdir(outpath):
        full_outpath = path.join(outpath, filename)
    else:
        exit("ERROR: '-draw-ind-plots' option must point to a directory!")

    # Draw the plot
    fig.savefig(full_outpath)
    plt.close(fig)


def draw_rona_plot(ronas, outpath):
    """
    Draws a RONA plot of the Nth most represented covariates.
    Plots the RONA+/-Stderr for each of the populations.
    """
    fig, ax = plt.subplots()
    width = 1 / (len(ronas) + 1)

    counter = 0
    axes = []
    names = []
    for rona in ronas:
        ind = np.arange((len(rona.pop_names)))
        rects = ax.bar(ind + width * counter, rona.avg_ronas, width,
                       color=COLOR_LIST[counter], yerr=rona.stderr_ronas,
                       alpha=0.5, error_kw=dict(ecolor='gray'))
        counter += 1
        axes.append(rects)
        names.append(rona.name)

    # add some text for labels, title and axes ticks
    ax.set_ylabel('RONA')

    # Define title:
    if len(ronas) > 1:
        title_msg = ("Average RONA per population for the most represented %s "
                     "covariates" % len(ronas))
    else:
        title_msg = "RONA per population for the most represented covariate"

    ax.set_title(title_msg)
    ax.set_xticks(ind + width)
    ax.set_xticklabels(rona.pop_names, rotation=45, va="top", ha="right")

    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.9, box.height])

    ax.legend(axes, names, loc='center left', bbox_to_anchor=(1, 0.5))

    # plt.tight_layout()
    plt.gcf().subplots_adjust(bottom=0.18)
    plt.gcf().subplots_adjust(right=0.85)

    fig.savefig(outpath)

    plt.close(fig)
