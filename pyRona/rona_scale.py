#!/usr/bin/python3

# Copyright 2018 Francisco Pina Martins <f.pinamartins@gmail.com>
# This file is part of pyRona.
# pyRona is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# pyRona is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with pyRona. If not, see <http://www.gnu.org/licenses/>.


from sys import argv
import matplotlib.pyplot as plt
import matplotlib.colors


def input_parser(values_filename):
    """
    Temporary function to parse a csv with the RONA value for each location.
    Takes a .csv file as input and retunrs a dict with {location:avg_RONA, ...}
    """
    ronas = {}
    infile = open(values_filename, 'r')
    infile.readline()  # Skip header
    for lines in infile:
        lines = lines.strip().split(",")
        ronas[lines[0]] = float(lines[1])

    return ronas


def values_to_color_scale(ronas):
    """
    Reads the RONA values and plots each location in a different color based on
    it.
    """
    locs = list(ronas.keys())
    vals = list(ronas.values())

    cmap = plt.cm.Wistia
    norm = matplotlib.colors.Normalize(vmin=min(vals), vmax=max(vals))

    fig, axe = plt.subplots()
    axe.bar(range(len(locs)), vals, color=cmap(norm(vals)))
    axe.set_xticks(range(len(locs)))
    axe.set_xticklabels(locs, rotation=90)

    scm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    scm.set_array([])
    fig.colorbar(scm)

    plt.show()


if __name__ == "__main__":
    RONAS = input_parser(argv[1])
    values_to_color_scale(RONAS)
