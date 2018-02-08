#!/usr/bin/python3

# Copyright 2018 Francisco Pina Martins <f.pinamartins@gmail.com>
# This file is part of geste2lfmm.
# geste2lfmm is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# geste2lfmm is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with geste2lfmm. If not, see <http://www.gnu.org/licenses/>.

# Usage: python3 geste2lfmm.py file.geste file.lfmm

from collections import OrderedDict


def parse_geste(infile_name):
    """
    Parses a GESTE file and retuns an OrderedDict with:
    {"Population_name":[Freq_ref_allele_on SNP_1,Freq_ref_allele_on SNP_2,...]}
    """
    infile = open(infile_name, "r")
    pop_freqs = OrderedDict()
    pop_starter = "[pop]="
    popname = ""
    for line in infile:
        # Neat trick to ignore data that is not SNP info
        # This code section should be very performant since it replaces most
        # if - else tests with try -> except statements
        line = line.split()
        try:
            int(line[0])
        except ValueError:  # In case it's a new section
            if line[0].startswith(pop_starter):
                popname = "Pop %s" % line[0].strip().replace(pop_starter, "")
                pop_freqs[popname] = []
            continue
        except IndexError:  # In case it's an empty line
            continue
        try:
            ref_frequency = round(int(line[3]) / int(line[1]), 3)
        except ZeroDivisionError:
            ref_frequency = 9
        pop_freqs[popname].append(ref_frequency)

    infile.close()

    return pop_freqs


def write_lfmm(pop_freqs, lfmm_filename):
    """
    Write a LFMM inpt file based on the OrderedDict extracted from the GESTE
    file.
    """
    outfile = open(lfmm_filename, 'w')
    for freqs in pop_freqs.values():
        outfile.write("\t".join(map(str, freqs)) + "\n")


if __name__ == "__main__":
    from sys import argv
    POP_FREQS = parse_geste(argv[1])
    write_lfmm(POP_FREQS, argv[2])
