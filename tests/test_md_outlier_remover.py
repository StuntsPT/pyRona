#!/usr/bin/python3
# Copyright 2016-2017 Francisco Pina Martins <f.pinamartins@gmail.com>
# and Joao Baptista <baptista.joao33@gmail.com>
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

import numpy as np
import pyRona.md_outlier_remover as mor

import pytest
assert pytest  # Hacky solution to stop the linters from complaining


def test_mhbs_dist_calculator():
    """
    Test the function mahalanobis_dist_calculator of md_outiler_remover.py.
    """

    x_coord = np.array([114, 71, 59, 117, 107, 86, 93, 87, 80, 76, 71, 122,
                        118, 88, 98, 83])
    y_coord = np.array([0.1957648, -0.49158111, 0.3691684, 0.27138479,
                        -0.00115721, 0.29038402, -0.31815359, 0.23805522,
                        0.15150419, -0.31698508, 0.4346372, 0.13334157,
                        -0.20163875, 0.8946341, -0.38784942, -0.55615943])

    test_dist = mor.mahalanobis_dist_calculator(x_coord, y_coord)

    control_dist = [1.2267382045449271, 1.758397269372304, 1.8952957763364617,
                    1.445191223272221, 0.79804140504730503,
                    0.69154415059165508, 0.919720382796942,
                    0.55031781984857764, 0.67436391809087115,
                    1.2473515530849255, 1.462153503620737, 1.5970060007168538,
                    1.4930142356624654, 2.1635163089688803, 1.1362365779638146,
                    1.5995900954742315]

    assert test_dist == control_dist


def test_md_remove_outliers():
    """
    Test the function md_remove_outliers of md_outiler_remover.py.
    """

    x_coords = np.array([114, 71, 59, 117, 107, 86, 93, 87, 80, 76, 71, 122,
                         118, 88, 98, 83])
    y_coords = np.array([0.1957648, -0.49158111, 0.3691684, 0.27138479,
                         -0.00115721, 0.29038402, -0.31815359, 0.23805522,
                         0.15150419, -0.31698508, 0.4346372, 0.13334157,
                         -0.20163875, 0.8946341, -0.38784942, -0.55615943])

    test_removed = mor.md_remove_outliers(x_coords, y_coords)

    control_removed = np.array([13])

    assert test_removed == control_removed
