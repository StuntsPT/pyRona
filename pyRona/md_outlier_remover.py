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

import numpy as np


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
    for _, item in enumerate(diff_xy):
        mh_dist.append(np.sqrt(np.dot(np.dot(np.transpose(item),
                                             inv_covariance_xy), item)))

    return mh_dist


def md_remove_outliers(x_coords, y_coords):
    """
    Removes outliers based on Mahalanobis Distance.
    Takes 2 np.array([]) as input which are used to calculate the MD distance.
    Returns an np.array([]) with the indices of the removed outliers.
    http://kldavenport.com/mahalanobis-distance-and-outliers/
    """
    mahalanobis_dists = mahalanobis_dist_calculator(x_coords, y_coords)
    threshold = np.mean(mahalanobis_dists) * 1.5  # adjust 1.5 accordingly

    # Single or no outliers approach
    # if outliers == 1:
    if max(mahalanobis_dists) >= threshold:
        outlier_indeces = [mahalanobis_dists.index(max(mahalanobis_dists))]
    else:
        outlier_indeces = []

    # # Multiple outlier approach
    # elif outliers == 2:
    #     n_x, n_y, outlier_indeces = [], [], []
    #     for i, item in enumerate(mahalanobis_dists):
    #         if item <= threshold:
    #             n_x.append(x_coords[i])
    #             n_y.append(y_coords[i])
    #         else:
    #             outlier_indeces.append(i)  # position of removed pair

    return np.array(outlier_indeces)
