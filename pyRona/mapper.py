#!/usr/bin/python3
# Copyright 2018 Francisco Pina Martins <f.pinamartins@gmail.com>
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

import cartopy.crs as ccrs
import matplotlib.pyplot as plt
from file_parser import parse_lfmm_envfile
import numpy as np


env_data = parse_lfmm_envfile("tests/data/LFMM_covars.txt")

# Define padding for the map edges
hpad = 0.10
vpad = 0.10

# Get map edges
max_lon = np.max(env_data[0])
max_lon = max_lon + abs(max_lon * vpad)
min_lon = np.min(env_data[0])
min_lon = min_lon - abs(min_lon * vpad)


max_lat = np.max(env_data[1])
max_lat = max_lat + abs(max_lat * hpad)
min_lat = np.min(env_data[1])
min_lat = min_lat - abs(min_lat * hpad)

ax = plt.axes(projection=ccrs.PlateCarree())
ax.stock_img()
ax.coastlines()

fig = plt.figure(figsize=(12, 12), facecolor="none")
ax  = plt.axes(projection=ccrs.Robinson())

    # make the map global rather than have it zoom in to
    # the extents of any plotted data

ax.set_extent([min_lat, max_lat, min_lon, max_lon])

ax.coastlines(resolution='50m')

ax.pcolormesh()

fig.savefig("/tmp/aa.png")
