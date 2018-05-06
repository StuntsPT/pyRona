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
import cartopy.feature as cfeature
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

fig = plt.figure(figsize=(22, 12), facecolor="none")
ax = plt.axes(projection=ccrs.Robinson())

ax.set_extent([min_lat, max_lat, min_lon, max_lon])

cfeature.BORDERS.scale = "50m"

ax.coastlines(resolution='50m')
ax.add_feature(cfeature.BORDERS)
max_rona = 0.35
cbat = np.arange(max_rona, 0, -0.01).reshape(int(max_rona * 100), 1)
im = ax.imshow(cbat, cmap='gist_earth')

# Plot dots
ax.plot(8.8, 40.0, 'bo', markersize=22, transform=ccrs.PlateCarree())
ax.plot(-8.4, 40.0, 'bo', markersize=22, transform=ccrs.Geodetic())

plt.colorbar(im, ax=ax, label='RONA')

fig.savefig("/tmp/aa.png")

# TODO: Sample the dot colour from the colourbar
# TODO: Eventually make an interpolation
# TODO: Paint the sea
