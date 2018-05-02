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

ax = plt.axes(projection=ccrs.PlateCarree())
ax.stock_img()
ax.coastlines()

fig = plt.figure(figsize=(12, 12), facecolor="none")
ax  = plt.axes(projection=ccrs.Mercator())

    # make the map global rather than have it zoom in to
    # the extents of any plotted data

ax.set_global()
ax.coastlines()
ax.stock_img()

fig.savefig("/tmp/aa.png")
