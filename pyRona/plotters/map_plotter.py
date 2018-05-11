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

try:
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
except ImportError:
    print("Error importing 'Cartopy'. Please look at pyRona's manual for "
          "information on how to install it on your system. Map plotting was "
          "not performed.")
    quit()
import matplotlib.pyplot as plt
import numpy as np


def map_plotter(ronas, latitudes, longitudes, out_filename):
    """
    Plots a map with each sampling site coloured per its averag RONA value
    """
    def _define_map_edges(latitudes, longitudes, padding=0.10):
        """
        Defines the edges of the map to be drawn.
        Takes a list of latitudes and longitudes as input and returns the map
         edges
        """
        # Define padding for the map edges
        hpad = padding
        vpad = padding

        # Get map edges
        max_lon = np.max(longitudes)
        max_lon = max_lon + abs(max_lon * vpad)
        min_lon = np.min(longitudes)
        min_lon = min_lon - abs(min_lon * vpad)

        max_lat = np.max(latitudes)
        max_lat = max_lat + abs(max_lat * hpad)
        min_lat = np.min(latitudes)
        min_lat = min_lat - abs(min_lat * hpad)

        return([min_lat, max_lat, min_lon, max_lon])

    max_ronas = []
    pop_names = ronas[0].pop_names
    for i, _ in enumerate(pop_names):
        max_ronas += [max([x.avg_ronas[i] for x in ronas])]

    fig = plt.figure(figsize=(22, 12), facecolor="none")
    map_area = plt.axes(projection=ccrs.Robinson())

    map_edges = _define_map_edges(latitudes, longitudes)

    map_area.set_extent(map_edges)
    map_area.coastlines(resolution='50m')
    cfeature.BORDERS.scale = "50m"
    map_area.add_feature(cfeature.BORDERS)

    dotplot = plt.scatter(latitudes, longitudes, c=max_ronas, s=700, vmin=0,
                          vmax=max(max_ronas), transform=ccrs.PlateCarree(),
                          cmap='autumn_r', zorder=2)

    for i, txt in enumerate(pop_names):
        map_area.annotate(txt, (latitudes[i], longitudes[i]))

    sidebar = fig.colorbar(dotplot)
    sidebar.ax.tick_params(labelsize=20)
    sidebar.set_label(label='RONA', size=30, weight='bold')

    fig.savefig(out_filename)

# TODO: Eventually make an interpolation
# TODO: Annotate the dots with the names (https://stackoverflow.com/questions/14432557/matplotlib-scatter-plot-with-different-text-at-each-data-point#14434334)
