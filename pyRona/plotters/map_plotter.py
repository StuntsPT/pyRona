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
    from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
except ImportError:
    print("\nError importing 'Cartopy'. Please look at pyRona's manual"
          "(http://pyrona.readthedocs.io/en/latest/install/#installing-cartopy)"
          " for information on how to install it on your system. Map plotting"
          " was **not** performed.")
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
    map_area = plt.axes(projection=ccrs.PlateCarree())

    map_edges = _define_map_edges(latitudes, longitudes)

    map_area.set_extent(map_edges)
    map_area.coastlines(resolution='50m')
    cfeature.BORDERS.scale = "50m"
    map_area.add_feature(cfeature.BORDERS)

    # Draw sampling sites
    dotplot = map_area.scatter(latitudes, longitudes, c=max_ronas, s=700,
                               vmin=0, vmax=max(max_ronas),
                               transform=ccrs.PlateCarree(),
                               cmap='autumn_r', zorder=2)

    # Label the locations
    for label, x, y in zip(pop_names, latitudes, longitudes):
        map_area.annotate(label.strip().replace("_", " "), xy=(x, y),
                          xytext=(0, -28), textcoords='offset points',
                          ha='center', va='bottom', fontsize=13,
                          bbox=dict(boxstyle='round,pad=0.1',
                                    fc='lightgrey', edgecolor="none"))

    # Control x and y ticks
    gridlines = map_area.gridlines(draw_labels=True)
    gridlines.xlines = False
    gridlines.ylines = False
    gridlines.ylabels_right = False
    gridlines.xlabels_top = False
    gridlines.xformatter = LONGITUDE_FORMATTER
    gridlines.yformatter = LATITUDE_FORMATTER
    gridlines.xlabel_style = {'size': 22}
    gridlines.ylabel_style = {'size': 22}

    # Control x and y labels
    map_area.text(-0.10, 0.55, 'Latitude', va='bottom', ha='center',
                  rotation='vertical', rotation_mode='anchor',
                  transform=map_area.transAxes, fontsize=28)
    map_area.text(0.5, -0.12, 'Longitude', va='bottom', ha='center',
                  rotation='horizontal', rotation_mode='anchor',
                  transform=map_area.transAxes, fontsize=28)

    # Sidebar settings
    sidebar = fig.colorbar(dotplot)
    sidebar.ax.tick_params(labelsize=20)
    sidebar.set_label(label='RONA', size=30, weight='bold')

    # Save the map
    fig.savefig(out_filename)

# TODO: Eventually make an interpolation
