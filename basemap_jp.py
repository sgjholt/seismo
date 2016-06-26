
from obspy import read
import numpy as np
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt

# llcrnrlat,llcrnrlon,urcrnrlat,urcrnrlon
# are the lat/lon values of the lower left and upper right corners
# of the map.
# lat_ts is the latitude of true scale.
# resolution = 'c' means use crude resolution coastlines.

path_to_kikdb = (
'/Volumes/J_Holt_HDD/MRes/Modules/Thesis/Data/sitepub_kik_en_full.csv')
path_to_jpdb = (
'/Volumes/J_Holt_HDD/MRes/Modules/Thesis/Data/quake_lat_lon.csv')

dat = np.genfromtxt(path_to_kikdb, dtype = str, delimiter=',')
jp = np.loadtxt(path_to_jpdb, dtype=float, delimiter=',')
plt.figure(figsize=(14,10))
m = Basemap(projection='merc',llcrnrlat=30,urcrnrlat=46,\
            llcrnrlon=129,urcrnrlon=146,resolution='h')
m.drawcoastlines()
#m.fillcontinents(color='coral',lake_color='aqua')
m.etopo()
lats = np.array(dat[:,1], dtype=float)
lons = np.array(dat[:,2], dtype=float)
names = dat[:,0]

eqlats = jp[:,0]
eqlons = jp[:,1]

x, y = m(lons, lats)
eqx, eqy = m(eqlons, eqlats)

m.scatter(x, y, 20, color="r", marker="v", edgecolor="k", zorder=3)
m.scatter(eqx, eqy, 250, color="y", marker="*", edgecolor="k", zorder=4)
#for i in range(int(len(names)/6)):
    #plt.text(x[i*6], y[i*6], names[i*6], va="top", family="monospace")
# draw parallels and meridians.
m.drawparallels(np.arange(30.,50.,5.), labels=[1,0,0,0])
m.drawmeridians(np.arange(130.,150.,5.), labels=[0,0,0,1])
#m.drawmapboundary(fill_color='aqua')

plt.title("Kik-Net Stations and Earthquake Locations", fontsize=16, y = 1.01)

plt.savefig('/Users/jamesholt/Dropbox/MRes/Modules/Thesis/poster/jap_map.pdf')
plt.show()
