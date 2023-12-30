""" Script to test / demonstrate the effect of doing statistics with sun-synchronous satellites 
    
    I'm calculating Poynting flux from AMPS and SWIPE models 
    along DMSP F13 orbit and averaging the result. The result is 
    highly asymmetric between hemispheres, even though AMPS and SWIPE 
    have been shown to be largely symmetric.
    
    The asymmetries are due to biased sampling that is persistent due
    to the sun-synchronous orbit. Satellites that drift in local time,
    like Swarm, would give better results with this approach if enough
    data was used (Swarm takes 5 years to collect enough data to have
    global coverage in all seasons)
"""


import tle
from ppigrf.ppigrf import geoc2geod # github.com/klaundal/ppigrf
import polplot # github.com/klaundal/polplot
import pyamps # github.com/klaundal/pyamps
import dipole #github.com/klaundal/dipole
from pyamps.mlt_utils import mlon_to_mlt 
from pyswipe.swipe import get_E # github.com/Dartspacephysiker/pyswipe
import apexpy
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

TIMERES = 10 # seconds between each sample
DLAT = 2 # latitude resolution in grid used for statistics
MU0 = np.pi * 4e-7

# make a set of sample times at which to  calcualte orbit
dates = pd.date_range('2000-1-1 00:00', '2000-12-31 23:59', freq = str(TIMERES) + 's')
#dates = pd.date_range('2000-1-1 00:00', '2000-1-31 23:59', freq = str(TIMERES) + 's')

# calculate orbit geocentric coordinates (TLE file taken from https://celestrak.org/NORAD/archives/)
r, theta, phi, v = tle.get_satellite_state(dates, 'dmspf13.txt', 'DMSPF13')
print('calculated orbit parameters')

# convert to geodetic coordinates (take into account geoid)
lat, height, v_north, v_down = geoc2geod(theta, r, -v[:, 1], v[:, 2])

# calcualte magnetic coordinates
mlat, mlon = apexpy.Apex(dates[0].year, refh = 110).geo2apex(lat, phi, height)
mlt = mlon_to_mlt(mlon, dates, dates[0].year)

# filter data to polar regions
iii = (mlat > 50) | (mlat < -50)
v = v[iii]
mlat = mlat[iii]
theta = theta[iii]
lat = lat[iii]
height = height[iii]
mlon = mlon[iii]
mlt = mlt[iii]
phi = phi[iii]
dates = dates[iii]


v_h  = v[:, :-1].T
v_h  = v_h / np.linalg.norm(v_h, axis = 0) # east, north-component of unit vector in satellite direction
v_ct = np.vstack((v_h[1], -v_h[0])) # east, north-components of unit vector in cross-track direction

# get dipole tilt angle:
d = dipole.Dipole(epoch = dates[0].year)
tilt = d.tilt(dates) # dipole tilt angle for each timestamp

# get solar wind parameters
omnidata = pd.read_csv('omni_data.csv', parse_dates = True, index_col = 0)
omnidata = omnidata.reindex(dates).interpolate(method = 'linear', limit = 60//TIMERES)

# calculate AMPS magnetic
print('calculating AMPS values...')
F107 = np.full_like(tilt, 100) # keep this constant for simplicity
Be, Bn, Bu = pyamps.get_B_space(lat, phi, height, dates, np.abs(omnidata['Vx'].values), omnidata['BY_GSM'].values, omnidata['BZ_GSM'].values, tilt, F107, epoch = dates[0].year)
Bx = v_h[ 0] * Be + v_h[ 1] * Bn # along-track
By = v_ct[0] * Be + v_ct[1] * Bn # cross-track
print('AMPS done')

print('calculating SWIPE values...')
Ee, En, Eu = get_E(lat, phi, height, dates, np.abs(omnidata['Vx'].values), omnidata['BY_GSM'].values, omnidata['BZ_GSM'].values, tilt, F107, epoch = dates[0].year, coords = 'geo')
Ex = v_h[ 0] * Ee + v_h[ 1] * En # along-track
Ey = v_ct[0] * Ee + v_ct[1] * En # cross-track
print('SWIPE done')



# Calculate Poynting flux using along-track and cross-track components
PF = pd.Series((Ex * By - Ey * Bx) * 1e-9 * 1e-3 / MU0, index = dates) * 1e3 # mW/m^2

# now we grid the data
sdgrid, mltres = polplot.grids.sdarngrid(dlat = DLAT, latmin = 50) # an almost equal-area grid used in superdarn community
grid_num_mag_north = polplot.grids.bin_number(sdgrid,  mlat, mlt)
grid_num_geo_north = polplot.grids.bin_number(sdgrid,  lat , (phi / 15) % 24)
grid_num_mag_south = polplot.grids.bin_number(sdgrid, -mlat, mlt)
grid_num_geo_south = polplot.grids.bin_number(sdgrid, -lat , (phi / 15) % 24)

# split the arrays into individual polar passes (this return slice objects):
north_pass_slices = np.ma.clump_unmasked(np.ma.masked_less_equal   (mlat,  60))
south_pass_slices = np.ma.clump_unmasked(np.ma.masked_greater_equal(mlat, -60))

# PLOT STATISTICS
fig, axes = plt.subplots(ncols = 2, nrows = 2, figsize = (14, 14))
levels = np.linspace(-0.25, 1.25, 22)

# northern hemisphere, magnetic:
avg_PF = PF.groupby(grid_num_mag_north).median().reindex(np.arange(len(mltres)))
pax = polplot.Polarplot(axes[0, 0])
print(avg_PF.min(), avg_PF.max())
pax.filled_cells(sdgrid[0], sdgrid[1], DLAT, mltres, avg_PF.values, levels = levels, cmap = plt.cm.bwr)

# southern hemisphere, magnetic:
avg_PF = PF.groupby(grid_num_mag_south).median().reindex(np.arange(len(mltres)))
pax = polplot.Polarplot(axes[0, 1])
print(avg_PF.min(), avg_PF.max())
pax.filled_cells(sdgrid[0], sdgrid[1], DLAT, mltres, avg_PF.values, levels = levels, cmap = plt.cm.bwr)

# northern hemisphere, geographic:
avg_PF = PF.groupby(grid_num_geo_north).median().reindex(np.arange(len(mltres)))
pax = polplot.Polarplot(axes[1, 0])
print(avg_PF.min(), avg_PF.max())
pax.filled_cells(sdgrid[0], sdgrid[1], DLAT, mltres, avg_PF.values, levels = levels, cmap = plt.cm.bwr)
pax.coastlines(color = 'grey', linewidth = 0.5, resolution = '110m')

# southern hemisphere, geographic:
avg_PF = PF.groupby(grid_num_geo_south).median().reindex(np.arange(len(mltres)))
pax = polplot.Polarplot(axes[1, 1])
print(avg_PF.min(), avg_PF.max())
pax.filled_cells(sdgrid[0], sdgrid[1], DLAT, mltres, avg_PF.values, levels = levels, cmap = plt.cm.bwr)
pax.coastlines(color = 'grey', linewidth = 0.5, resolution = '110m', north = False)

plt.subplots_adjust(left = 0.01, wspace= 0.01)

fig.savefig('PF_maps.png', dpi = 200)


# PLOT ORBITS IN MAGNETIC COORDS
fig, axes = plt.subplots(ncols = 2)

paxes = [polplot.Polarplot(x, minlat = 60, linestyle = ':', color = 'black', linewidth = .5) for x in  axes]
paxes = np.array(paxes)

# plot a random selection of orbits
N = 500 # number of lines to plot

# NORTH
iii = np.random.choice(np.arange(len(north_pass_slices)), N)
xs, ys = [], []
for i in range(N):
    x, y = paxes[0]._latlt2xy(mlat[north_pass_slices[iii[i]]], mlt[north_pass_slices[iii[i]]])
    xs.append(x)
    ys.append(y)

max_l = np.max([len(x) for x in xs])
xarr = np.full((len(xs), max_l), np.nan) # shape: number of lines x max_l
yarr = np.full((len(xs), max_l), np.nan) # shape: number of lines x max_l
for i in range(N):
    xarr[i][:len(xs[i])] = xs[i]
    yarr[i][:len(ys[i])] = ys[i]

segs = np.dstack((xarr, yarr))
paxes[0].ax.add_collection(LineCollection(segs, linewidths=0.3, colors='lightgrey', linestyle='solid'))

# SOUTH
iii = np.random.choice(np.arange(len(south_pass_slices)), N)
xs, ys = [], []
for i in range(N):
    x, y = paxes[1]._latlt2xy(mlat[south_pass_slices[iii[i]]], mlt[south_pass_slices[iii[i]]])
    xs.append(x)
    ys.append(y)

max_l = np.max([len(x) for x in xs])
xarr = np.full((len(xs), max_l), np.nan) # shape: number of lines x max_l
yarr = np.full((len(xs), max_l), np.nan) # shape: number of lines x max_l
for i in range(N):
    xarr[i][:len(xs[i])] = xs[i]
    yarr[i][:len(ys[i])] = ys[i]

segs = np.dstack((xarr, yarr))
paxes[1].ax.add_collection(LineCollection(segs, linewidths=0.3, colors='lightgrey', linestyle='solid'))

plt.subplots_adjust(left = 0.01, wspace= 0.01)

fig.savefig('bias_demo.png', dpi = 200)
plt.show()
