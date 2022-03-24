""" global plots (by M. Madelaire) """

#%% Import
import numpy as np
import scipy
from scipy import interpolate
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import cartopy.crs as ccrs

#%%
def global_scatterplot(lat, lon, x, title='', subtitle='', cmap='PuOr_r', clabel='', vmin=None, vmax=None, coast=True, MLT=False):
    if subtitle == '':
        subtitle = 'Data-points: {}'.format(int(len(x)))
    
    fig = plt.figure(figsize = (12, 12))
    gs=GridSpec(7,6)    
    
    # North Pole
    proj = ccrs.NorthPolarStereo()
    ax=fig.add_subplot(gs[0:3,0:3], projection=proj)
    if coast:
        ax.coastlines(resolution='50m')
    ax.gridlines()
    ax.set_extent([-180, 180, 50, 90], ccrs.PlateCarree())
    theta = np.linspace(0, 2*np.pi, 100)
    center, radius = [0.5, 0.5], 0.5
    verts = np.vstack([np.sin(theta), np.cos(theta)]).T
    circle = matplotlib.path.Path(verts * radius + center)
    ax.set_boundary(circle, transform=ax.transAxes)
    ax.scatter(lon, lat, 5, x, 'o', cmap=cmap, vmin=vmin, vmax=vmax, transform=ccrs.PlateCarree())
    ax.text(-4*10**6, -5*10**6, 'North Pole')
    if MLT:
        ax.text(-1*10**6, 5*10**6, '12 MLT')
    
    # South Pole
    proj = ccrs.SouthPolarStereo()
    ax=fig.add_subplot(gs[0:3, 3:6], projection=proj)
    if coast:
        ax.coastlines(resolution='50m')
    ax.gridlines()
    ax.set_extent([-180, 180, -90, -50], ccrs.PlateCarree())
    ax.set_boundary(circle, transform=ax.transAxes)
    ax.scatter(-1*((lon+180)%360)+360, lat, 5, x, 'o', cmap=cmap, vmin=vmin, vmax=vmax, transform=ccrs.PlateCarree())
    ax.text(1*10**6, -5*10**6, 'South Pole')
    if MLT:
        ax.text(-1*10**6, 5*10**6, '12 MLT')

    # Hammer
    proj = ccrs.Mollweide()
    ax=fig.add_subplot(gs[3::,:], projection=proj)
    if coast:
        ax.coastlines(resolution='50m')
    ax.gridlines()
    cs = ax.scatter(lon, lat, 5, x, 'o', cmap=cmap, vmin=vmin, vmax=vmax, transform=ccrs.PlateCarree())
    ax.text(-4*10**6, 0.95*10**7, subtitle)
    if not (title == ''):
        ax.text(-3*10**6, 1.05*10**7, title)
    if MLT:
        ax.text(-1*10**6, 0, '0 MLT')
    plt.tight_layout()

    # Colorbar
    cax,kw = matplotlib.colorbar.make_axes(ax,location='bottom',pad=0.05,shrink=0.5)
    out=fig.colorbar(cs,cax=cax,extend='both',**kw)
    if not( clabel == ''):
        out.set_label(clabel, size=12)
    return fig

#%%
def global_contour(lat, lon, x, res=1, n_levels=30, title='', subtitle='', cmap='PuOr_r', clabel='', vmin=None, vmax=None, coast=True, MLT=False, log10=False):    
    # Create default subtitle
    if subtitle == '':
        subtitle = 'Grid: {} degrees'.format(res)
    
    # Define amount of levels in contour
    levels = np.linspace(vmin, vmax, n_levels+1)

    # Start figure
    fig = plt.figure(figsize = (12, 12))
    gs=GridSpec(7,6)    

    # North Pole
        # Create meshgrid 
    lat_grid = np.arange(40, 90+res, res)
    lon_grid = np.arange(0, 360+res, res)
    lat_grid, lon_grid = np.meshgrid(lat_grid, lon_grid)
        # Interpolate to meshgrid
    var = interpolate.griddata((lat, lon), x, (lat_grid, lon_grid))
    if log10:
        flag_neg = var < 0
        x = np.log10(abs(var))
        for j in range(0, var.shape[1]):
            var[flag_neg[:, j], j] *= -1
        # Remove white patches
    var[var<levels[0]] = levels[0]
    var[var>levels[-1]] = levels[-1]
    
        # Plot
    proj = ccrs.NorthPolarStereo()
    ax=fig.add_subplot(gs[0:3,0:3], projection=proj)
    ax.set_global()
    if coast:
        ax.coastlines(resolution='50m')
    ax.gridlines()
    ax.set_extent([-180, 180, 50, 90], ccrs.PlateCarree())
    theta = np.linspace(0, 2*np.pi, 100)
    center, radius = [0.5, 0.5], 0.5
    verts = np.vstack([np.sin(theta), np.cos(theta)]).T
    circle = matplotlib.path.Path(verts * radius + center)
    ax.set_boundary(circle, transform=ax.transAxes)
    cs = ax.contourf(lon_grid, lat_grid, var, cmap=cmap, transform=ccrs.PlateCarree(), levels=levels)
    ax.text(-4*10**6, -5*10**6, 'North Pole')
    if MLT:
        ax.text(-1*10**6, 5*10**6, '12 MLT')
    
    # South Pole
        # Create meshgrid 
    lat_grid = np.arange(-90, -40+res, res)
    lon_grid = np.arange(0, 360+res, res)
    lat_grid, lon_grid = np.meshgrid(lat_grid, lon_grid)
        # Interpolate to meshgrid
    var = interpolate.griddata((lat, -1*((lon+180)%360)+360), x, (lat_grid, lon_grid))
    #lon_grid = -1*((lon_grid+180)%360)+360
    if log10:
        flag_neg = var < 0
        var = np.log10(abs(var))
        for j in range(0, var.shape[1]):
            var[flag_neg[:, j], j] *= -1
            
        # Remove white patches
    var[var<levels[0]] = levels[0]
    var[var>levels[-1]] = levels[-1]
        # Plot
    proj = ccrs.SouthPolarStereo()
    ax=fig.add_subplot(gs[0:3, 3:6], projection=proj)
    ax.set_global()
    if coast:
        ax.coastlines(resolution='50m')
    ax.gridlines()
    ax.set_extent([-180, 180, -90, -50], ccrs.PlateCarree())
    ax.set_boundary(circle, transform=ax.transAxes)
    cs = ax.contourf(lon_grid, lat_grid, var, cmap=cmap, transform=ccrs.PlateCarree(), levels=levels)
    ax.text(1*10**6, -5*10**6, 'South Pole')
    if MLT:
        ax.text(-1*10**6, 5*10**6, '12 MLT')

    # Hammer
        # Create meshgrid
    lat_grid = np.arange(-90, 90+res, res)
    lon_grid = np.arange(0, 360+res, res)
    lat_grid, lon_grid = np.meshgrid(lat_grid, lon_grid)
        # Interpolate to meshgrid
    var = interpolate.griddata((lat, lon), x, (lat_grid, lon_grid))
    if log10:
        flag_neg = x < 0
        var = np.log10(abs(var))
        for j in range(0, var.shape[1]):
            var[flag_neg[:, j], j] *= -1

        # Remove white patches
    var[var<levels[0]] = levels[0]
    var[var>levels[-1]] = levels[-1]
        # Plot
    proj = ccrs.Mollweide()
    ax=fig.add_subplot(gs[3::,:], projection=proj)
    ax.set_global()
    if coast:
        ax.coastlines(resolution='50m')
    ax.gridlines()
    cs = ax.contourf(lon_grid, lat_grid, var, cmap=cmap, transform=ccrs.PlateCarree(), levels=levels)
    ax.text(-4*10**6, 0.95*10**7, subtitle)
    if not (title == ''):
        ax.text(-3*10**6, 1.05*10**7, title)
    if MLT:
        ax.text(-1*10**6, 0, '0 MLT')
    #fig.canvas.draw()
    plt.tight_layout()

    # Colorbar
    cax,kw = matplotlib.colorbar.make_axes(ax,location='bottom',pad=0.05,shrink=0.5)
    out=fig.colorbar(cs,cax=cax,extend='both',**kw)
    if not( clabel == ''):
        out.set_label(clabel, size=12)
    return fig