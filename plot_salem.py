from sys import platform as sys_pf
import os
import geopandas as gpd
from geopandas import GeoDataFrame
import pandas as pd
import matplotlib
if sys_pf == 'darwin':
    matplotlib.use("TkAgg")

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.legend_handler import HandlerLine2D
import matplotlib.colors as colors
from matplotlib import ticker
from shapely import wkt
from shapely.geometry import Point
import rasterio
import rasterio.plot
import shapely
import xarray as xr
import rasterio
import salem
import numpy as np
import numpy.matlib
from matplotlib.mlab import griddata
import pickle
import scipy
import warnings
warnings.filterwarnings("ignore")

SCATTER = False
mask_path='mask.pkl'

# TODO: labels of lakes, criteria on/off

def create_mask(xi,yi,zi, outline):
    try:
        mask = pickle.load(open(mask_path,'rb'))
        if mask.shape != zi.shape:
            raise
        else:
            return mask
    except:
        mask = np.zeros_like(zi, dtype=bool)
        geom = outline.iloc[0].geometry

        xn = np.matlib.repmat(xi, 1, len(yi))[0]
        yn = np.repeat(yi, len(xi))

        zn = gpd.GeoDataFrame()
        zn.loc[:, 'xn'] = xn
        zn.loc[:, 'yn'] = yn

        zn.loc[:, 'p'] = zn.apply(lambda x: Point(x.xn, x.yn), axis=1)
        zn.loc[:, 'm'] = ~zn.p.apply(lambda x: geom.contains(x))

        mask = np.reshape(zn.m.values, zi.shape)
        pickle.dump(mask, open(mask_path,'wb'))
        return mask


def plot_data(fname_csv, fname_png, cb='jet', cmin=0, cmax=100, cb_extend='neither', cb_label='', cb_ticks=[], EPSG='epsg:3031', res=1000, scale=1, log=False, symlog=False):
    print(plot_name)
    df = pd.read_csv(fname_csv, header=None, names=['x','y','z'])
    geometry = [Point(xy) for xy in zip(df.x, df.y)]
    crs = {'init': EPSG} 
    data_df = GeoDataFrame(df, crs=crs, geometry=geometry)
    fname = fname_png
    data_df.z = data_df.z.divide(scale)
    da = salem.DataLevels(data_df.z)
    da.set_cmap(cb)
    da.set_vmin(cmin)
    da.set_vmax(cmax)
    da.set_extend(cb_extend)
    x, y = sm.grid.transform(data_df.x.values, data_df.y.values, crs=data_df.crs)
    if SCATTER:
        ax.scatter(x, y, color=da.to_rgb(), s=1, linewidths=1)
        
    else:
        outline = gpd.read_file('rignot+grl_1.shp', crs=data_df.crs)
        sm.set_shapefile('rignot+grl_1.shp', color='0.2', linewidth=1)
        x = data_df.x
        y = data_df.y

        xi = np.linspace(x.min(), x.max(), res)
        yi = np.linspace(y.min(), y.max(), res)

        zi = griddata(x, y, data_df.z.values, xi, yi, interp='linear')

        mask = create_mask(xi, yi, zi, outline)
        
        zi = np.ma.array(zi, mask=mask)
        xi,yi = sm.grid.transform(xi, yi, crs=data_df.crs)
        if log==True:
            #cmin = cmin + 0.0001
            zi[zi<1e-5] = 1e-5
            norm = colors.LogNorm(vmin=cmin, vmax=cmax, clip=True)
            step = 0.01
            a = np.log10(cmin)
            b = np.log10(cmax)
            lev_exp = np.arange(a, b+step, step)
            levs = np.power(10, lev_exp)
            levs = np.concatenate(([zi.min()],levs,[zi.max()]), axis=None)
            cs = ax.contourf(xi, yi, zi, levels=levs, norm=norm, cmap=cb, alpha=0.6, antialiased=True) 
            cbar = fig.colorbar(cs, ax=ax, label=cb_label, pad=0.018, extend=cb_extend, aspect=50, shrink=0.823, fraction=0.017, drawedges=False, ticks=cb_ticks)
            cbar.ax.set_yticklabels(cb_ticks)
            cbar.solids.set_edgecolor("face")
        elif symlog==True:
            norm = colors.SymLogNorm(linthresh=0.03, linscale=0.03, vmin=cmin, vmax=cmax)
            #a = np.ceil(np.log10(cmax))
            #b = np.ceil(np.log10(-cmin))
            a = np.log10(cmax)
            b = np.log10(-cmin)
            step = 0.01
            A = np.arange(-1,a+step,step)
            B = np.arange(b,-(1+step),-step)
            levs_A = np.power(10, A)
            levs_B = np.power(10, B)
            levs = np.concatenate((-levs_B,levs_A), axis=None)
            #lev_exp = np.arange(cmin,cmax,1)
            #levs = np.power(10, lev_exp)
            #levs=lev_exp
            cs = ax.contourf(xi, yi, zi, levels=levs, norm=norm, cmap=cb, alpha=0.6, antialiased=True) 
            cbar = fig.colorbar(cs, ax=ax, label=cb_label, pad=0.018, extend=cb_extend, aspect=50, shrink=0.823, fraction=0.017, drawedges=False, ticks=cb_ticks)
            cbar.ax.set_yticklabels(cb_ticks)
            cbar.solids.set_edgecolor("face")
        else:
            step=1
            levs = np.arange(cmin, cmax, step)
            levs = np.concatenate(([zi.min()],levs,[zi.max()]), axis=None)
            cs = ax.contourf(xi, yi, zi, levels=levs, cmap=cb, vmin=cmin, vmax=cmax, alpha=0.6, antialiased=True)
            da.append_colorbar(ax, label=cb_label, pad=0.1, size='1.5%') #ticks=cb_ticks


for i in range(5,6):#12
    fig, ax = plt.subplots(dpi=300, figsize=(8,6))

    # radar raster
    ds = salem.open_xr_dataset("amm1mos_400m_recisl.tif")
    sm = ds.salem.get_map(cmap='Greys', countries=False, factor=1)
    sm.set_data(ds['data'])

    # lakes
    sm.set_shapefile('lakes.shp', linewidth=2, label='lakes')

    # shearmargin
    sm.set_shapefile('shearmargin.shp', linewidth=3., linestyle='dashed')

    # sinks
    sm.set_shapefile('RECISL_sinks_10km_v7_contours.shp', color='w', linewidth=2)


    # model data plot
    if i==0:
        # alpha     
        plot_name = 'alpha.png'
        plot_data('refined_new2_alpha.csv',plot_name,'Blues', 0, 150, 'max', '$[(Pa\ a/m)^{(1/2)}]$', [0,50,100,150], 'epsg:3031', 1000)
        
    elif i==1:
        # k
        plot_name = 'k.png'
        plot_data('refined_new2_k.csv',plot_name,'Blues', 0, 150, 'max', '$[(s/m)^{(1/2)}]$', [0,50,100,150], 'epsg:3031', 1000)

    elif i==2:
        # N_eff
        plot_name = 'N_eff.png'
        plot_data('refined_new2_Neff.csv',plot_name,'RdYlBu_r', 0, 30, 'max', '$[MPa]$', [0,10,20,30], 'epsg:3031', 1000, scale=1000000)
    elif i==3:
        # tau_b
        plot_name = 'tau_b.png'
        plot_data('refined_new2_tau_b.csv',plot_name,'YlGnBu', 5, 200, 'max', '$[kPa]$', [5,10,25,50,100,200], 'epsg:3031', 1000, scale=1000, log=True)
    elif i==4:
        #tau_d
        plot_name = 'tau_d.png'
        plot_data('refined_new2_tau_d.csv',plot_name,'YlGnBu', 5, 200, 'max', '$[kPa]$', [5,10,25,50,100,200], 'epsg:3031', 1000, scale=1000, log=True)
    elif i==5:
        # tau_d-tau_b
        plot_name = 'tau_d-tau_b.png'
        plot_data('refined_new2_tau_d-tau_b.csv',plot_name,'RdBu_r', -100, 100, 'both', '$[kPa]$', [-100,-50,-25,0,25,50,100], 'epsg:3031', 1000, scale=1000)
    elif i==6:
        # v_obs
        plot_name = 'v_obs.png'
        #plot_data('refined_new2_vel_obs.csv',plot_name,'Spectral_r', 0, 500, 'max', '$[m/a]$', [0,100,200,300,400,500], 'epsg:3031', 1000, log=True)
        plot_data('refined_new2_vel_obs.csv',plot_name,cb='Spectral_r', cmin=1, cmax=1000, cb_extend='max', cb_label='$[m/a]$', cb_ticks=[0.1,1,5,10,25,50,100,250,500,800],EPSG='epsg:3031', res=1000, log=True)
    elif i==7:
        # v_sim
        plot_name = 'v_sim.png'
        #plot_data('refined_new2_vel_sim_surf.csv',plot_name,'Spectral_r', 0, 500, 'max', '$[m/a]$', [0,100,200,300,400,500], 'epsg:3031', 1000)
        plot_data('refined_new2_vel_sim_surf.csv',plot_name,cb='Spectral_r', cmin=1, cmax=1000, cb_extend='max', cb_label='$[m/a]$', cb_ticks=[0.1,1,5,10,25,50,100,250,500,800],EPSG='epsg:3031', res=1000, log=True)
    elif i==8:
        # v_sim-v_obs
        plot_name = 'v_sim-v_obs.png'
        plot_data('refined_new2_v_sim-v_obs.csv',plot_name,'RdBu_r', -250, 250, 'both', '$[m/a]$', [-250,-100,-50,-25,-10,-5,-1,0,1,5,10,25,50,100,250], 'epsg:3031', 1000, symlog=True)
    elif i==9:
        # v_base
        plot_name = 'v_b.png'
        plot_data('refined_new2_vel_base.csv',plot_name,'Spectral_r', cmin=1, cmax=1000, cb_extend='max', cb_label='$[m/a]$', cb_ticks=[0.1,1,5,10,25,50,100,250,500,800], EPSG='epsg:3031', res=1000, log=True)
    elif i==10:
        # v_base/v_surf
        # TODO colorbar von Thomas, input data other format -> three columns
        plot_name = 'v_b_v_sim_ratio.png'
        plot_data('refined_new2_vel_surf_base_ratio.csv',plot_name,'winter', 0, 10, 'max', ' ', [0,1,2,5], 'epsg:3031', 100)
    elif i==11:
        # k_diff interp
        # TODO
        plot_name = 'k_diff_interp.png'
        plot_data('refined_new2_k.csv',plot_name,'Blues_r', 0, 150, 'max', '$[(s/m)^{(1/2)}]$', [0,50,100,150], 'epsg:3031', 100)
    elif i==12:
        # k_diff filter 
        # TODO
        plot_name = 'k_diff_filter.png'
        plot_data('refined_new2_k.csv',plot_name,'Blues_r', 0, 150, 'max', '$[(s/m)^{(1/2)}]$', [0,50,100,150], 'epsg:3031', 100)
    else:
        fname='other.png'
        print(i)


    # peak
    peak_df = gpd.read_file('peak.shp')
    peak_df.loc[:, 'x'] = peak_df.geometry.apply(lambda y: y.x)
    peak_df.loc[:, 'y'] = peak_df.geometry.apply(lambda y: y.y)

    x, y = sm.grid.transform(peak_df.x.values, peak_df.y.values, crs=peak_df.crs)

    ax.scatter(x, y, color='0.5', s=10, linewidths=1)

    # power
    power_df = gpd.read_file('power.shp')
    power_df.loc[:, 'x'] = power_df.geometry.apply(lambda y: y.x)
    power_df.loc[:, 'y'] = power_df.geometry.apply(lambda y: y.y)

    x, y = sm.grid.transform(power_df.x.values, power_df.y.values, crs=power_df.crs)

    ax.scatter(x, y, color='purple', s=10, linewidths=1)


    # lonlat grid
    sm.set_lonlat_contours(add_ytick_labels=True, xinterval=10, yinterval=2, linewidths=0.8,
                         linestyles='dashed', colors='0.2')
    sm.visualize(ax=ax, addcbar=False)

    #legend
    custom_lines = [Line2D([0], [0], color='k', lw=2),
                    Line2D([0], [0], color='k', lw=2,linestyle='dashed'),
                    Line2D([0], [0], color='w', lw=2),
		    Line2D([0], [0], color='w', linewidth=0,marker='o', markerfacecolor='0.5',lw=2,linestyle='',markeredgecolor='0.5'),
		    Line2D([0], [0], color='w', linewidth=0,marker='o', markerfacecolor='purple',lw=2, linestyle='',markeredgecolor='purple')
		    ]


    ax.legend(custom_lines, ['lakes', 'shear margin', 'sinks','peak','power'],framealpha=0.7,facecolor='0.9',loc=3)

    #plt.show()
    plt.savefig(plot_name, dpi=300, format='png')
