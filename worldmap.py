# -*- coding: utf-8 -*-
"""
Class to plot maps

Author: Milon Miah
"""

import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import numpy as np
import xarray as xr
from matplotlib.animation import FuncAnimation

class WorldMap(object):
    
    #Initializor
    #Define projection type, central longitude, 
    def __init__(self, proj_type = 'PC', 
                 central_longitude = 0, 
                 central_latitude = 0, 
                 num_fig = None, 
                 gridlines = True, 
                 real_color = True, 
                 zoom = None,
                 fig_size = None,
                 ax = None,
                 fig = None):
        
        #Save variables for internal use
        self.Proj_type = proj_type
        self.Central_longitude = central_longitude
        self.Central_latitude = central_latitude
        self.Num_fig = num_fig
        self.Gridlines = gridlines
        self.Real_color = real_color
        self.Zoom = zoom
        self.Fig_size = fig_size

        #Initialize projection type 
        proj_obj = [ccrs.PlateCarree(central_longitude = central_longitude),
                    ccrs.Robinson(central_longitude = central_longitude),
                    ccrs.Mercator(central_longitude = central_longitude),
                    ccrs.Orthographic(central_longitude = central_longitude, central_latitude = central_latitude),
                    ccrs.InterruptedGoodeHomolosine(central_longitude = central_longitude),
                    ccrs.AlbersEqualArea(central_longitude = central_longitude, central_latitude = central_latitude)]
        
        proj_types = ['PC','R','M','O','IGM', 'AEA']
        
        #Plot map according to projection
        if (proj_type in proj_types) == True:
            for i in range(len(proj_types)):
                if proj_type == proj_types[i]:
                    
                    if ax == None:
                        self.Fig = plt.figure(figsize = self.Fig_size)
                        self.Ax = self.Fig.add_subplot(1, 1, 1, projection = proj_obj[i])
                    else:
                        self.Ax = ax
                        
                    #Add optionally gridlines 
                    if gridlines == True:
                        self.Ax.gridlines(draw_labels = True)
                    
                    #Add optionally realistic colors  
                    if real_color == True:
                        self.Ax.stock_img()
                        
                    #Zoom in optionally
                    if zoom != None:
                        self.Ax.set_extent(zoom)
                    
                    #Plot coastlines
                    self.Ax.coastlines()
                    
                    #Label axes
                    self.Ax.set_xlabel('Longitude')
                    self.Ax.set_ylabel('Latitude')
                    
                    break 
                
        else:
            raise Exception('Invalid projection type!')
            
    def plot_colormesh(self, x, y, z, vmin = None, vmax = None):

       
        """
         Add a colormesh plot to WorldMap
         - Optionally algorithm with autoscale could be used (Takes a lot of computational capacities)
        
        """
        #
        # if (type(x) or type(y)) != np.ndarray:
        #     x = x.to_numpy()
        #     y = y.to_numpy()
        
        #Find minimum and maximum of x and y
        # min_x = x[0][0]
        # max_x = x[0][0]
        # min_y = y[0][0]
        # max_y = y[0][0]
        
        # for i in range(len(x)):
        #     for j in range(len(x[i])):
        #         if x[i][j] < min_x:
        #             min_x = x[i][j]
        #         if x[i][j] > max_x:
        #             max_x = x[i][j]
                    
        # for i in range(len(y)):
        #     for j in range(len(y[i])):
        #         if y[i][j] < min_y:
        #             min_y = y[i][j]
        #         if y[i][j] > max_y:
        #             max_y = y[i][j]
                            
        # self.Ax.set_extent([min_x, max_x, min_y, max_y])
        
        plot2D = self.Ax.pcolormesh(x, y, z, transform = ccrs.PlateCarree(), vmin = vmin, vmax = vmax)
        self.Ax.coastlines()
        
        return plot2D
        
    def plot_quiver(self, x, y, u, v, skip = None):
        
        """
        Add quiver plot to WorldMap

        """
        if skip != None:
            to_skip = (slice(None, None, skip), slice(None, None, skip))
            x = x[to_skip]
            y = y[to_skip]
            u = u[to_skip]
            v = v[to_skip]
        
        plot2D = self.Ax.quiver(x.to_numpy(), y.to_numpy(), u.to_numpy(), v.to_numpy(), transform = ccrs.PlateCarree())
        self.Ax.coastlines()

        return plot2D

    def plotByGroup_T(self, dataset, group, choice, subplot = True, vmin = None, vmax = None):
        
        """
        Function to plot grouped data (e.g. for seasonal cycle, here group = 'month')
        - Choice: Chosen years
        - For temperature plot 
        - Not very general function
        """

        dataset_grouped = list(dataset.groupby('{}'.format(group)))
        
        #Convert year into index
        init_var = dataset_grouped[0][0]
        dataset_choice = dataset_grouped[choice-init_var][1] 
        
        if subplot == True:
            
            fig = plt.figure(figsize = (30,20), constrained_layout=True)
            axs = []
            for i in range(1, 13):
                ax = fig.add_subplot(3, 4, i, projection = ccrs.PlateCarree(self.Central_longitude))
                ax.set_title('{}/'.format(i) + '{}'.format(choice))
                atl_group = WorldMap(proj_type = self.Proj_type,
                                     central_longitude = self.Central_longitude, 
                                     central_latitude = self.Central_latitude, 
                                     num_fig = self.Num_fig,
                                     real_color = self.Real_color,
                                     gridlines = self.Gridlines, 
                                     fig_size = self.Fig_size,
                                     zoom = self.Zoom,
                                     ax = ax)
                
                color_plot = atl_group.plot_colormesh(dataset_choice.isel(time_counter = i-1).nav_lon,
                                                  dataset_choice.isel(time_counter = i-1).nav_lat,
                                                  dataset_choice.isel(time_counter = i-1).votemper,
                                                  vmin = vmin, vmax = vmax)
                if i % 4 == 0:
                    axs.append(ax)
            plt.colorbar(color_plot, ax = axs, label = 'Sea Surface Temperature [°C]', shrink = 0.8, location='right')
            plt.subplots_adjust(wspace = 1, hspace = 1)                             
            
        else:
            for i in range(len(dataset_choice.time_counter)):  
                
                atl_group = WorldMap(proj_type = self.Proj_type,
                                     central_longitude = self.Central_longitude, 
                                     central_latitude = self.Central_latitude, 
                                     num_fig = self.Num_fig,
                                     real_color = self.Real_color,
                                     gridlines = self.Gridlines, 
                                     fig_size = self.Fig_size,
                                     zoom = self.Zoom)
                
        
                color_plot = atl_group.plot_colormesh(dataset_choice.isel(time_counter = i).nav_lon,
                                                 dataset_choice.isel(time_counter = i).nav_lat,
                                                 dataset_choice.isel(time_counter = i).votemper,
                                                 vmin = vmin, vmax = vmax)
                            
                
                plt.colorbar(color_plot, label = 'Sea Surface Temperature [°C]')
                plt.title('{}-{}'.format(i+1, choice))
                
    def plotByGroup_S(self, dataset, group, choice, vmin = None, vmax = None, subplot = True):
        
        """
        Function to plot grouped data (e.g. for seasonal cycle, here group = 'month')
        - Choice: Chosen years
        - Salinity plots
        - Not very general function

        """
        dataset_grouped = list(dataset.groupby('{}'.format(group)))
        
        #Convert year into index
        init_var = dataset_grouped[0][0]
        dataset_choice = dataset_grouped[choice-init_var][1] 
        
        if subplot == True:
            fig = plt.figure(figsize = (30,20), constrained_layout = True)
            axs = []
            for i in range(1, 13):
                ax = fig.add_subplot(3, 4, i, projection = ccrs.PlateCarree(self.Central_longitude))
                ax.set_title('{}/'.format(i) + '{}'.format(choice))
                atl_group = WorldMap(proj_type = self.Proj_type,
                                     central_longitude = self.Central_longitude, 
                                     central_latitude = self.Central_latitude, 
                                     num_fig = self.Num_fig,
                                     real_color = self.Real_color,
                                     gridlines = self.Gridlines, 
                                     fig_size = self.Fig_size,
                                     zoom = self.Zoom,
                                     ax = ax)
                
                color_plot = atl_group.plot_colormesh(dataset_choice.isel(time_counter = i-1).nav_lon,
                                                  dataset_choice.isel(time_counter = i-1).nav_lat,
                                                  dataset_choice.isel(time_counter = i-1).vosaline,
                                                  vmin = vmin, vmax = vmax)
                
                if i % 4 == 0:
                    axs.append(ax)
                
            plt.colorbar(color_plot, ax = axs, label = 'Sea Surface Salinity [g/kg]', shrink = 0.8, location='right')
            plt.subplots_adjust(wspace = 1, hspace = 1)   
        else:
        
            for i in range(len(dataset_choice.time_counter)):
                atl_group = WorldMap(proj_type = self.Proj_type,
                                     central_longitude = self.Central_longitude, 
                                     central_latitude = self.Central_latitude, 
                                     num_fig = self.Num_fig,
                                     real_color = self.Real_color,
                                     gridlines = self.Gridlines, 
                                     fig_size = self.Fig_size,
                                     zoom = self.Zoom)
                
                color_plot = atl_group.plot_colormesh(dataset_choice.isel(time_counter = i).nav_lon,
                                                 dataset_choice.isel(time_counter = i).nav_lat,
                                                 dataset_choice.isel(time_counter = i).vosaline,
                                                 vmin = vmin, vmax = vmax)
                
                plt.colorbar(color_plot, label = 'Sea Surface Salinity [g/kg]')
                plt.title('{}-{}'.format(i+1, choice))
            
    def time_mean_field(self, dataset, z, vmin = None, vmax = None, subplot = True):
        
        """
        Function to plot time-mean field
        - z: Property to plot

        """

        dataset_month = list(dataset.groupby('time_counter.month'))
        
        if z == 'votemper':
            label = 'Temperature [°C]'
        else:
            label = 'Salinity [g/kg]'

        if subplot == True:
            fig = plt.figure(figsize = (30,20), constrained_layout = True)
            axs = []
            
            for i in range(1, 13):
                
                data_month_ave = dataset_month[i-1][1].mean(dim = 'time_counter', skipna = True)
                ax = fig.add_subplot(3, 4, i, projection = ccrs.PlateCarree(self.Central_longitude))
                ax.set_title('Month {}'.format(i))
                atl_month = WorldMap(proj_type = self.Proj_type,
                                     central_longitude = self.Central_longitude, 
                                     central_latitude = self.Central_latitude, 
                                     num_fig = self.Num_fig,
                                     real_color = self.Real_color,
                                     gridlines = self.Gridlines, 
                                     fig_size = self.Fig_size,
                                     zoom = self.Zoom,
                                     ax = ax)
                
                color_plot = atl_month.plot_colormesh(data_month_ave.nav_lon,
                                                 data_month_ave.nav_lat,
                                                 data_month_ave['{}'.format(z)],
                                                 vmin = vmin, vmax = vmax)
                
                if i % 4 == 0:
                    axs.append(ax)
                
            plt.colorbar(color_plot, ax = axs, label = 'Sea Surface {}'.format(label), shrink = 0.8, location='right')
            plt.subplots_adjust(wspace = 1, hspace = 1)   
                
        else:
            for i in range(len(dataset_month)):
                data_month_ave = dataset_month[i][1].mean(dim = 'time_counter', skipna = True)
                
                atl_month = WorldMap(proj_type = self.Proj_type,
                                     central_longitude = self.Central_longitude, 
                                     central_latitude = self.Central_latitude, 
                                     num_fig = self.Num_fig,
                                     real_color = self.Real_color,
                                     gridlines = self.Gridlines, 
                                     fig_size = self.Fig_size,
                                     zoom = self.Zoom)
                
                color_plot = atl_month.plot_colormesh(data_month_ave.nav_lon,
                                                 data_month_ave.nav_lat,
                                                 data_month_ave['{}'.format(z)],
                                                 vmin = vmin, vmax = vmax)
                
                plt.colorbar(color_plot, label = 'Sea Surface {}'.format(label))
                plt.title('{}'.format(i+1)) 
                
    def add_bathymetry(self, z):
        
        z = float(z)
        bmty = xr.open_dataset('DV_NAshelf_NEMO//Data//bathymetry_NWA.nc')
        
        bthy = self.Ax.contour(bmty.x, bmty.y, bmty.z, levels = [z], colors = 'black', label = '{} m'.format(-1 * z))
        self.Ax.clabel(bthy, inline = 2, fontsize = 15)

        return bthy
    
