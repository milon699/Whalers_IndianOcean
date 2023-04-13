#!/usr/bin/env python
# coding: utf-8

#%%

#Import packages

import pickle
import numpy as np
from worldmap import WorldMap
import matplotlib.pyplot as plt
import cmocean as cmo
import cartopy.crs as ccrs
import matplotlib.patches as mpatches
from windrose import WindroseAxes
from matplotlib.animation import FuncAnimation
import cartopy.feature as cfeature
from matplotlib import colors

#Ignore warnings for now

import warnings
warnings.filterwarnings('ignore')

#Set Plot parameters
plt.close('all')
plt.rcParams['font.size'] = 20.0
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rcParams['font.weight'] = 'bold'
plt.rcParams['axes.labelsize'] = 'medium'
plt.rcParams['axes.labelweight'] = 'bold'
plt.rcParams['axes.linewidth'] = 1
plt.rcParams['lines.linewidth'] = 1

plt.close('all')

#%%
#Load data - Data Update 13-01-2023
with open('Shiplog Data/logentries-export-2023-01-13-clean_tier1.pkl', 'rb') as file:
    data = pickle.load(file)

#%%
#Drop not usable latitudes and longitudes, wind speed and wind directions
data_ubl = data.drop(data[data.usable_WD == False].index)

#%%
#To begin with plot world map with all the datapoints
worldmap = WorldMap(fig_size = (30,15))
worldmap.Ax.scatter(data_ubl.Longitude, data_ubl.Latitude, color = 'r', marker = 'x', s = 1)
worldmap.Fig.savefig('Plots/all_data_points.png')

#%%
#Plot and safe global wind speeds on new world map
worldmap2 = WorldMap(fig_size = (30,15), real_color = False)
sctt = worldmap2.Ax.scatter(data_ubl.Longitude, data_ubl.Latitude, s = 1, c = data_ubl['Wind Speed/Force'], cmap = 'cmo.speed')
worldmap2.Fig.colorbar(sctt, ax = worldmap2.Ax, label = 'Wind Speed (Beaufort)', location = 'right', fraction = 0.02, pad = 0.1)
worldmap2.Fig.savefig('Plots/windspeed_global.png')

#%% 
#Account for the number of data points within the specified region
data_ubl['usable_WD'] = data_ubl["Wind Direction"].notna()
print('There are {} usable data points in the Indian Ocean\n'.format(len(data_ubl.query('-60 <= Latitude <= 20 & 20 <= Longitude <= 130').loc[(data_ubl['usable_WD']==True), "Entry Date Time"])))
print('This are {:.2f} % of all data points'.format(len(data_ubl.query('-60 <= Latitude <= 20 & 20 <= Longitude <= 130').loc[(data_ubl['usable_WD']==True), "Entry Date Time"])/len(data_ubl.loc[(data_ubl['usable_WD']==True), "Entry Date Time"])*100))

#%% 
#Plot wind vectors in Indian Ocean in specific time periods
plot_yr = True
if plot_yr == True:
    
    #Set True, if plots should be created seperately for each year or False, if 3x4 subplots are desired
    onebyone = False
    if onebyone == True:
        
        #Define observation period
        time_begin = 1835
        time_end = 1913
        
        #Cut out data of Indian Ocean
        data_IO = data_ubl.query('-60 <= Latitude <= 20 & 20 <= Longitude <= 130').groupby(data_ubl['Entry Date Time'].dt.year)
        
        for i in range(time_begin, time_end):
            
            #Create figure
            fig = plt.figure(figsize = (30,15))
            ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
            
            plt.title(i)
            ax.add_feature(cfeature.COASTLINE)
            ax.set_extent([20, 130, -60, 20])
            ax.gridlines(draw_labels = True)
               
            #If there is data, carry out the quiver plot, otherwise just continue with loop
            try:
                quiv = ax.quiver(data_IO.get_group(i).Longitude, data_IO.get_group(i).Latitude, 
                                 -data_IO.get_group(i)['Wind Speed/Force']*np.sin(data_IO.get_group(i)['Wind Direction'].to_numpy().astype(float)*2*np.pi/360), 
                                 -data_IO.get_group(i)['Wind Speed/Force']*np.cos(data_IO.get_group(i)['Wind Direction'].to_numpy().astype(float)*2*np.pi/360))
            except KeyError:
                continue
        
    else:  
        
        #Set starting year
        time_begin = 1840
        time_end = time_begin + 12
        
        #Cut out data of Indian Ocean
        data_IO = data_ubl.query('-60 <= Latitude <= 20 & 20 <= Longitude <= 130').groupby(data_ubl['Entry Date Time'].dt.year)
        
        fig = plt.figure(figsize = (50,25))
        
        for i in range(time_begin, time_end):
           
            #Add subplot to figure
            ax = fig.add_subplot(3, 4, i-time_begin+1, projection=ccrs.PlateCarree())
            ax.add_feature(cfeature.COASTLINE)
            ax.set_extent([20, 130, -60, 20])
            ax.gridlines(draw_labels = True)
            ax.title.set_text(i)
            
            #Draw boxes for east and west of Indian Ocean
            rect_left=mpatches.Rectangle((35,-25), 25, 45, fill = False, color = "red", linewidth = 3)
            rect_right=mpatches.Rectangle((100,-35), 20, 35, fill = False, color = "red", linewidth = 3)
            ax.add_patch(rect_left)
            ax.add_patch(rect_right)
            
            #If there is data, carry out the quiver plot, otherwise just continue with loop
            try:
                quiv = ax.quiver(data_IO.get_group(i).Longitude, data_IO.get_group(i).Latitude, 
                                  -data_IO.get_group(i)['Wind Speed/Force']*np.sin(data_IO.get_group(i)['Wind Direction'].to_numpy().astype(float)*2*np.pi/360), 
                                  -data_IO.get_group(i)['Wind Speed/Force']*np.cos(data_IO.get_group(i)['Wind Direction'].to_numpy().astype(float)*2*np.pi/360))
            except KeyError:
                continue
        
        #Save plot with 12 panels
        fig.savefig('Plots/Quiver_Plot_{}_{}.png'.format(time_begin, time_end))
                
#%%
#Create quiver plot of wind vectors in entire Indian Ocean
quiver_plot = True
if quiver_plot == True:

    #Create new worldmap
    worldmap3 = WorldMap(fig_size = (30, 15), real_color = False, zoom = [20, 130, -60, 20])
    
    #Plot boxes for east and west Indian Ocean
    #Specify the location with (left, bottom), width, height
    rect_left=mpatches.Rectangle((35,-25), 25, 45, fill = False, color = "red", linewidth = 3)
    rect_right=mpatches.Rectangle((100,-35), 20, 35, fill = False, color = "red", linewidth = 3)
    
    worldmap3.Ax.add_patch(rect_left)
    worldmap3.Ax.add_patch(rect_right)
    
    #Create quiver plot for all data points in Indian ocean
    quiv = worldmap3.Ax.quiver(data_ubl.Longitude, data_ubl.Latitude, -data_ubl['Wind Speed/Force']*np.sin(data_ubl['Wind Direction'].to_numpy().astype(float)*2*np.pi/360),
                               -data_ubl['Wind Speed/Force']*np.cos(data_ubl['Wind Direction'].to_numpy().astype(float)*2*np.pi/360), 
                               scale = 300)
    
    #Save figure
    worldmap3.Fig.savefig('Plots/IO_windvectors.png')

#%%
#Create quiver plot for all data points in a certain year
year_spe = True
if year_spe == True:
    year_input = int(input('Which year would you like to plot?\n'))
    valid_year = False
    while valid_year == False:
        try:
            data_yr = data_ubl.query('-60 <= Latitude <= 20 & 20 <= Longitude <= 130').groupby(data_ubl['Entry Date Time'].dt.year).get_group(year_input)
            valid_year = True
            print('There are {} data points in {} in the Indian Ocean'.format(len(data_yr), year_input))
        except KeyError:
            year_input= int(input('No data points in Indian Ocean in this year, try another one:\n'))
    
    
    worldmap_yr = WorldMap(fig_size = (30, 15), real_color = False, zoom = [20, 130, -60, 20])
    plt.title('Data Points Indian Ocean {}'.format(year_input))
    rect_left=mpatches.Rectangle((35,-25), 25, 45, fill = False, color = "red", linewidth = 3)
    rect_right=mpatches.Rectangle((100,-35), 20, 35, fill = False, color = "red", linewidth = 3)

    worldmap_yr.Ax.add_patch(rect_left)
    worldmap_yr.Ax.add_patch(rect_right)
    quiv_yr = worldmap_yr.Ax.quiver(data_yr.Longitude, data_yr.Latitude, -data_yr['Wind Speed/Force']*np.sin(data_yr['Wind Direction'].to_numpy().astype(float)*2*np.pi/360),
                               -data_yr['Wind Speed/Force']*np.cos(data_yr['Wind Direction'].to_numpy().astype(float)*2*np.pi/360), 
                               scale = 300)

    worldmap_yr.Fig.savefig('Plots/IO_windvectors_{}.png'.format(year_input))
    
#%%
#Create grid over Indian Ocean to assess data points per grid cell
assess_data = True
if assess_data == True:
    
    #Define starting year and do the plot later on for a certain decade
    start_year = 1880
    end_year = start_year + 10
    
    #Cut out data in time column
    data_dp = data_ubl[(data_ubl['Entry Date Time'].dt.year >= start_year) & (data_ubl['Entry Date Time'].dt.year <= end_year)]
    
    #Set boundary parameters and stepsize for the grid
    lat_min = -60
    lat_max = 20
    lon_min = 20
    lon_max = 130
    stepsize = 1
    lon = np.arange(lon_min, lon_max, stepsize)
    lat = np.arange(lat_min, lat_max, stepsize)        
    
    #Create corresponding grid
    grid = np.zeros((lat_max-lat_min, lon_max-lon_min))
    for i in range(lon_min, lon_max):
        for j in range(lat_min, lat_max):
            grid[j-lat_min][i-lon_min] = len(data_dp[(data_dp["Latitude"] >= j) & (data_dp["Latitude"] < j+1) 
                                      & (data_dp['Longitude'] >= i) & (data_dp['Longitude'] < i+1)])
     
    #Create colormesh with amount of data per grid cell
    fig = plt.figure(figsize = (30, 15))
    ax = fig.add_subplot(1,1,1, projection = ccrs.PlateCarree())
    ax.set_title('{}-{}'.format(start_year, end_year))
    ax.set_extent([20, 130, -60, 20])
    ax.gridlines(draw_labels = True)
    cmap = plt.get_cmap('binary')
    bounds = np.arange(0, 11, 1)
    norm = colors.BoundaryNorm(bounds, cmap.N)
    cmesh = ax.pcolormesh(lon, lat, grid, transform = ccrs.PlateCarree(), norm = norm, cmap = cmap)
    plt.colorbar(cmesh, ax = ax, location = 'right', label = 'Data Points', pad = 0.05)
    ax.coastlines()
    
    fig.savefig('Plots/data_points_{}_{}.png'.format(start_year, end_year))

#Optionally an animation can be created as well
create_animation = True
if create_animation == True:

    #Define observation period
    start_year = 1830
    end_year = 1920
    
    #Cut out data on time column
    data_dp = data_ubl[(data_ubl['Entry Date Time'].dt.year >= start_year) & (data_ubl['Entry Date Time'].dt.year <= (start_year + 1))]
    
    #Set boundary parameters and stepsize of grid
    lat_min = -60
    lat_max = 20
    lon_min = 20
    lon_max = 130
    stepsize = 1
    lon = np.arange(lon_min, lon_max, stepsize)
    lat = np.arange(lat_min, lat_max, stepsize)
    
    #Create actual grid
    grid = np.zeros((lat_max-lat_min, lon_max-lon_min))
    for i in range(lon_min, lon_max):
        for j in range(lat_min, lat_max):
            grid[j-lat_min][i-lon_min] = len(data_dp[(data_dp["Latitude"] >= j) & (data_dp["Latitude"] < j+1) 
                                      & (data_dp['Longitude'] >= i) & (data_dp['Longitude'] < i+1)])
    
    #Create initial colormesh plot
    fig = plt.figure(figsize = (20,20))
    ax = fig.add_subplot(1,1,1, projection = ccrs.PlateCarree())
    ax.set_extent([20, 130, -60, 20]) 
    ax.gridlines(draw_labels = True)
    
    cmap = plt.get_cmap('binary')
    bounds = np.arange(0, 11, 1)
    norm = colors.BoundaryNorm(bounds, cmap.N)
    cmesh = ax.pcolormesh(lon, lat, grid, transform = ccrs.PlateCarree(), norm = norm, cmap = cmap)
    plt.colorbar(cmesh, ax = ax, location = 'right', label = 'Data Points', pad = 0.1)
    ax.coastlines()
    
    #Function to update the colormesh
    def animate(i):
        
        plt.title(start_year + i)
        data_dp = data_ubl[(data_ubl['Entry Date Time'].dt.year >= (start_year + i)) & (data_ubl['Entry Date Time'].dt.year <= (start_year + i + 1))]
        for i in range(lon_min, lon_max):
            for j in range(lat_min, lat_max):
                grid[j-lat_min][i-lon_min] = len(data_dp[(data_dp["Latitude"] >= j) & (data_dp["Latitude"] < j+1) 
                                          & (data_dp['Longitude'] >= i) & (data_dp['Longitude'] < i+1)])
        cmesh.set_array(grid)
       
    #Create and save animation with corresponding fps
    anim = FuncAnimation(fig, animate, frames = end_year-start_year, interval = 5, blit = False)
    anim.save('Plots/data_points_anim.gif', writer = 'imagemagick', fps = 1)

#%%
#Create datasets for both boxes, west and east of Indian Ocean 
data_ubl_west = data_ubl.query('-25 <= Latitude <= 20 & 35 <= Longitude <= 60')  
data_ubl_east = data_ubl.query('-35 <= Latitude <= -5 & 100 <= Longitude <= 120')  
 
#Try to group data in climatological cycle in both boxes  
group_data = True
if group_data == True:
    
    #Group data by months
    gr_lt = data_ubl_west.groupby(data_ubl_west['Entry Date Time'].dt.month)
    gr_rt = data_ubl_east.groupby(data_ubl_east['Entry Date Time'].dt.month)
    gr = data_ubl.groupby(data_ubl['Entry Date Time'].dt.month)
    
    #Create list of windroses axes and maps
    ax_lt = []
    ax_rt = []
    worldmaps = []
    months = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Okt', 'Nov', 'Dec'] 
    
    for i in range(12):
        
        #Create windrose west
        ax_lt.append(WindroseAxes.from_ax())
        plt.title('West IO: {}'.format(months[i]))
        ax_lt[i].bar(gr_lt['Wind Direction'].get_group(i+1), gr_lt['Wind Speed/Force'].get_group(i+1), 
                      normed=True, opening=0.8, edgecolor='white')
        ax_lt[i].set_legend()
        fig_lt = plt.gcf()
        
        #Create windrose east
        ax_rt.append(WindroseAxes.from_ax())
        plt.title('East IO: {}'.format(months[i]))
        ax_rt[i].bar(gr_rt['Wind Direction'].get_group(i+1), gr_rt['Wind Speed/Force'].get_group(i+1), 
                      normed=True, opening=0.8, edgecolor='white')
        ax_rt[i].set_legend()
        fig_rt = plt.gcf()
        
        #Create map of Indian Ocean with boxes and add data for corresponding time period
        worldmaps.append(WorldMap(fig_size = (30, 15), real_color = False, zoom = [20, 130, -60, 20]))
        rect_left=mpatches.Rectangle((35,-25), 25, 45, fill = False, color = "red", linewidth = 3)
        rect_right=mpatches.Rectangle((100,-35), 20, 35, fill = False, color = "red", linewidth = 3)
        
        quiv = worldmaps[i].Ax.quiver(gr['Longitude'].get_group(i+1), gr['Latitude'].get_group(i+1), -gr['Wind Speed/Force'].get_group(i+1)*np.sin(gr['Wind Direction'].get_group(i+1).to_numpy().astype(float)*2*np.pi/360),
                                    -gr['Wind Speed/Force'].get_group(i+1)*np.cos(gr['Wind Direction'].get_group(i+1).to_numpy().astype(float)*2*np.pi/360), 
                                    scale = 300)
        
        plt.title('{}'.format(months[i]))
        worldmaps[i].Ax.add_patch(rect_left)
        worldmaps[i].Ax.add_patch(rect_right)

        #Save figures
        fig_lt.savefig('Plots/WestIO_windrose_{}.png'.format(months[i]))
        fig_rt.savefig('Plots/EastIO_windrose_{}.png'.format(months[i]))
        worldmaps[i].Fig.savefig('Plots/Windvectors_{}.png'.format(months[i]))
        
#%%
#Try to map Mascarene High seasonally and for certain decades with data
worldmap = WorldMap(fig_size = (30, 15), real_color = False, zoom = [20, 130, -55, -0])

#Plot boxes for Mascarene High (south, north, east, west)
# specify the location with (left, bottom), width, height
w = 45
h = 15
s = 58
rect_south=mpatches.Rectangle((s,-50), w, h, fill = False, color = "red", linewidth = 3)
rect_north=mpatches.Rectangle((s,-35), w, h, fill = False, color = "red", linewidth = 3)
rect_east=mpatches.Rectangle((s+w,-50), 8, 2*h, fill = False, color = "red", linewidth = 3)
rect_west=mpatches.Rectangle((s-8,-50), 8, 2*h, fill = False, color = "red", linewidth = 3)

#Add boxes to plot
worldmap.Ax.add_patch(rect_south)
worldmap.Ax.add_patch(rect_north)
worldmap.Ax.add_patch(rect_east)
worldmap.Ax.add_patch(rect_west)

#Add all data points in area
quiv = worldmap.Ax.quiver(data_ubl.Longitude, data_ubl.Latitude, -data_ubl['Wind Speed/Force']*np.sin(data_ubl['Wind Direction'].to_numpy().astype(float)*2*np.pi/360),
                           -data_ubl['Wind Speed/Force']*np.cos(data_ubl['Wind Direction'].to_numpy().astype(float)*2*np.pi/360), 
                           scale = 300)

#Save figure
worldmap.Fig.savefig('Plots/MH_all_data_pount.png')

#%%
#Split data for the different boxes
data_MH_south = data_ubl.query('-50 <= Latitude <= -35 & 60 <= Longitude <= 105')
data_MH_north = data_ubl.query('-35 <= Latitude <= -20 & 60 <= Longitude <= 105')
data_MH_west = data_ubl.query('-50 <= Latitude <= -20 & 50 <= Longitude <= 60')
data_MH_east = data_ubl.query('-50 <= Latitude <= -20 & 105 <= Longitude <= 115')

#Create windroses of data for all years
ax_south = WindroseAxes.from_ax()
plt.title('Mascarene High Southern Box')
ax_south.bar(data_MH_south['Wind Direction'], data_MH_south['Wind Speed/Force'], normed=True, opening=0.8, edgecolor='white')
ax_south.set_legend()
ax_south.set_yticks(np.arange(4, 24, step=4))
ax_south.set_yticklabels(np.arange(4, 24, step=4))
fig_S = plt.gcf()

ax_north = WindroseAxes.from_ax()
plt.title('Mascarene High Northern Box')
ax_north.bar(data_MH_north['Wind Direction'], data_MH_north['Wind Speed/Force'], normed=True, opening=0.8, edgecolor='white')
ax_north.set_legend()
ax_north.set_yticks(np.arange(4, 24, step=4))
ax_north.set_yticklabels(np.arange(4, 24, step=4))
fig_N = plt.gcf()

ax_east = WindroseAxes.from_ax()
plt.title('Mascarene High Eastern Box')
ax_east.bar(data_MH_east['Wind Direction'], data_MH_east['Wind Speed/Force'], normed=True, opening=0.8, edgecolor='white')
ax_east.set_legend()
ax_east.set_yticks(np.arange(4, 24, step=4))
ax_east.set_yticklabels(np.arange(4, 24, step=4))
fig_E = plt.gcf()

ax_west = WindroseAxes.from_ax()
plt.title('Mascarene High Western Box')
ax_west.bar(data_MH_west['Wind Direction'], data_MH_west['Wind Speed/Force'], normed=True, opening=0.8, edgecolor='white')
ax_west.set_legend()
ax_west.set_yticks(np.arange(4, 24, step=4))
ax_west.set_yticklabels(np.arange(4, 24, step=4))
fig_W = plt.gcf()

#Save figures
fig_S.savefig('Plots/windrose_MH_south.png')
fig_N.savefig('Plots/windrose_MH_north.png')
fig_E.savefig('Plots/windrose_MH_east.png')
fig_W.savefig('Plots/windrose_MH_west.png')

#%%
#Try to also do a climatological analysis on data, group data by either months or quarters
gr_south = data_MH_south.groupby(data_MH_south['Entry Date Time'].dt.month)
gr_north = data_MH_north.groupby(data_MH_north['Entry Date Time'].dt.month)
gr_west = data_MH_west.groupby(data_MH_west['Entry Date Time'].dt.month)
gr_east = data_MH_east.groupby(data_MH_east['Entry Date Time'].dt.month)
gr = data_ubl.groupby(data_ubl['Entry Date Time'].dt.month)

#Either analyze on monthly or quarterly base
three_month_interval = True
if three_month_interval == False:
    months = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Okt', 'Nov', 'Dec'] 
    N = 12
    for i in range(N):
        
        #Create new worldmap and add Mascarene High boxes
        worldmap_gr = WorldMap(fig_size = (30, 15), real_color = False, zoom = [20, 130, -55, -0])
        plt.title('Month ' + months[i], fontsize = 20)
        w = 45
        h = 15
        s = 58
        rect_south=mpatches.Rectangle((s,-50), w, h, fill = False, color = "red", linewidth = 3)
        rect_north=mpatches.Rectangle((s,-35), w, h, fill = False, color = "red", linewidth = 3)
        rect_east=mpatches.Rectangle((s+w,-50), 8, 2*h, fill = False, color = "red", linewidth = 3)
        rect_west=mpatches.Rectangle((s-8,-50), 8, 2*h, fill = False, color = "red", linewidth = 3)
        worldmap_gr.Ax.add_patch(rect_south)
        worldmap_gr.Ax.add_patch(rect_north)
        worldmap_gr.Ax.add_patch(rect_east)
        worldmap_gr.Ax.add_patch(rect_west)
        
        #Create quiver plot for data points in corresponding month
        quiv = worldmap_gr.Ax.quiver(gr.Longitude.get_group(i+1), gr.Latitude.get_group(i+1), -gr['Wind Speed/Force'].get_group(i+1)*np.sin(gr['Wind Direction'].get_group(i+1).to_numpy().astype(float)*2*np.pi/360),
                                     -gr['Wind Speed/Force'].get_group(i+1)*np.cos(gr['Wind Direction'].get_group(i+1).to_numpy().astype(float)*2*np.pi/360),
                                     scale = 300)
            
        #Check if there are actually data points for each month (otherwise quiver plot will lead to error)
        if i+1 in gr_south.groups:
            ax_south = WindroseAxes.from_ax()
            plt.title('Mascarene High Southern Box ' + months[i])
            ax_south.bar(gr_south['Wind Direction'].get_group(i+1), gr_south['Wind Speed/Force'].get_group(i+1), normed=True, opening=0.8, edgecolor='white')
            ax_south.set_legend()
            fig_S = plt.gcf()
            fig_S.savefig('Plots/Windrose_MH_south_{}.png'.format(months[i]))
            
        if i+1 in gr_north.groups:
            ax_north = WindroseAxes.from_ax()
            plt.title('Mascarene High Northern Box ' + months[i])
            ax_north.bar(gr_north['Wind Direction'].get_group(i+1), gr_north['Wind Speed/Force'].get_group(i+1), normed=True, opening=0.8, edgecolor='white')
            ax_north.set_legend()
            fig_N = plt.gcf()
            fig_N.savefig('Plots/Windrose_MH_north_{}.png'.format(months[i]))
        
        if i+1 in gr_east.groups:
            ax_east = WindroseAxes.from_ax()
            plt.title('Mascarene High Eastern Box ' + months[i])
            ax_east.bar(gr_east['Wind Direction'].get_group(i+1), gr_east['Wind Speed/Force'].get_group(i+1), normed=True, opening=0.8, edgecolor='white')
            ax_east.set_legend()
            fig_E = plt.gcf()
            fig_E.savefig('Plots/Windrose_MH_east_{}.png'.format(months[i]))
            
        if i+1 in gr_west.groups:
            ax_west = WindroseAxes.from_ax()
            plt.title('Mascarene High Western Box ' + months[i])
            ax_west.bar(gr_west['Wind Direction'].get_group(i+1), gr_west['Wind Speed/Force'].get_group(i+1), normed=True, opening=0.8, edgecolor='white')
            ax_west.set_legend()
            fig_W = plt.gcf()            
            fig_W.savefig('Plots/Windrose_MH_west_{}.png'.format(months[i]))
        
        #Save Figures
        worldmap_gr.Fig.savefig('Plots/Windvectors_MH_{}.png'.format(months[i]))
        
else:
    
    #Set order to DJF, MAM, JJA, SON
    order = [12,1,2,3,4,5,6,7,8,9,10,11]
    months = ['DJF', 'MAM', 'JJA', 'SON'] 
    N = 4
    for i in range(N):
        
        #Create map with Mascarene High boxes
        worldmap_gr = WorldMap(fig_size = (30, 15), real_color = False, zoom = [20, 130, -55, -0])       
        plt.title('Months ' + months[i], fontsize = 20)
        w = 45
        h = 15
        s = 58
        rect_south=mpatches.Rectangle((s,-50), w, h, fill = False, color = "red", linewidth = 3)
        rect_north=mpatches.Rectangle((s,-35), w, h, fill = False, color = "red", linewidth = 3)
        rect_east=mpatches.Rectangle((s+w,-50), 8, 2*h, fill = False, color = "red", linewidth = 3)
        rect_west=mpatches.Rectangle((s-8,-50), 8, 2*h, fill = False, color = "red", linewidth = 3)
        worldmap_gr.Ax.add_patch(rect_south)
        worldmap_gr.Ax.add_patch(rect_north)
        worldmap_gr.Ax.add_patch(rect_east)
        worldmap_gr.Ax.add_patch(rect_west)
        
        ax_south = WindroseAxes.from_ax()
        ax_north = WindroseAxes.from_ax()
        ax_east = WindroseAxes.from_ax()
        ax_west = WindroseAxes.from_ax()
        
        #Create quiver plot for months in interval (loop over three months)
        for j in range(3):
            
            quiv = worldmap_gr.Ax.quiver(gr.Longitude.get_group(order[3*i+j]), gr.Latitude.get_group(order[3*i+j]), -gr['Wind Speed/Force'].get_group(order[3*i+j])*np.sin(gr['Wind Direction'].get_group(order[3*i+j]).to_numpy().astype(float)*2*np.pi/360),
                                         -gr['Wind Speed/Force'].get_group(order[3*i+j])*np.cos(gr['Wind Direction'].get_group(order[3*i+j]).to_numpy().astype(float)*2*np.pi/360),
                                         scale = 300)
            
            #Create corrrsponding wind roses
            if order[3*i+j] in gr_south.groups:
                
                ax_south.bar(gr_south['Wind Direction'].get_group(order[3*i+j]), gr_south['Wind Speed/Force'].get_group(order[3*i+j]), normed=True, opening=0.8, edgecolor='white')
                ax_south.set_legend()
                ax_south.set_title('Mascarene High Southern Box ' + months[i])
                fig_S = plt.gcf()
                
            if order[3*i+j] in gr_north.groups:
                
                ax_north.bar(gr_north['Wind Direction'].get_group(order[3*i+j]), gr_north['Wind Speed/Force'].get_group(order[3*i+j]), normed=True, opening=0.8, edgecolor='white')
                ax_north.set_legend()
                ax_north.set_title('Mascarene High Northern Box ' + months[i])
                fig_N = plt.gcf()
                
            if order[3*i+j] in gr_east.groups:
                
                ax_east.bar(gr_east['Wind Direction'].get_group(order[3*i+j]), gr_east['Wind Speed/Force'].get_group(order[3*i+j]), normed=True, opening=0.8, edgecolor='white')
                ax_east.set_legend()
                ax_east.set_title('Mascarene High Eastern Box ' + months[i])
                fig_E = plt.gcf()
                
            if order[3*i+j] in gr_west.groups:
                
                ax_west.bar(gr_west['Wind Direction'].get_group(order[3*i+j]), gr_west['Wind Speed/Force'].get_group(order[3*i+j]), normed=True, opening=0.8, edgecolor='white')
                ax_west.set_legend()
                ax_west.set_title('Mascarene High Western Box ' + months[i])
                fig_W = plt.gcf() 
                
        #Save Figures
        worldmap_gr.Fig.savefig('Plots/Windvectors_MH_{}.png'.format(months[i]))
        fig_S.savefig('Plots/Windrose_MH_south_{}.png'.format(months[i]))
        fig_N.savefig('Plots/Windrose_MH_north_{}.png'.format(months[i]))
        fig_E.savefig('Plots/Windrose_MH_east_{}.png'.format(months[i]))
        fig_W.savefig('Plots/Windrose_MH_west_{}.png'.format(months[i]))
        
#%%
#Try to create histograms of wind speed to find location of main positions of flow (does not work well with data)
gr_south_mean_WD = np.zeros([12,15])*np.nan
gr_south_mean_WS = np.zeros([12,15])*np.nan
lats_south = np.arange(-50, -35, 1)

gr_north_mean_WD = np.zeros([12,15])*np.nan
gr_north_mean_WS = np.zeros([12,15])*np.nan
lats_north = np.arange(-35, -20, 1)

for i in range(12):
    
    if i+1 in gr_south.groups:
        for j in range(len(lats_south)):
            gr_south_mean_WD[i][j] = gr_south.get_group(i+1)[(gr_south.get_group(i+1).Latitude >= lats_south[j]) & (gr_south.get_group(i+1).Latitude < (lats_south[j] + 1))].mean()['Wind Direction']
            gr_south_mean_WS[i][j] = gr_south.get_group(i+1)[(gr_south.get_group(i+1).Latitude >= lats_south[j]) & (gr_south.get_group(i+1).Latitude < (lats_south[j] + 1))].mean()['Wind Speed/Force']

    if i+1 in gr_north.groups:
        for j in range(len(lats_north)):
            gr_north_mean_WD[i][j] = gr_north.get_group(i+1)[(gr_north.get_group(i+1).Latitude >= lats_north[j]) & (gr_north.get_group(i+1).Latitude < (lats_north[j] + 1))].mean()['Wind Direction']
            gr_north_mean_WS[i][j] = gr_north.get_group(i+1)[(gr_north.get_group(i+1).Latitude >= lats_north[j]) & (gr_north.get_group(i+1).Latitude < (lats_north[j] + 1))].mean()['Wind Speed/Force']

# #Plot histograms of data
# for i in range(12):
#     plt.figure()
#     plt.title('Southern Box Mascerene High ' + months[i])   
#     plt.bar(lats_south, gr_south_mean_WS[i], align = 'edge') 
#     plt.xlim(lats_south[0], lats_south[-1])
#     plt.ylim(0,10)
#     plt.grid()
#     plt.xlabel('Latitude')
#     plt.ylabel('Wind Speed (Beaufort)')
    
#     plt.figure()
#     plt.title('Northern Box Mascerene High ' + months[i])   
#     plt.bar(lats_north, gr_north_mean_WS[i], align = 'edge') 
#     plt.xlim(lats_north[0], lats_north[-1])
#     plt.ylim(0,10)
#     plt.grid()
#     plt.xlabel('Latitude')
#     plt.ylabel('Wind Speed (Beaufort)')
    
#%%
#Repeat analysis for a certain decade instead od climatological cycle
start_year = 1840
end_year = start_year + 10
data_dp = data_ubl[(data_ubl['Entry Date Time'].dt.year >= start_year) & (data_ubl['Entry Date Time'].dt.year <= end_year)]

#Split data into MH boxes
data_MH_south = data_dp.query('-50 <= Latitude <= -35 & 60 <= Longitude <= 105')
data_MH_north = data_dp.query('-35 <= Latitude <= -20 & 60 <= Longitude <= 105')
data_MH_west = data_dp.query('-50 <= Latitude <= -20 & 50 <= Longitude <= 60')
data_MH_east = data_dp.query('-50 <= Latitude <= -20 & 105 <= Longitude <= 115')
  
#Create worldmap and add MH boxes
worldmap_dec = WorldMap(fig_size = (30, 15), real_color = False, zoom = [20, 130, -55, -0])
plt.title('Decade {}-{} '.format(start_year, end_year), fontsize = 20)
w = 45
h = 15
s = 58
rect_south=mpatches.Rectangle((s,-50), w, h, fill = False, color = "red", linewidth = 3)
rect_north=mpatches.Rectangle((s,-35), w, h, fill = False, color = "red", linewidth = 3)
rect_east=mpatches.Rectangle((s+w,-50), 8, 2*h, fill = False, color = "red", linewidth = 3)
rect_west=mpatches.Rectangle((s-8,-50), 8, 2*h, fill = False, color = "red", linewidth = 3)    
worldmap_dec.Ax.add_patch(rect_south)
worldmap_dec.Ax.add_patch(rect_north)
worldmap_dec.Ax.add_patch(rect_east)
worldmap_dec.Ax.add_patch(rect_west)

#Create quiver plot for all data points in decade the MH region
quiv = worldmap_dec.Ax.quiver(data_dp.Longitude, data_dp.Latitude, -data_dp['Wind Speed/Force']*np.sin(data_dp['Wind Direction'].to_numpy().astype(float)*2*np.pi/360),
                             -data_dp['Wind Speed/Force']*np.cos(data_dp['Wind Direction'].to_numpy().astype(float)*2*np.pi/360),
                             scale = 300)

#Create windrose for each box
ax_south = WindroseAxes.from_ax()
plt.title('Mascarene High Southern Box Decade {}-{} '.format(start_year, end_year), fontsize = 20)
ax_south.bar(data_MH_south['Wind Direction'], data_MH_south['Wind Speed/Force'], normed=True, opening=0.8, edgecolor='white')
ax_south.set_legend()
fig_S = plt.gcf()
        
ax_north = WindroseAxes.from_ax()
plt.title('Mascarene High northern Box Decade {}-{} '.format(start_year, end_year), fontsize = 20)
ax_north.bar(data_MH_north['Wind Direction'], data_MH_north['Wind Speed/Force'], normed=True, opening=0.8, edgecolor='white')
ax_north.set_legend()
fig_N = plt.gcf()

ax_east = WindroseAxes.from_ax()
plt.title('Mascarene High eastern Box Decade {}-{} '.format(start_year, end_year), fontsize = 20)
ax_east.bar(data_MH_east['Wind Direction'], data_MH_east['Wind Speed/Force'], normed=True, opening=0.8, edgecolor='white')
ax_east.set_legend()
fig_E = plt.gcf()

ax_west = WindroseAxes.from_ax()
plt.title('Mascarene High western Box Decade {}-{} '.format(start_year, end_year), fontsize = 20)
ax_west.bar(data_MH_west['Wind Direction'], data_MH_west['Wind Speed/Force'], normed=True, opening=0.8, edgecolor='white')
ax_west.set_legend()
fig.W = plt.gcf()

#Save Figures
worldmap_dec.Fig.savefig('Plots/Windvectors_MH_{}_{}.png'.format(start_year, end_year))
fig_S.savefig('Plots/Windrose_MH_south_{}_{}.png'.format(start_year, end_year))
fig_N.savefig('Plots/Windrose_MH_north_{}_{}.png'.format(start_year, end_year))
fig_E.savefig('Plots/Windrose_MH_east_{}_{}.png'.format(start_year, end_year))
fig_W.savefig('Plots/Windrose_MH_west_{}_{}.png'.format(start_year, end_year))