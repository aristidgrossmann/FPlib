# region IMPORTS 
import numpy as np
import os
import sys
from matplotlib import pyplot as plt
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[3]))
from fplib import fpplot
# endregion



# region DATA LOADING
os.chdir('examples/plots/with-zoom-window/data')

data_no_absorber = np.genfromtxt('no absorber.txt')

os.chdir('..')
# endregion



############################   CLOSEUP ZOOM WINDOW  #################################

#########  SPECIFICATIONS   ##########
xdata = np.arange(0, len(data_no_absorber))
ydata = data_no_absorber 

title = r'$\gamma-$radiation absorption in lead: No absorber'
xlabel = 'Absorber thickness [mm]'
ylabel = 'Counts [-]'
legend_loc = 'upper right'
file_name = "zoomWindow"

zoom_xlim = (320, 400)  # Zoomed x range
zoom_ylim = (-100, 250)  # Zoomed y range

# Zoom window size and position
scaling_factor = 2    #scales the size of the zoomed window
zoom_window_position = "center right"  # Position of the zoom window

fpplot.plot_raw_data_no_uncertainty_zoom_window(xdata=xdata, ydata=ydata, title=title, xlabel=xlabel, ylabel=ylabel, 
            legend_loc=legend_loc, file_name=file_name, zoom_xlim=zoom_xlim, zoom_ylim=zoom_ylim, 
            scaling_factor=scaling_factor, zoom_window_position=zoom_window_position)




############################   USING THE ZOOM WINDOW AS A MARKER  #################################

#########  SPECIFICATIONS   ##########
xdata = np.arange(0, len(data_no_absorber))
ydata = data_no_absorber 

title = r'$\gamma-$radiation absorption in lead: No absorber'
xlabel = 'Absorber thickness [mm]'
ylabel = 'Counts [-]'
legend_loc = 'upper right'
file_name = "BoundingBox"

zoom_xlim = (320, 400)  # Zoomed x range
zoom_ylim = (-100, 250)  # Zoomed y range

# Zoom window size and position
scaling_factor = 0   #scales the size of the zoomed window
zoom_window_position = "center right"  # Position of the zoom window

fpplot.plot_raw_data_no_uncertainty_zoom_window(xdata=xdata, ydata=ydata, title=title, xlabel=xlabel, ylabel=ylabel, 
            legend_loc=legend_loc, file_name=file_name, zoom_xlim=zoom_xlim, zoom_ylim=zoom_ylim, 
            scaling_factor=scaling_factor, zoom_window_position=zoom_window_position)


