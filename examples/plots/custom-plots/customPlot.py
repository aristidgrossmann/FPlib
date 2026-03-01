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
os.chdir('examples/plots/custom-plots/data')

data_no_absorber = np.genfromtxt('no absorber.txt')

os.chdir('..')
# endregion





############################   CUSTOMIZE SIMPLE PLOT #################################

#########  SPECIFICATIONS   ##########
xdata = np.arange(0, len(data_no_absorber))
xerr = np.zeros_like(xdata)

ydata = data_no_absorber 
yerr = np.zeros_like(ydata)

title = r'$\gamma-$radiation absorption in lead: No absorber'
xlabel = 'Absorber thickness [mm]'
ylabel = 'Counts [-]'
legend_loc = 'upper right'
file_name = "customPlotSimple"


custom_plot = True  #returns the axes object of the plot

ax = fpplot.plot_raw_data_with_uncertainty(xdata=xdata, ydata=ydata, xerr=xerr, yerr=yerr, title=title, xlabel=xlabel, ylabel=ylabel, 
            legend_loc=legend_loc, file_name=file_name, custom_plot=custom_plot)

#add a horizontal line
ax.axhline(y=2000, xmin=0, xmax=1, color = 'black', linewidth = 3)

plt.savefig(file_name + '.pdf')

plt.show()




############################   CUSTOMIZE ZOOM WINDOW PLOT #################################

#########  SPECIFICATIONS   ##########
xdata = np.arange(0, len(data_no_absorber))
ydata = data_no_absorber 

title = r'$\gamma-$radiation absorption in lead: No absorber'
xlabel = 'Absorber thickness [mm]'
ylabel = 'Counts [-]'
legend_loc = 'upper right'
file_name = "customPlotZoomWindow"

zoom_xlim = (320, 400)  # Zoomed x range
zoom_ylim = (-100, 250)  # Zoomed y range

# Zoom window size and position
scaling_factor = 1    #scales the size of the zoomed window
zoom_window_position = "center right"  # Position of the zoom window

custom_plot = True  #returns all axes objects (in case of the zoom window, 2 objects)

(ax, insetax) = fpplot.plot_raw_data_no_uncertainty_zoom_window(xdata=xdata, ydata=ydata, title=title, xlabel=xlabel, ylabel=ylabel, 
            legend_loc=legend_loc, file_name=file_name, zoom_xlim=zoom_xlim, zoom_ylim=zoom_ylim, 
            scaling_factor=scaling_factor, zoom_window_position=zoom_window_position, custom_plot=custom_plot)

#add a horizontal line
ax.axhline(y=3500, xmin=0, xmax=1, color = 'black', linewidth = 3)

plt.savefig(file_name + '.pdf')

plt.show()





