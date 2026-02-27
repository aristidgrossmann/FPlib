############################   IMPORTS   #######################
import numpy as np
import os
import sys
from matplotlib import pyplot as plt
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[3]))
from src import fp_library



############################   DATA LOADING   #######################
os.chdir('examples/fits/simple-Gaussian/data')

data_no_absorber = np.genfromtxt('semi 35,70cm.txt')
channels = np.arange(0, len(data_no_absorber))

os.chdir('..')



############################   RAW DATA PLOTTING   #######################

#########  SPECIFICATIONS   ##########
xdata = channels
ydata = data_no_absorber 

title = r'Peak, $35.70$cm'
xlabel = 'Channel [-]'
ylabel = 'Counts [-]'
legend_loc = 'upper left'
file_name = "rawData"

zoom_xlim = (1450, 1720)  # Zoomed x range
zoom_ylim = (20, 160)  # Zoomed y range

# Zoom window size and position
scaling_factor = 0   #scales the size of the zoomed window
zoom_window_position = "center right"  # Position of the zoom window (not needed here because scaling_factor = 0)

custom_plot = True

(ax, insetax) = fp_library.plot_raw_data_no_uncertainty_zoom_window(xdata=xdata, ydata=ydata, title=title, xlabel=xlabel, ylabel=ylabel, 
            legend_loc=legend_loc, file_name=file_name, zoom_xlim=zoom_xlim, zoom_ylim=zoom_ylim, 
            scaling_factor=scaling_factor, zoom_window_position=zoom_window_position, custom_plot= custom_plot)

# add some custom text in figure (not supported by the library function)
ax.text(1400, 140, "Fit Gaussian distrib to this peak", ha='right') 



############################   FIT SIMPLE GAUSSIAN   #######################
channel_low = 1450
channel_high = 1720

xdata = channels
ydata = data_no_absorber

#########  SPECIFICATIONS   ##########

xdata = channels[channel_low:channel_high]
ydata = data_no_absorber[channel_low:channel_high]
xerr = np.ones_like(channels[channel_low:channel_high])/np.sqrt(12)
yerr = np.sqrt(data_no_absorber[channel_low:channel_high])

model = fp_library.gaussian_fit  #name of the model function (supported: linear_fit, exponential_fit, inverse_exponential_fit, gaussian_fit, double_gaussian_fit)
p0 = [130,1600, 100]  #starting guess

#if true, a curvefit without uncertainties is performed. the result is taken as the new starting guess. Aims at improving convergence
optimize_starting_guess = True  

plot_uncertainties = False

fit1_label = 'Gaussian fit'

xlabels = 'Channel [-]'

plot1_ylabel = 'Counts [-]'
plot2_ylabel = 'Residuum [-]'

plot1_title = 'Peak 2: Gaussian fit'
plot2_title = 'Residuals of Gaussian fit'

plot1_legend_loc = 'upper right'
plot2_legend_loc = 'lower right'

file_name = "Peak2_singlegauss_fit"


popt, popt_std, pcorr = fp_library.general_curve_fit(xdata=xdata, ydata=ydata, xerr=xerr, yerr=yerr, model=model, p0=p0, 
                                                    optimize_starting_guess=optimize_starting_guess,  
                                                    fit1_label=fit1_label, plot_uncertainties=plot_uncertainties, xlabels=xlabels, 
                                                    plot1_ylabel=plot1_ylabel, plot2_ylabel=plot2_ylabel, plot1_title=plot1_title,  
                                                    plot2_title=plot2_title, plot1_legend_loc=plot1_legend_loc, plot2_legend_loc=plot2_legend_loc, 
                                                    file_name=file_name)

peak_2_channel_singlegauss = np.array([popt[1], popt_std[1]])   #read out optimized parameters


