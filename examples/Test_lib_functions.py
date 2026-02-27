import src.fp_library as fp_library
import os
import numpy as np


os.chdir('C:/Users/arist/Desktop/WS24_25/Fortgeschrittenenpraktikum/Template Code')

os.chdir('Daten/2_2 Scintillator gamma Led')

data_no_absorber = np.genfromtxt('no absorber.txt')
data_1mm = np.genfromtxt('1mm.txt')
data_2mm = np.genfromtxt('2mm.txt')
data_3mm = np.genfromtxt('3mm.txt')
data_5mm = np.genfromtxt('5mm.txt')
data_10mm = np.genfromtxt('10mm.txt')
data_20mm = np.genfromtxt('20mm.txt')

lower_channel = 320
upper_channel = 410
channels = np.arange(lower_channel, upper_channel)

data_full_spectrum = data_no_absorber
data_no_absorber = data_no_absorber[lower_channel-1:upper_channel-1]

data_1mm = data_1mm[lower_channel-1:upper_channel-1]
data_2mm = data_2mm[lower_channel-1:upper_channel-1]
data_3mm = data_3mm[lower_channel-1:upper_channel-1]
data_5mm = data_5mm[lower_channel-1:upper_channel-1]
data_10mm = data_10mm[lower_channel-1:upper_channel-1]
data_20mm = data_20mm[lower_channel-1:upper_channel-1]



thickness = np.array([0,1,2,3,5,10,20])
thickness_uncertainties = np.ones_like(thickness)*1/np.sqrt(12)
noof_counts = np.array([np.sum(data_no_absorber), np.sum(data_1mm), np.sum(data_2mm), np.sum(data_3mm), \
                        np.sum(data_5mm), np.sum(data_10mm), np.sum(data_20mm)])

os.chdir('..')
os.chdir('..')
os.chdir('Figs')


#############################################################################  TEST: plot_raw_data_no_uncertainty
#########  SPECIFICATIONS   ##########
xdata = np.arange(0, len(data_full_spectrum))
ydata = data_full_spectrum 

title = r'$\gamma-$radiation absorption in lead: inverse exponential fit'
xlabel = 'Channel [-]'
ylabel = 'Counts [-]'
legend_loc = 'upper right'
file_name = "whole_spectrum.svg"

fp_library.plot_raw_data_no_uncertainty(xdata, ydata, title, xlabel, ylabel, legend_loc, file_name)
fp_library.plot_raw_data_no_uncertainty(xdata, ydata, title, xlabel, ylabel, legend_loc, file_name, xlims = (0, 500), ylims = (0, 1000))  #example of setting ylims



#############################################################################  TEST: plot_raw_data_no_uncertainty_zoom_window

#########  SPECIFICATIONS   ##########
xdata = np.arange(0, len(data_full_spectrum))
ydata = data_full_spectrum 

# Plot labels and title
title = r'$\gamma-$radiation absorption in lead'
xlabel = 'Absorber thickness [mm]'
ylabel = 'Counts [-]'
legend_loc = 'upper right'
file_name = "counts_over_thickness_raw_zoom.svg"

# Zoom limits
zoom_xlim = (300, 400)  # Zoomed x range
zoom_ylim = (-100, 250)  # Zoomed y range

# Cutout size and position
scaling_factor = 2    #scales the size of the zoomed window
zoom_window_position = "center right"  # Position of the zoomed cutout

fp_library.plot_raw_data_no_uncertainty_zoom_window(xdata, ydata, title, xlabel, ylabel, legend_loc, file_name, #data and plot
                                             zoom_xlim, zoom_ylim, scaling_factor, zoom_window_position)



#############################################################################  TEST: plot_raw_data_with_uncertainty
#########  SPECIFICATIONS   ##########
xdata = thickness
ydata = noof_counts 
xerr = thickness_uncertainties
yerr = np.zeros_like(thickness_uncertainties)

title = r'$\gamma-$radiation absorption in lead'
xlabel = 'Absorber thickness [mm]'
ylabel = 'Counts [-]'
legend_loc = 'upper right'
file_name = "counts_over_thickness_raw.svg"

fp_library.plot_raw_data_with_uncertainty(xdata, ydata, xerr, yerr, title, xlabel, ylabel, legend_loc, file_name)



#############################################################################  TEST: plot_raw_data_with_uncertainty_zoom_window
#########  SPECIFICATIONS   ##########
xdata = thickness
ydata = noof_counts 
xerr = thickness_uncertainties
yerr = np.zeros_like(thickness_uncertainties)

# Plot labels and title
title = r'$\gamma-$radiation absorption in lead'
xlabel = 'Absorber thickness [mm]'
ylabel = 'Counts [-]'
legend_loc = 'upper right'
file_name = "counts_over_thickness_raw_zoom.svg"

# Zoom limits
zoom_xlim = (0, 3)  # Zoomed x range
zoom_ylim = (3200, 4100)  # Zoomed y range

# Cutout size and position
scaling_factor = 2    #scales the size of the zoomed window
zoom_window_position = "center right"  # Position of the zoomed cutout

fp_library.plot_raw_data_with_uncertainty_zoom_window(xdata, ydata, xerr, yerr, title, xlabel, ylabel, legend_loc, file_name, #data and plot
                                             zoom_xlim, zoom_ylim, scaling_factor, zoom_window_position)





#############################################################################  TEST: general_curve_fit
#########  SPECIFICATIONS   ##########
xdata = noof_counts
ydata = thickness
xerr = np.ones_like(thickness_uncertainties)*100
yerr = thickness_uncertainties

model = fp_library.inverse_exponential_fit  #name of the model function (supported: linear_fit, exponential_fit, inverse_exponential_fit, gaussian_fit, double_gaussian_fit)
p0 = [1/5, 4000]  #starting guess

#if true, a curvefit without uncertainties is performed. the result is taken as the new starting guess. Aims at improving convergence
optimize_starting_guess = True   

fit1_label = 'Inverse exponential fit'

plot_uncertainties = True


xlabels = 'Counts [-]'
plot1_ylabel = 'Absorber thickness [mm]'
plot2_ylabel = 'Residuum [mm]'

plot1_title = r'$\gamma-$radiation absorption in lead: inverse exponential fit'
plot2_title = 'Residuals of inverse exponential fit'

plot1_legend_loc = 'upper right'
plot2_legend_loc = 'lower right'

file_name = "counts_residuum.svg"


popt, popt_std, pcorr = fp_library.general_curve_fit(xdata, ydata, xerr, yerr, model, p0, optimize_starting_guess,   # data and optimization settings
                                                    fit1_label, plot_uncertainties, xlabels, plot1_ylabel, plot2_ylabel, plot1_title,  #plot settings
                                                    plot2_title, plot1_legend_loc, plot2_legend_loc, file_name)

print(popt, popt_std, pcorr)