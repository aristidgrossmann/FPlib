############################   IMPORTS   #######################
import numpy as np
import os
import sys
from matplotlib import pyplot as plt
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[3]))
from src import fp_library



############################   DATA LOADING   #######################
os.chdir('examples/fits/double-Gaussian/data')

anorganic_data = np.genfromtxt('Anorganic Szin.txt')
channels = np.arange(1, len(anorganic_data)+1)

os.chdir('..')



############################   RAW DATA PLOTTING   #######################

#########  SPECIFICATIONS   ##########
xdata = channels[1000:2000]
ydata = anorganic_data[1000:2000]

title = r'$\gamma-$spectrum: Photo peak; anorganic scintillator'
xlabel = 'Channel [-]'
ylabel = 'Counts [-]'
legend_loc = 'upper right'
file_name = "PhotoPeakRawData"

fp_library.plot_raw_data_no_uncertainty(xdata=xdata, ydata=ydata, title=title, xlabel=xlabel, ylabel=ylabel, 
            legend_loc=legend_loc, file_name=file_name)





############################   FIT DOUBLE GAUSSIAN (MAIN PEAK + BACKGROUND PEAK) #######################
#########  SPECIFICATIONS   ##########
xdata = channels[1000:2000]
ydata = anorganic_data[1000:2000]

xerr = np.ones_like(xdata)/np.sqrt(12)
yerr = np.sqrt(ydata)

model = fp_library.double_gaussian_fit  #name of the model function
p0 = [1300, 1300,  184.84817149, 500, 800, 300] #starting guess

#if true, a curvefit without uncertainties is performed. the result is taken as the new starting guess. Aims at improving convergence
optimize_starting_guess = True 

fit1_label = 'Double Gauss'

plot_uncertainties = False

xlabels = 'Channel [-]'
plot1_ylabel = 'Counts [-]'
plot2_ylabel = 'Residuum [-]'

plot1_title = r'$\gamma-$spectrum: Photo peak; anorganic scintillator: Double Gauss fit'
plot2_title = 'Residuals of double Gauss fit'

plot1_legend_loc = 'upper right'
plot2_legend_loc = 'lower right'

file_name = "mainAndBackgroundPeakFit"

peak_index = 0  #set to either 0 or 1 (to change which peak is the main and which is the background)
# if peak index is not set, it will lable them that the higher peak is the main peak and the smaller one is the background

popt, popt_std, popt_corr = fp_library.general_curve_fit(xdata=xdata, ydata=ydata, 
                        xerr=xerr, yerr=yerr, model=model, p0=p0, optimize_starting_guess=optimize_starting_guess,  
                     fit1_label=fit1_label, plot_uncertainties=plot_uncertainties, xlabels=xlabels, plot1_ylabel=plot1_ylabel, 
                     plot2_ylabel=plot2_ylabel, plot1_title=plot1_title, 
                     plot2_title=plot2_title, plot1_legend_loc=plot1_legend_loc, plot2_legend_loc=plot2_legend_loc, 
                     file_name=file_name, peak_index=peak_index)





### ALTERNATIVELY: specify peak labels and pass them to the function

peak1_label = 'Peak 1'
peak2_label = 'Peak 2'

file_name = "Peak1and2"

popt, popt_std, popt_corr = fp_library.general_curve_fit(xdata=xdata, ydata=ydata, 
                        xerr=xerr, yerr=yerr, model=model, p0=p0, optimize_starting_guess=optimize_starting_guess,  
                     fit1_label=fit1_label, plot_uncertainties=plot_uncertainties, xlabels=xlabels, plot1_ylabel=plot1_ylabel, 
                     plot2_ylabel=plot2_ylabel, plot1_title=plot1_title, 
                     plot2_title=plot2_title, plot1_legend_loc=plot1_legend_loc, plot2_legend_loc=plot2_legend_loc, 
                     file_name=file_name, peak1_label=peak1_label, peak2_label=peak2_label)



#### OR: simply dont specify anything. In this case, the higher peak will be the main and the smaller peak the background
file_name = "AutomaticMainAndBackground"

popt, popt_std, popt_corr = fp_library.general_curve_fit(xdata=xdata, ydata=ydata, 
                        xerr=xerr, yerr=yerr, model=model, p0=p0, optimize_starting_guess=optimize_starting_guess,  
                     fit1_label=fit1_label, plot_uncertainties=plot_uncertainties, xlabels=xlabels, plot1_ylabel=plot1_ylabel, 
                     plot2_ylabel=plot2_ylabel, plot1_title=plot1_title, 
                     plot2_title=plot2_title, plot1_legend_loc=plot1_legend_loc, plot2_legend_loc=plot2_legend_loc, 
                     file_name=file_name)



#### OR: if you dont wish them to be labeled at all, use the general N-gaussian-peaks model

model = fp_library.N_gaussian_fit  #name of the model function
file_name = "NGaussianFit-no-labels"

popt, popt_std, popt_corr = fp_library.general_curve_fit(xdata=xdata, ydata=ydata, 
                        xerr=xerr, yerr=yerr, model=model, p0=p0, optimize_starting_guess=optimize_starting_guess,  
                     fit1_label=fit1_label, plot_uncertainties=plot_uncertainties, xlabels=xlabels, plot1_ylabel=plot1_ylabel, 
                     plot2_ylabel=plot2_ylabel, plot1_title=plot1_title, 
                     plot2_title=plot2_title, plot1_legend_loc=plot1_legend_loc, plot2_legend_loc=plot2_legend_loc, 
                     file_name=file_name)



