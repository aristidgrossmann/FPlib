############################   IMPORTS   #######################
import numpy as np
import os
import sys
from matplotlib import pyplot as plt
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[3]))

from fplib import fpplot, fpfit
from fplib.models.GaussianCCDF import GaussianCCDF




############################   DATA LOADING   #######################
os.chdir('examples/fits/Gaussian-CCDF/data')

data = np.genfromtxt('distance_voltage_measurements.txt', skip_header = 1)
voltage_mean_and_uncertainty = np.ones_like(data)
offset = -0.005
for i in range(len(data)):
    voltage_mean_and_uncertainty[i][0] = data[i][0]
    voltage_mean_and_uncertainty[i][1] = np.mean(data[i][1]) - offset
    voltage_mean_and_uncertainty[i][2] = np.sqrt(((data[i][2] - data[i][1])/2)**2 + (0.001/np.sqrt(12))**2)

position = voltage_mean_and_uncertainty[:, 0] 
position = np.ones_like(position)*position[0] - position #distance from grating
position_std = np.ones_like(position)*0.01/np.sqrt(12)

current = voltage_mean_and_uncertainty[:, 1]*100/10    #in mücro A
current_std = voltage_mean_and_uncertainty[:, 2]*100/10

os.chdir('..')



############################   PLOT THE DATA   #######################################

#SPECIFICATIONS
xdata = position
ydata = current 
xerr = position_std
yerr = current_std

# Plot labels and title
title = r'$\alpha-$radiation: ionization current'
xlabel = 'Distance [cm]'
ylabel = r'Current [$\mu$ A]'
legend_loc = 'upper right'
file_name = "rawData"

fpplot.plot_raw_data_with_uncertainty(xdata, ydata, xerr, yerr, title, xlabel, ylabel, legend_loc, file_name)




cutoff_214_Po = 10

#SPECIFICATIONS
xdata = position[cutoff_214_Po:]
ydata = current[cutoff_214_Po:]
xerr = position_std[cutoff_214_Po:]
yerr = current_std[cutoff_214_Po:]

# Plot labels and title
title = r'$\alpha-$radiation: ionization current'
xlabel = 'Distance [cm]'
ylabel = r'Current [$\mu$ A]'
legend_loc = 'upper right'
file_name = "PolloniumPeak"

fpplot.plot_raw_data_with_uncertainty(xdata, ydata, xerr, yerr, title, xlabel, ylabel, legend_loc, file_name)



############################   FIT GAUSSIAN CCDF   #######################################

#SPECIFICATIONS
xdata = position[cutoff_214_Po:]
ydata = current[cutoff_214_Po:]
xerr = position_std[cutoff_214_Po:]
yerr = current_std[cutoff_214_Po:]

model = GaussianCCDF  #name of the model function (supported: linear_fit, exponential_fit, inverse_exponential_fit, gaussian_fit, double_gaussian_fit)
p0 = [0.25, 3, 1]  #starting guess

#if true, a curvefit without uncertainties is performed. the result is taken as the new starting guess. Aims at improving convergence
optimize_starting_guess = False   

fit1_label = 'gaussian ccdf fit'

plot_uncertainties = True


xlabels = 'Distance [cm]'
plot1_ylabel = r'Current [$\mu$ A]'
plot2_ylabel = r'Residuum [$\mu$ A]'

plot1_title = r'$\alpha-$radiation: ionization current: gaussian ccdf fit'

plot1_legend_loc = 'lower left'

file_name = "current_Po214_gaussian_ccdf_fit.pdf"


popt, popt_std, pcorr = fpfit.general_curve_fit(xdata=xdata, ydata=ydata, xerr=xerr, yerr=yerr, model=model, p0=p0, 
                                                optimize_starting_guess=optimize_starting_guess,  
                                                fit1_label=fit1_label, plot_uncertainties=plot_uncertainties, xlabels=xlabels, 
                                                plot1_ylabel=plot1_ylabel, plot2_ylabel=plot2_ylabel, plot1_title=plot1_title,  
                                                plot1_legend_loc=plot1_legend_loc, file_name=file_name)



