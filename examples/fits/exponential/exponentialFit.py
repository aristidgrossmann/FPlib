############################   IMPORTS   #######################
import numpy as np
import os
import sys
from matplotlib import pyplot as plt
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[3]))
from fplib import fpplot, fpfit
from fplib.models.Exponential import Exponential



############################   DATA LOADING   #######################
os.chdir('examples/fits/exponential/data')

data_no_absorber = np.genfromtxt('no absorber.txt')
data_1mm = np.genfromtxt('1mm.txt')
data_2mm = np.genfromtxt('2mm.txt')
data_3mm = np.genfromtxt('3mm.txt')
data_5mm = np.genfromtxt('5mm.txt')
data_10mm = np.genfromtxt('10mm.txt')
data_20mm = np.genfromtxt('20mm.txt')

lower_channel = 250
upper_channel = 450
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





############################   PLOT RAW DATA   #######################################
#########  SPECIFICATIONS   ##########
xdata = thickness
ydata = noof_counts 
xerr = thickness_uncertainties
yerr = np.sqrt(noof_counts)

title = r'$\gamma-$radiation absorption in lead'
xlabel = 'Absorber thickness [mm]'
ylabel = 'Counts [-]'
legend_loc = 'upper right'
file_name = "rawCounts"

fpplot.plot_raw_data_with_uncertainty(xdata=xdata, ydata=ydata, xerr=xerr, yerr=yerr, title=title, xlabel=xlabel, 
                                          ylabel=ylabel, legend_loc=legend_loc, file_name=file_name)






############################   FIT EXPONTENTIAL FUNCTION   #######################################
#########  SPECIFICATIONS   ##########
xdata = thickness
ydata = noof_counts 
xerr = thickness_uncertainties
yerr = np.sqrt(noof_counts)

model = Exponential  #name of the model function (supported: linear_fit, exponential_fit, inverse_exponential_fit, gaussian_fit, double_gaussian_fit)
p0 = [5000, 0.2]  #starting guess

#if true, a curvefit without uncertainties is performed. the result is taken as the new starting guess. Aims at improving convergence
optimize_starting_guess = True   

fit1_label = 'exponential fit'

plot_uncertainties = True


xlabels = 'Absorber thickness [mm]'
plot1_ylabel = 'Counts [-]'
plot2_ylabel = 'Residuum [-]'

plot1_title = r'$\gamma-$radiation absorption in lead: exponential fit'

plot1_legend_loc = 'upper right'

file_name = "ExponentialFit"

popt, popt_std, pcorr = fpfit.general_curve_fit(xdata=xdata, ydata=ydata, xerr=xerr, yerr=yerr, model=model, p0=p0, 
                                                optimize_starting_guess=optimize_starting_guess,  
                                                fit1_label=fit1_label, plot_uncertainties=plot_uncertainties, xlabels=xlabels, 
                                                plot1_ylabel=plot1_ylabel, plot2_ylabel=plot2_ylabel, plot1_title=plot1_title,  
                                                plot1_legend_loc=plot1_legend_loc, file_name=file_name)