############################   IMPORTS   #######################
import numpy as np
import os
import sys
from matplotlib import pyplot as plt
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[3]))
from fplib import fpplot, fpfit
from fplib.models.NVoigt import NVoigt



############################   DATA LOADING   #######################
os.chdir('examples/fits/N-Voigt/data')

background_noise_spectrum = np.genfromtxt('Rauschmessung.txt', skip_header = 0)/(15*60)
Z_without_absorber = np.genfromtxt('without_absorber_without_velocity.txt')/(30*60)

os.chdir('..')





############################   FIT SUM OF N VOIGT PROFILES   #######################################

lb = 0
ub = 700

#########  SPECIFICATIONS   ##########
ydata = Z_without_absorber[0:len(background_noise_spectrum)] - background_noise_spectrum*2 
xdata = np.arange(1, len(ydata)+1)
ydata = ydata[lb:ub]
xdata = xdata[lb:ub]
xerr = np.ones_like(xdata)*1/np.sqrt(12)
yerr = np.sqrt(ydata)*1/np.sqrt(30*60)

model = NVoigt #name of the model function (supported: linear_fit, exponential_fit, inverse_exponential_fit, gaussian_fit, double_gaussian_fit)
p0 = [ 2e4, 1.93498813e+02, 10,  3.37366915e+01, 
      7e3, 4.44602681e+01, 10, 1.38779250e+01,  
      1e4,  4.30750490e+02, 10, 7.15511361e+01, 
     1e3, 600, 1, 20]  #starting guess

#if true, a curvefit without uncertainties is performed. the result is taken as the new starting guess. Aims at improving convergence
optimize_starting_guess = False   

fit1_label = 'fit'

plot_uncertainties = False
exclude_zero_count_data_points = True

xlabels = 'Channel [-]'
plot1_ylabel = 'Count rate [1/s]'
plot2_ylabel = 'Residuum [1/s]'

plot1_title = 'Z(without absorber): Voigt profile fit'

plot1_legend_loc = 'upper right'

file_name = "4VoigtsFit"


popt_voigt, popt_std_voigt, pcorr = fpfit.general_curve_fit(xdata=xdata, ydata=ydata, xerr=xerr, yerr=yerr, model=model, p0=p0, 
                                                optimize_starting_guess=optimize_starting_guess,  
                                                fit1_label=fit1_label, plot_uncertainties=plot_uncertainties, xlabels=xlabels, 
                                                plot1_ylabel=plot1_ylabel, plot2_ylabel=plot2_ylabel, plot1_title=plot1_title,  
                                                plot1_legend_loc=plot1_legend_loc, 
                                                file_name=file_name, exclude_zero_count_data_points=exclude_zero_count_data_points)
