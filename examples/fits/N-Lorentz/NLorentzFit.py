############################   IMPORTS   #######################
import numpy as np
import os
import sys
from matplotlib import pyplot as plt
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[3]))
from src import fp_library



############################   DATA LOADING   #######################
os.chdir('examples/fits/N-Lorentz/data')

background_noise_spectrum = np.genfromtxt('Rauschmessung.txt', skip_header = 0)/(15*60)
Z_without_absorber = np.genfromtxt('without_absorber_without_velocity.txt')/(30*60)

os.chdir('..')





############################   FIT SUM OF N LORENTZ DISTRIBUTIONS   #######################################

lb = 0
ub = 700

#########  SPECIFICATIONS   ##########
ydata = Z_without_absorber[0:len(background_noise_spectrum)] - background_noise_spectrum*2 
xdata = np.arange(1, len(ydata)+1)
ydata = ydata[lb:ub]
xdata = xdata[lb:ub]
xerr = np.ones_like(xdata)*1/np.sqrt(12)
yerr = np.sqrt(ydata)*1

model = fp_library.N_lorentz_distribution_fit #name of the model function (supported: linear_fit, exponential_fit, inverse_exponential_fit, gaussian_fit, double_gaussian_fit)
p0 = [ 1.89191629e+04, 1.93498813e+02,  3.37366915e+01, 
      7.12745398e+03, 4.44602681e+01, 1.38779250e+01,  
      1.39511247e+04,  4.30750490e+02, 7.15511361e+01, 
     2000, 620, 100]  #starting guess (4 Lorentz distibs)

#if true, a curvefit without uncertainties is performed. the result is taken as the new starting guess. Aims at improving convergence
optimize_starting_guess = True   

fit1_label = 'fit'

plot_uncertainties = False


xlabels = 'Channel [-]'
plot1_ylabel = 'Count rate [1/s]'
plot2_ylabel = 'Residuum [1/s]'

plot1_title = 'Z(without absorber): Lorentz fit'
plot2_title = 'Residuals of Lorentz fit'

plot1_legend_loc = 'upper right'
plot2_legend_loc = 'lower right'

file_name = "4LorentzFit"


popt_lorentz, popt_std_lorentz, pcorr = fp_library.general_curve_fit(xdata, ydata, xerr, yerr, model, p0, optimize_starting_guess,   # data and optimization settings
                                                    fit1_label, plot_uncertainties, xlabels, plot1_ylabel, plot2_ylabel, plot1_title,  #plot settings
                                                    plot2_title, plot1_legend_loc, plot2_legend_loc, file_name, exclude_0_err = True)
