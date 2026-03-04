############################   IMPORTS   #######################
import numpy as np
import os
import sys
from matplotlib import pyplot as plt
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[3]))

from fplib import fpplot, fpfit
from fplib.models.DampedOscillation import DampedOscillation



############################   DATA CREATION   #######################
os.chdir('examples/fits/damped-oscillation')

np.random.seed(42)

t = np.linspace(0, 10, 1000)

# Damped oscillation parameters
A = 1.0          # initial amplitude
tau = 2      # damping coefficient
omega = 2.0*np.pi  # angular frequency (1 Hz)
phi = 0.0        # phase

y_clean = A * np.exp(- t/tau) * np.cos(omega * t + phi)

# Add Gaussian noise
noise_std = 0.1
noise = np.random.normal(0, noise_std, size=t.shape)
y_noisy = y_clean + noise





############################   FIT DAMPED OSCILLATION   #######################################
#########  SPECIFICATIONS   ##########
xdata = t
ydata = y_noisy 
xerr = np.zeros_like(xdata)
yerr = np.ones_like(ydata)*0.1

model = DampedOscillation  #name of the model function 
p0 = [1, 0.3, 6, 0]  #starting guess

optimize_starting_guess = True   

fit_label = 'fit'

plot_uncertainties = False
compressed_Latex_output = True

xlabels = 'x [-]'
curveFitPlot_ylabel = 'y [-]'
residualPlot_ylabel = 'Residuum [-]'

plot_title = 'Damped Oscillation Fit'

legend_loc = 'upper right'

file_name = "DampedOscillationFit"

popt, popt_std, pcorr = fpfit.general_curve_fit(xdata=xdata, ydata=ydata, xerr=xerr, yerr=yerr, model=model, p0=p0, 
                                                optimize_starting_guess=optimize_starting_guess,  
                                                fit_label=fit_label, plot_uncertainties=plot_uncertainties, xlabels=xlabels, 
                                                curveFitPlot_ylabel=curveFitPlot_ylabel, residualPlot_ylabel=residualPlot_ylabel, plot_title=plot_title,  
                                                legend_loc=legend_loc, file_name=file_name, 
                                                compressed_Latex_output = compressed_Latex_output)








