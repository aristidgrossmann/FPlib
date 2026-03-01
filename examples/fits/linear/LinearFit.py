############################   IMPORTS   #######################
import numpy as np
import os
import sys
from matplotlib import pyplot as plt
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[3]))

from fplib import fpplot, fpfit
from fplib.models.Linear import Linear



############################   DATA LOADING   #######################
os.chdir('examples/fits/linear/data')

x_exemplary = np.array([1,2,3,4,5,6,7,8])
y_exemplary = np.array([1,2,2.5, 3, 5.5, 6, 7, 8.2])

os.chdir('..')





############################   PLOT RAW DATA   #######################################
#########  SPECIFICATIONS   ##########
xdata = x_exemplary
ydata = y_exemplary 
xerr = np.ones_like(xdata)*np.sqrt(1/12)
yerr = xerr

title = 'exemplary affine linear fit'
xlabel = 'x [-]'
ylabel = 'y [-]'
legend_loc = 'upper left'
file_name = "rawData"

fpplot.plot_raw_data_with_uncertainty(xdata=xdata, ydata=ydata, xerr=xerr, yerr=yerr, title=title, xlabel=xlabel, 
                                          ylabel=ylabel, legend_loc=legend_loc, file_name=file_name)






############################   FIT EXPONTENTIAL FUNCTION   #######################################
#########  SPECIFICATIONS   ##########
xdata = x_exemplary
ydata = y_exemplary 
xerr = np.ones_like(xdata)*np.sqrt(1/12)
yerr = xerr

model = Linear  #name of the model function (supported: linear_fit, exponential_fit, inverse_exponential_fit, gaussian_fit, double_gaussian_fit)
p0 = [1]  #starting guess

#if true, a curvefit without uncertainties is performed. the result is taken as the new starting guess. Aims at improving convergence
optimize_starting_guess = True   

fit1_label = 'linear fit'

plot_uncertainties = True


xlabels = 'x [-]'
plot1_ylabel = 'y [-]'
plot2_ylabel = 'Residuum [-]'

plot1_title = 'Linear fit'

plot1_legend_loc = 'upper left'

file_name = "LinearFit"

popt, popt_std, pcorr = fpfit.general_curve_fit(xdata=xdata, ydata=ydata, xerr=xerr, yerr=yerr, model=model, p0=p0, 
                                                optimize_starting_guess=optimize_starting_guess,  
                                                fit1_label=fit1_label, plot_uncertainties=plot_uncertainties, xlabels=xlabels, 
                                                plot1_ylabel=plot1_ylabel, plot2_ylabel=plot2_ylabel, plot1_title=plot1_title,  
                                                plot1_legend_loc=plot1_legend_loc, file_name=file_name)