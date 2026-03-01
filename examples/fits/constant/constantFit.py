############################   IMPORTS   #######################
import numpy as np
import os
import sys
from matplotlib import pyplot as plt
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[3]))

from fplib import fpplot, fpfit
from fplib.models.Constant import Constant



############################   DATA LOADING   #######################
os.chdir('examples/fits/constant/data')

x_exemplary = np.array([1,2,3,4,5,6,7])
y_exemplary = np.array([1,2,1,1,1.4, 2, 1])

os.chdir('..')




############################   PLOT RAW DATA   #######################################
#########  SPECIFICATIONS   ##########
xdata = x_exemplary
ydata = y_exemplary 
xerr = np.ones_like(xdata)/10
yerr = xerr/2

title = 'exemplary constant fit'
xlabel = 'x [-]'
ylabel = 'y [-]'
legend_loc = 'upper center'
file_name = "rawData"

fpplot.plot_raw_data_with_uncertainty(xdata=xdata, ydata=ydata, xerr=xerr, yerr=yerr, title=title, xlabel=xlabel, 
                                          ylabel=ylabel, legend_loc=legend_loc, file_name=file_name)






############################   FIT EXPONTENTIAL FUNCTION   #######################################
#########  SPECIFICATIONS   ##########
xdata = x_exemplary
ydata = y_exemplary 
xerr = np.ones_like(xdata)*np.sqrt(1/12)
yerr = xerr

model = Constant  #name of the model function 
p0 = [1]  #starting guess

#if true, a curvefit without uncertainties is performed. the result is taken as the new starting guess. Aims at improving convergence
optimize_starting_guess = True   

fit1_label = 'constant fit'

plot_uncertainties = True
compressed_Latex_output = True

xlabels = 'x [-]'
plot1_ylabel = 'y [-]'
plot2_ylabel = 'Residuum [-]'

plot1_title = 'exemplary constant fit'

plot1_legend_loc = 'upper center'

file_name = "constantFit"

popt, popt_std, pcorr = fpfit.general_curve_fit(xdata=xdata, ydata=ydata, xerr=xerr, yerr=yerr, model=model, p0=p0, 
                                                optimize_starting_guess=optimize_starting_guess,  
                                                fit1_label=fit1_label, plot_uncertainties=plot_uncertainties, xlabels=xlabels, 
                                                plot1_ylabel=plot1_ylabel, plot2_ylabel=plot2_ylabel, plot1_title=plot1_title,  
                                                plot1_legend_loc=plot1_legend_loc, file_name=file_name, 
                                                compressed_Latex_output = compressed_Latex_output)