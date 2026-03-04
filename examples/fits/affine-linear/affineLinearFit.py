############################   IMPORTS   #######################
import numpy as np
import os
import sys
from matplotlib import pyplot as plt
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[3]))

from fplib import fpplot, fpfit
from fplib.models.AffineLinear import AffineLinear



############################   DATA LOADING   #######################
os.chdir('examples/fits/affine-linear/data')

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






############################   FIT ax + b FUNCTION   #######################################
#########  SPECIFICATIONS   ##########
xdata = x_exemplary
ydata = y_exemplary 
xerr = np.ones_like(xdata)*np.sqrt(1/12)
yerr = xerr

model = AffineLinear  #name of the model function 
p0 = [0, 1]  #starting guess

optimize_starting_guess = True   

fit_label = 'linear fit'

plot_uncertainties = True
compressed_Latex_output = True

xlabels = 'x [-]'
curveFitPlot_ylabel = 'y [-]'
residualPlot_ylabel = 'Residuum [-]'

plot_title = 'affine linear fit'

legend_loc = 'upper left'

file_name = "affineLinearFit"

popt, popt_std, pcorr = fpfit.general_curve_fit(xdata=xdata, ydata=ydata, xerr=xerr, yerr=yerr, model=model, p0=p0, 
                                                optimize_starting_guess=optimize_starting_guess,  
                                                fit_label=fit_label, plot_uncertainties=plot_uncertainties, xlabels=xlabels, 
                                                curveFitPlot_ylabel=curveFitPlot_ylabel, residualPlot_ylabel=residualPlot_ylabel, plot_title=plot_title,  
                                                legend_loc=legend_loc, file_name=file_name, 
                                                compressed_Latex_output = compressed_Latex_output)