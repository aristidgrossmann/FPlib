############################   IMPORTS   #######################
import numpy as np
import os
import sys
from matplotlib import pyplot as plt
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[3]))

from fplib import fpplot, fpfit
from fplib.models.AbsCosine import AbsCosine



############################   DATA LOADING   #######################
os.chdir('examples/fits/abs-cosine/data')

steel_single_line_velo_calibration = np.genfromtxt('kalibration_stahl_2_komma_25_mms.txt')

os.chdir('..')





############################   RAW DATA PLOTTING   #######################
#SPECIFICATIONS
ydata = steel_single_line_velo_calibration 
xdata = np.arange(1, len(ydata)+1)

title = 'Steel: Velocity calibration'
xlabel = 'Channel [-]'
ylabel = 'Number of fringes [-]'
legend_loc = 'lower center'
file_name = "rawData"

fpplot.plot_raw_data_no_uncertainty(xdata=xdata, ydata=ydata, title=title, xlabel=xlabel, 
                                     ylabel=ylabel, legend_loc=legend_loc, file_name=file_name)



############################   MODEL FITTING   #######################
#SPECIFICATIONS  
ydata = steel_single_line_velo_calibration
xdata = np.arange(1, len(ydata)+1)
xerr = np.ones_like(xdata)*1/np.sqrt(12)
yerr = ydata*np.sqrt((0.5e-2)**2 + (3.2e-6)**2 + (4e-4)**2) 

model = AbsCosine  #name of the model function 
p0 =  [875.5, 0.006128, 0.06646]  #starting guess

optimize_starting_guess = True   

fit1_label = 'cos fit'

plot_uncertainties = False
exclude_zero_count_data_points = True
compressed_Latex_output = True

xlims = (0, 1000)


xlabels = 'Channel [-]'
plot1_ylabel = 'Number of fringes [-]'
plot2_ylabel = 'Residuum [-]'

plot1_title = 'Steel: Velocity calibration: curve fit'
plot2_title = 'Residuals of cos fit'

plot1_legend_loc = 'lower center'
plot2_legend_loc = 'lower left'

file_name = "velocity_calibration_fit"


popt_cal, popt_cal_std, pcorr = fpfit.general_curve_fit(xdata=xdata, ydata=ydata, xerr=xerr, yerr=yerr, model=model, p0=p0, 
                                                optimize_starting_guess=optimize_starting_guess,  
                                                fit1_label=fit1_label, plot_uncertainties=plot_uncertainties, xlabels=xlabels, 
                                                plot1_ylabel=plot1_ylabel, plot2_ylabel=plot2_ylabel, plot1_title=plot1_title,  
                                                plot1_legend_loc=plot1_legend_loc, file_name=file_name, 
                                                exclude_zero_count_data_points=exclude_zero_count_data_points,
                                                compressed_Latex_output=compressed_Latex_output, xlims=xlims)
