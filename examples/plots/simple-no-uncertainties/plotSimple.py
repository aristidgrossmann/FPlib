# region IMPORTS 
import numpy as np
import os
import sys
from matplotlib import pyplot as plt
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[3]))
from src import fp_library
# endregion



# region DATA LOADING
os.chdir('examples/plots/simple-no-uncertainties/data')

organic_data = np.genfromtxt('Organic Szin.txt')
channels = np.arange(1, len(organic_data)+1)


os.chdir('..')




##############################   PLOTTING WITHOUT ERRORBARS   #############################

#note that using the plotting function with uncertainties and setting xerr, yerr = np.zeros(...) is also possible

#########  SPECIFICATIONS   ##########
xdata = channels
ydata = organic_data 

title = r'$\gamma-$spectrum, organic scintillator'
xlabel = 'Channel [-]'
ylabel = 'Counts [-]'
legend_loc = 'upper right'
file_name = "simplePlot"



fp_library.plot_raw_data_no_uncertainty(xdata=xdata, ydata=ydata, title=title, xlabel=xlabel, ylabel=ylabel, 
                                        legend_loc=legend_loc, file_name=file_name)


