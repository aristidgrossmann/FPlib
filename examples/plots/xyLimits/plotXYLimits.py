# region IMPORTS 
import numpy as np
import os
import sys
from matplotlib import pyplot as plt
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[3]))

from fplib import fpplot
# endregion



# region DATA LOADING
os.chdir('examples/plots/simple-no-uncertainties/data')

organic_data = np.genfromtxt('Organic Szin.txt')
channels = np.arange(1, len(organic_data)+1)


os.chdir('../..')
os.chdir('xyLimits')




##############################   PLOTTING WITH X & Y LIMITS   #############################

#note that using the plotting function with uncertainties and setting xerr, yerr = np.zeros(...) is also possible

#########  SPECIFICATIONS   ##########
xdata = channels
ydata = organic_data 

xlims = (2000, 4000)
ylims = (50, 200)

title = r'$\gamma-$spectrum, organic scintillator'
xlabel = 'Channel [-]'
ylabel = 'Counts [-]'
legend_loc = 'upper right'
file_name = "plotWithXYLims"



fpplot.plot_raw_data_no_uncertainty(xdata=xdata, ydata=ydata, title=title, xlabel=xlabel, ylabel=ylabel, 
                                        legend_loc=legend_loc, xlims=xlims, ylims=ylims, file_name=file_name)


