from matplotlib import pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.patches import Rectangle

from fplib.utils import utils_for_fpplot


####################    PLOTTERS    ####################
#This file contains plotting functions. supported are: 
#
# plot_raw_data_no_uncertainty()
# plot_raw_data_no_uncertainty_zoom_window()
# plot_raw_data_with_uncertainty()
# plot_raw_data_with_uncertainty_zoom_window()



# These 4 functions are also unified into the two functions (though i never used them); 
# they can be found at the bottom (commented out)
# plot_general_raw_data()
# plot_general_raw_data_zoom_window()





def plot_raw_data_no_uncertainty(xdata, ydata, title, xlabel, ylabel, legend_loc, file_name, custom_plot = False, xlims = None, ylims = None):


    fig, ax = plt.subplots(figsize = (6, 4.5))
    ax.scatter(xdata, ydata, s=8, label = 'Data', color = '#1f77b4')
    ax.set_xlabel(xlabel, fontsize = 12)
    ax.set_ylabel(ylabel, fontsize = 12)
    ax.tick_params(axis='both', labelsize=12)
    ax.set_title(title, fontsize = 13)
    ax.legend(loc = legend_loc, fontsize = 12)
    ax.grid(True)

    if xlims is not None:
        ax.set_xlim(xlims)
    if ylims is not None:
        ax.set_ylim(ylims)

    file_name += '.pdf'

    if custom_plot == False:
        plt.savefig(file_name, format='pdf', bbox_inches='tight')

        plt.show()
    else:
        return ax


def plot_raw_data_no_uncertainty_zoom_window(xdata, ydata, title, xlabel, ylabel, legend_loc, file_name, #data and plot
                                             zoom_xlim, zoom_ylim, scaling_factor, zoom_window_position, custom_plot = False, xlims = None, ylims = None):

    fig, ax = plt.subplots(figsize=(6, 4.5))

    # Plot the main data with error bars
    ax.scatter(xdata, ydata, s=8, label = 'Data', color = '#1f77b4')
    ax.set_xlabel(xlabel, fontsize = 12)
    ax.set_ylabel(ylabel, fontsize = 12)
    ax.tick_params(axis='both', labelsize=12)
    ax.set_title(title, fontsize = 13)
    ax.legend(loc = legend_loc, fontsize = 12)
    ax.grid(True)

    if xlims is not None:
        ax.set_xlim(xlims)
    if ylims is not None:
        ax.set_ylim(ylims)

    # Draw a dotted rectangle to indicate the zoomed area
    rect_width = zoom_xlim[1] - zoom_xlim[0]
    rect_height = zoom_ylim[1] - zoom_ylim[0]

    rect = Rectangle(
        (zoom_xlim[0], zoom_ylim[0]),  # Bottom-left corner
        rect_width,                    # Width
        rect_height,                   # Height
        edgecolor = 'red',
        facecolor='none', # Transparent fill
        linestyle=':',    # Dotted border
        linewidth=2       # Border thickness
    )
    ax.add_patch(rect)

    if scaling_factor > 0:
        # Add a zoomed inset axes
        x_range = max(xdata)- min(xdata)
        y_range = max(ydata) - min(ydata)
        inset_width = rect_width/x_range*4.0/3.0
        inset_height = rect_height/y_range
        max_length = max(inset_width, inset_height)
        inset_width *= scaling_factor/max_length
        inset_height *= scaling_factor/max_length

        inset_ax = inset_axes(ax, width=inset_width, height=inset_height, loc=zoom_window_position)

        # Plot the same data in the inset axes but zoomed
        inset_ax.scatter(xdata, ydata, s=8, label = 'Data', color = '#1f77b4')
        inset_ax.set_xlim(zoom_xlim)
        inset_ax.set_ylim(zoom_ylim)

        # Remove tick labels from the inset for clarity
        inset_ax.set_xticklabels([])
        inset_ax.set_yticklabels([])

        # Ensure the inset stays within the bounds of the figure
        plt.tight_layout()

                ### draw connecting line 

        # Convert the zoom window position to real coordinates
        zoom_window_point = utils_for_fpplot.position_to_coords(zoom_window_position, (min(xdata), max(xdata)), (min(ydata), max(ydata)))  # Convert the zoom window position to real coordinates

        # Find the nearest corner of the rectangle to the zoom window point
        closest_point = utils_for_fpplot.closest_point_on_rectangle(zoom_window_point, rect)

        # Draw a dotted line from the nearest corner to the zoom window position
        ax.plot([closest_point[0], zoom_window_point[0]], [closest_point[1], zoom_window_point[1]], '--', color = 'red')

    else: 
        inset_ax = None




    file_name += '.pdf'

    if custom_plot == False:
        plt.savefig(file_name, format='pdf', bbox_inches='tight')
        plt.show()
    else:
        return ax, inset_ax
    




def plot_raw_data_with_uncertainty(xdata, ydata, xerr, yerr, title, xlabel, ylabel, legend_loc, file_name, custom_plot = False, xlims = None, ylims = None):
    fig, ax = plt.subplots(figsize = (6, 4.5))

    if(np.all(xerr == 0)):
        ax.errorbar(xdata, ydata, yerr = yerr, fmt='.',elinewidth = 1.5,  markersize = 10, capsize=0, label = 'Data', color = '#1f77b4')
    elif(np.all(yerr == 0)):
        ax.errorbar(xdata, ydata, xerr = xerr, fmt='.',elinewidth = 1.5,  markersize = 10, capsize=0, label = 'Data', color = '#1f77b4')
    else:
        ax.errorbar(xdata, ydata, xerr = xerr, yerr = yerr, fmt='.',elinewidth = 1.5,  markersize = 10, capsize=0, label = 'Data', color = '#1f77b4')
    ax.set_xlabel(xlabel, fontsize = 12)
    ax.set_ylabel(ylabel, fontsize = 12)
    ax.tick_params(axis='both', labelsize=12)
    ax.set_title(title, fontsize = 13)
    ax.legend(loc = legend_loc, fontsize = 12)
    ax.grid(True)

    if xlims is not None:
        ax.set_xlim(xlims)
    if ylims is not None:
        ax.set_ylim(ylims)
    
    file_name += '.pdf'

    if custom_plot == False:
        plt.savefig(file_name, format='pdf', bbox_inches='tight')
        plt.show()
    else:
        return ax




def plot_raw_data_with_uncertainty_zoom_window(xdata, ydata, xerr, yerr, title, xlabel, ylabel, legend_loc, file_name, #data and plot
                                             zoom_xlim, zoom_ylim, scaling_factor, zoom_window_position, custom_plot = False, xlims = None, ylims = None):
    
    fig, ax = plt.subplots(figsize=(6, 4.5))

    # Plot the main data with error bars
    if(np.all(xerr == 0)):
        ax.errorbar(xdata, ydata, yerr=yerr, fmt='.', elinewidth=1.5,  markersize = 10, capsize=0, color = '#1f77b4', label = 'Data')
    elif(np.all(yerr == 0)):
        ax.errorbar(xdata, ydata, xerr=xerr, fmt='.', elinewidth=1.5,  markersize = 10, capsize=0, color = '#1f77b4', label = 'Data')
    else:
        ax.errorbar(xdata, ydata, xerr=xerr, yerr=yerr, fmt='.', elinewidth=1.5,  markersize = 10, capsize=0, color = '#1f77b4', label = 'Data')
    ax.set_xlabel(xlabel, fontsize = 12)
    ax.set_ylabel(ylabel, fontsize = 12)
    ax.tick_params(axis='both', labelsize=12)
    ax.set_title(title, fontsize = 13)
    ax.legend(loc = legend_loc, fontsize = 12)
    ax.grid(True)

    if xlims is not None:
        ax.set_xlim(xlims)
    if ylims is not None:
        ax.set_ylim(ylims)

    # Draw a dotted rectangle to indicate the zoomed area
    rect_width = zoom_xlim[1] - zoom_xlim[0]
    rect_height = zoom_ylim[1] - zoom_ylim[0]

    rect = Rectangle(
        (zoom_xlim[0], zoom_ylim[0]),  # Bottom-left corner
        rect_width,                    # Width
        rect_height,                   # Height
        edgecolor = 'red',
        facecolor='none', # Transparent fill
        linestyle=':',    # Dotted border
        linewidth=2       # Border thickness
    )
    ax.add_patch(rect)

    if scaling_factor > 0:
        # Add a zoomed inset axes
        x_range = max(xdata)- min(xdata)
        y_range = max(ydata) - min(ydata)
        inset_width = rect_width/x_range*4.0/3.0
        inset_height = rect_height/y_range
        max_length = max(inset_width, inset_height)
        inset_width *= scaling_factor/max_length
        inset_height *= scaling_factor/max_length

        inset_ax = inset_axes(ax, width=inset_width, height=inset_height, loc=zoom_window_position)

        # Plot the same data in the inset axes but zoomed
        inset_ax.errorbar(xdata, ydata, xerr=xerr, yerr=yerr, fmt='.', elinewidth=1,  markersize = 10, capsize=0, color = '#1f77b4')
        inset_ax.set_xlim(zoom_xlim)
        inset_ax.set_ylim(zoom_ylim)

        # Remove tick labels from the inset for clarity
        inset_ax.set_xticklabels([])
        inset_ax.set_yticklabels([])

        # Ensure the inset stays within the bounds of the figure
        plt.tight_layout()

        ### draw connecting line 

        # Convert the zoom window position to real coordinates
        zoom_window_point = utils_for_fpplot.position_to_coords(zoom_window_position, (min(xdata), max(xdata)), (min(ydata), max(ydata)))  # Convert the zoom window position to real coordinates

        # Find the nearest corner of the rectangle to the zoom window point
        closest_point = utils_for_fpplot.closest_point_on_rectangle(zoom_window_point, rect)

        # Draw a dotted line from the nearest corner to the zoom window position
        ax.plot([closest_point[0], zoom_window_point[0]], [closest_point[1], zoom_window_point[1]], '--', color = 'red')

    else: 
        inset_ax = None




    file_name += '.pdf'

    if custom_plot == False:
        plt.savefig(file_name, format='pdf', bbox_inches='tight')
        plt.show()
    else:
        return ax, inset_ax


'''

def plot_general_raw_data(xdata, ydata, xerr, yerr, title, xlabel, ylabel, legend_loc, file_name, plot_uncertainties = None,  xlims = None, ylims = None):
    if plot_uncertainties == False:
        plot_raw_data_no_uncertainty(xdata, ydata, title, xlabel, ylabel, legend_loc, file_name, xlims = None, ylims = None)
    else:
        plot_raw_data_with_uncertainty(xdata, ydata, xerr, yerr, title, xlabel, ylabel, legend_loc, file_name, xlims = None, ylims = None)
    return None

def plot_general_raw_data_zoom_window(xdata, ydata, xerr, yerr, title, xlabel, ylabel, legend_loc, file_name, #data and plot
                                             zoom_xlim, zoom_ylim, scaling_factor, zoom_window_position, plot_uncertainties = None, xlims = None, ylims = None):

    if plot_uncertainties == False: 
        plot_raw_data_no_uncertainty_zoom_window(xdata, ydata, title, xlabel, ylabel, legend_loc, file_name, #data and plot
                                             zoom_xlim, zoom_ylim, scaling_factor, zoom_window_position, xlims = None, ylims = None)
    else:
        plot_raw_data_with_uncertainty_zoom_window(xdata, ydata, xerr, yerr, title, xlabel, ylabel, legend_loc, file_name, #data and plot
                                             zoom_xlim, zoom_ylim, scaling_factor, zoom_window_position, xlims = None, ylims = None)
    return None

'''

