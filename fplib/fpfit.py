import numpy as np
import scipy
from scipy.odr import ODR, Model, RealData
from matplotlib import pyplot as plt
from contextlib import redirect_stdout


from fplib.utils import utils_for_fpfit
from fplib.models.ModelTemplate import ModelTemplate
from fplib.models.Gaussian import Gaussian
from fplib.models.DoubleGaussian import DoubleGaussian
from fplib.models.LogLinear import LogLinear




def general_curve_fit(xdata, ydata, xerr, yerr, model:type[ModelTemplate], p0, optimize_starting_guess,  # data and optimization settings
                    fit1_label, plot_uncertainties, xlabels, plot1_ylabel, plot2_ylabel, plot1_title, plot1_legend_loc, file_name, 
                    xlims = None, ylims = None, #limits on the plot
                    custom_plot = False,   # if true, the plots axis object is returned 
                    peak_index = None, peak1_label = None, peak2_label = None,  #when fitting double gaussian
                    exclude_indices = None, #indices of the data points that are not accounted for when fitting
                    exclude_zero_count_data_points = False  # exclude 0 counts (they have variance 0 according to Poisson distrib)
                    ):
    
    ###################################################################################
    # perform chi2 curvefit (xerr = 0, yerr) or orthigonal distance regression (xerr, yerr)
    #get better starting values by fitting with (xerr = 0, yerr = 0)
    
    

    if np.any(exclude_indices):

        excluded_xdata = np.delete(xdata, exclude_indices)
        excluded_xerr = np.delete(xerr, exclude_indices)
        excluded_ydata = np.delete(ydata, exclude_indices)
        excluded_yerr = np.delete(yerr, exclude_indices)

        sorted_indices = np.argsort(excluded_xdata)

        excluded_xdata = excluded_xdata[sorted_indices]
        excluded_xerr = excluded_xerr[sorted_indices]
        excluded_ydata = excluded_ydata[sorted_indices]
        excluded_yerr = excluded_yerr[sorted_indices]


        xdata = xdata[exclude_indices]
        xerr = xerr[exclude_indices]
        ydata = ydata[exclude_indices]
        yerr = yerr[exclude_indices]
    



    if optimize_starting_guess:
        p0 = scipy.optimize.curve_fit(model.modelFunction, xdata= xdata, ydata = ydata,  p0 = p0)[0]
    if(np.all(xerr) == 0):
        algorithm = 'curve\\_fit'

        #does not account for data points with zero error 
        zero_error_data_points = []
        for i in range(len(xdata)):
            if yerr[i] == 0:
                zero_error_data_points.append(i)

        xdata = np.delete(xdata, zero_error_data_points)
        ydata = np.delete(ydata, zero_error_data_points)
        xerr = np.delete(xerr, zero_error_data_points)
        yerr = np.delete(yerr, zero_error_data_points)

        popt, pcov = scipy.optimize.curve_fit(model.modelFunction, xdata= xdata, ydata = ydata,  p0 = p0, \
                                           sigma = yerr, absolute_sigma = True)



    else: 

        if exclude_zero_count_data_points:
            zero_count_data_points = []
            for i in range(len(xdata)):
                if yerr[i] == 0:
                    zero_count_data_points.append(i)

            xdata = np.delete(xdata, zero_count_data_points)
            ydata = np.delete(ydata, zero_count_data_points)
            xerr = np.delete(xerr, zero_count_data_points)
            yerr = np.delete(yerr, zero_count_data_points)


        algorithm = 'ODR'

        data = RealData(xdata, ydata, sx=xerr, sy=yerr)

        # Create a Model object
        ODR_model = Model(utils_for_fpfit.curve_fit_to_odr_wrapper(model.modelFunction))

        # Create ODR object
        odr = ODR(data, ODR_model, beta0=p0)
        result = odr.run()
        popt, pcov = result.beta, result.cov_beta 

    #evaluate results
    popt_std = np.sqrt(np.diag(pcov))
    pcorr = pcov / np.outer(popt_std, popt_std)
    residuum = ydata - model.modelFunction(xdata, *popt)
    chi2_dof = utils_for_fpfit.chi_squared_dof(residuum, xerr*model.modelFunctionDerivative(xdata, *popt), yerr, len(residuum)-len(popt))

    #print results
    print('-----   curve fit: ', algorithm, '-------')
    print('Fitted Function:  ', model.getName())

    p0_str = '['
    for i in range(len(p0)):
        if (i == len(p0)-1):
            p0_str += f"{p0[i]:.4g}" + ']'
        else: 
            p0_str += f"{p0[i]:.4g}" + ', '
    print('Starting guess:  ', p0_str)

    # print parameters and uncertainties
    for i in range(len(popt)):
        print('parameter ', i+1, ': ', f"{popt[i]:.4g}", ' $\\pm$ ', f"{popt_std[i]:.4g}")

    #print correlation
    print('Correlation:')
    for i in range(len(popt)):
        corr_str = ''
        for j in range(len(popt)):
            if (j == len(popt)-1):
                corr_str += f"{pcorr[i][j]:.4g}"
            else:
                corr_str += f"{pcorr[i][j]:.4g}" + ' & '
        print(corr_str)

    #print chi2/dof
    print('chi2/dof:   ', f"{chi2_dof:.4g}")

    print('--------------------------')


    # Create the figure and subplots
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(4.5,4), sharex=True, gridspec_kw={'height_ratios': [2.2, 1]})

   
    xmodel = np.linspace(np.min(xdata), np.max(xdata), 1000)
    ymodel = model.modelFunction(xmodel, *popt)


    #########################################################    TOP PLOT
    #PLOT DATA
    if plot_uncertainties:  # WHEN PLOTTING UNCERTAINTIES
        if(np.all(xerr == 0)):
            ax1.errorbar(xdata, ydata, yerr=yerr, fmt='.', elinewidth=1.5,  markersize = 10, capsize=0, color = '#1f77b4', label = 'Data')
            if np.any(exclude_indices):
                ax1.errorbar(excluded_xdata, excluded_ydata, yerr=excluded_yerr, fmt='.', elinewidth=1.5,  markersize = 10, capsize=0, color = '#1f77b4')
        else:
            ax1.errorbar(xdata, ydata, xerr=xerr, yerr=yerr, fmt='.', elinewidth=1.5,  markersize = 10, capsize=0, color = '#1f77b4', label = 'Data')
            if np.any(exclude_indices):
                ax1.errorbar(excluded_xdata, excluded_ydata, xerr=excluded_xerr, yerr=excluded_yerr, fmt='.', elinewidth=1.5,  markersize = 10, capsize=0, color = '#1f77b4')


    else:
        ax1.scatter(xdata, ydata, s=8, color = '#1f77b4', label = 'Data')

        if np.any(exclude_indices):
            ax1.scatter(excluded_xdata, excluded_ydata, s=8, color = '#1f77b4')


    # PLOT FIT


    if(model.getName() == DoubleGaussian.getName()):  # if double gaussian, plot individual distributions and their sum

        #determine which gaussian is background gaussian:
        if peak1_label == None and peak2_label == None:
            if peak_index == None:  #background index not specified
                peak_index = np.argmin(np.array([popt[2], popt[5]])) #determine peak index via smaller variance
                if peak_index == 0:
                    ax1.plot(xmodel, Gaussian.modelFunction(xmodel, *popt[3:6]), color = '#006400', linewidth = 2, label = 'Background')
                    ax1.plot(xmodel, Gaussian.modelFunction(xmodel, *popt[0:3]), color = '#9467bd', linewidth = 2, label = 'Peak')
                    
                else:
                    ax1.plot(xmodel, Gaussian.modelFunction(xmodel, *popt[0:3]), color = '#006400', linewidth = 2, label = 'Background')
                    ax1.plot(xmodel, Gaussian.modelFunction(xmodel, *popt[3:6]), color = '#9467bd', linewidth = 2, label = 'Peak')


            else:  #peak index is specified
                if peak_index == 0:
                    ax1.plot(xmodel, Gaussian.modelFunction(xmodel, *popt[3:6]), color = '#006400', linewidth = 2, label = 'Background')
                    ax1.plot(xmodel, Gaussian.modelFunction(xmodel, *popt[0:3]), color = '#9467bd', linewidth = 2, label = 'Peak')
                    
                else:
                    ax1.plot(xmodel, Gaussian.modelFunction(xmodel, *popt[0:3]), color = '#006400', linewidth = 2, label = 'Background')
                    ax1.plot(xmodel, Gaussian.modelFunction(xmodel, *popt[3:6]), color = '#9467bd', linewidth = 2, label = 'Peak')
        else: 
            ax1.plot(xmodel, Gaussian.modelFunction(xmodel, *popt[0:3]), color = '#006400', linewidth = 2, label = peak1_label)
            ax1.plot(xmodel, Gaussian.modelFunction(xmodel, *popt[3:6]), color = '#9467bd', linewidth = 2, label = peak2_label)

        ax1.plot(xmodel, ymodel, color = 'red', linewidth = 2, label = 'Total')

    else: #for all other models, just plot one line
        ax1.plot(xmodel, ymodel, color = 'red', linewidth = 2, label = fit1_label)



    ax1.set_ylabel(plot1_ylabel, fontsize = 10)
    ax1.tick_params(axis='both', labelsize=10)
    ax1.legend(fontsize = 10, loc = plot1_legend_loc)
    ax1.set_title(plot1_title + r' ($\chi^2/N_{\mathrm{dof}}$ = ' + f"{chi2_dof:.4g}" + ')', fontsize = 11)
    ax1.grid(True)

    # if loglinear fit (exponential fit, but logarithmic y axis)
    ax1.set_xscale(model.getXscale())
    ax1.set_yscale(model.getYscale())


    if xlims is not None:
        ax1.set_xlim(xlims)
    if ylims is not None:
        ax1.set_ylim(ylims)



    ############################# Bottom plot: Residuals
    if(np.all(xerr == 0)):
        residuum_std = yerr
    else:  #calulate combined reisdual errors (sqrt((delY/delX)**2 *deltaX**2 + deltaY**2)
        residuum_std = np.sqrt((model.modelFunctionDerivative(xdata, *popt)*xerr)**2 + yerr**2)



    if plot_uncertainties:  # WHEN PLOTTING UNCERTAINTIES
        ax2.errorbar(xdata, residuum, yerr=residuum_std, fmt='.', elinewidth=1.5,  markersize = 10, capsize=0, label = 'Residuals', color = '#9467bd')
    else: 
            ax2.scatter(xdata, residuum, s = 8, label = 'Residuals', color = '#9467bd')

    ax2.axhline(y = 0, linestyle = '--', linewidth = 1.5, color = 'black')
    ax2.set_ylabel(plot2_ylabel, fontsize = 10)
    ax2.set_xlabel(xlabels, fontsize = 10)
    ax2.tick_params(axis='both', labelsize=10)
    ax2.grid(True)

    plt.tight_layout()

    #print results to latex
    with open(file_name + "-LatexFigure.txt", "w") as latexOutputFile:
        with redirect_stdout(latexOutputFile):
            utils_for_fpfit.print_fit_data_to_latex(model, algorithm, p0, file_name, chi2_dof, popt, popt_std,  pcorr, peak_index)


    if custom_plot == False:
        plt.savefig(file_name + '.pdf', format='pdf', bbox_inches='tight')
        plt.show()
        return popt, popt_std, pcorr
    else:   #return the ax1 object in addition to allow edits
        return popt, popt_std, pcorr, ax1