####  FP LIBRARY ###
#this file contains functions for the FP

#### import native python libraries
import numpy as np
import scipy
from scipy.odr import ODR, Model, RealData
from scipy.signal import butter, filtfilt
from matplotlib import pyplot as plt
import scipy.special




###############################################################################################################################   FIT FUNCTIONS




def inverse_exponential_fit(data, const, char_length):
    return -1/char_length*np.log(data/const)


def linear_m1_bias_fit(data, bias):
    return data + bias*np.ones_like(data)





def cos_fit(data, A, omega, phi):
    return A*np.cos(omega*data + phi)

def lorentz_distribution_plus_C_fit(data, C, A, x0, Gamma):
    return C+ A/(2*np.pi)*Gamma/((data - x0)**2 + Gamma**2/4)

def N_lorentz_distribution_plus_C_fit(data, *params):

    if (len(params)-1)%3 != 0:
        raise ValueError("Wrong number of parameters")
    sum = params[0]

    for i in range(int((len(params)-1)/3)):
        sum += lorentz_distribution_fit(data, *params[3*i+1:3*i+4])
    return sum

def N_lorentz_distribution_plus_c_times_1_over_sin_2(data, *params):

    if (len(params)-2)%3 != 0:
        raise ValueError("Wrong number of parameters")

    v_max = np.max(np.abs(data))
    return (1/(1 + params[0]*np.sqrt(1-(data/v_max)**2)))*N_lorentz_distribution_plus_C_fit(data, *params[1:])


def abs_cos_fixed_frequency_fit(data, A, phi):
    omega = 2*np.pi/len(data)
    return A*np.abs(np.cos(data*omega + phi))

def abs_cos_variable_frequency_fit(data, A, omega, phi):
    return A*np.abs(np.cos(data*omega + phi))






def reflektometrie_fit(data, C, mu, omega_1, phi_1, omega_2, phi_2):                 # used in Roetngen: reflektrometrie
    return C*np.exp(data*mu)*np.cos(data*omega_1+ phi_1)*np.cos(data*omega_2 + phi_2)

def N_voigt_profile_fit_plus_exp_background(data, *params):         # used in Roetngen: reflektrometrie
    return  exponential_fit(data, params[0], params[1]) +  N_voigt_profile_fit(data, *params[2:])

def N_voigt_profile_fit_plus_polynomial_0_deg_exp_background(data, *params):            # used in Roetngen: reflektrometrie
    polynomial_0_deg_background = params[0] 
    return polynomial_0_deg_background + exponential_fit(data, params[1], params[2]) +  N_voigt_profile_fit(data, *params[3:])

def N_voigt_profile_fit_plus_polynomial_3rd_deg_exp_background(data, *params):      # used in Roetngen: reflektrometrie
    polynomial_3rd_deg_background = params[0] + params[1]*data + params[2]*data**2 + params[3]*data**3
    return polynomial_3rd_deg_background + exponential_fit(data, params[4], params[5]) +  N_voigt_profile_fit(data, *params[6:])


def real_beam_box(data, mu, I_0,  p, D):    # used in Stern-Gerlach
    ydata = np.zeros_like(data)
    for i in range(len(data)):
        if (data[i] - mu) > -D and (data[i] - mu) < -p:
            ydata[i] = D + (data[i] - mu)
        elif (data[i] - mu) >= -p and (data[i] - mu) <= p:
            ydata[i] = D -0.5*p -0.5*((data[i] - mu))**2/p
        elif (data[i] - mu) > p and (data[i] - mu) < D:
            ydata[i] = D - (data[i] - mu)
    return I_0*ydata

def Malus_law(data, N_min, N_max, phi_0):          # used in EPS
    return (N_max - N_min)*np.cos((data - phi_0)*np.pi/180)**2 + N_min




###############################################################################################################################   FIT FUNCTION DERIVATIVES

def gain_fit_derivative(data, a):
    return np.ones_like(data)*a

def const_fit_derivative(data, C):
    return np.zeros_like(data)


def inverse_exponential_fit_derivative(data, const, char_length):
    """Derivative of inverse_exponential_fit with respect to data."""
    return -1 / (char_length * data)  # Derivative of -1/char_length * log(data/const)



def gaussian_ccdf_fit_derivative(data, const, mu, sigma):
    return -gaussian_fit(data, const, mu, sigma)

def linear_m1_bias_fit_derivative(data, bias):
    return np.ones_like(data)



def cos_fit_derivative(data, A, omega, phi):
    return -A*omega*np.sin(data*omega + phi)

def lorentz_distribution_plus_C_fit_derivative(data, C, A, x0, Gamma):
    return - A*Gamma*(data-x0)/(np.pi*((data - x0)**2 + Gamma**2/4)**2)


def N_lorentz_distribution_plus_C_fit_derivative(data, *params):
    return N_lorentz_distribution_fit_derivative(data, *params[1:])

def N_lorentz_distribution_plus_c_times_1_over_sin_2_derivative(data, *params):

    v_max = np.max(np.abs(data))

    derivative_term1 = (params[0]*data)/(v_max**2*np.sqrt(1-(data/v_max)**2)*(1 + params[0]*np.sqrt(1-(data/v_max)**2))**2)*N_lorentz_distribution_plus_C_fit(data, *params[1:])
    derivative_term2 = (1/(1 + params[0]*np.sqrt(1-(data/v_max)**2)))*N_lorentz_distribution_plus_C_fit_derivative(data, *params[1:])
    return derivative_term1 + derivative_term2

def abs_cos_fixed_frequency_fit_derivative(data, A, phi):
    omega = 2*np.pi/len(data)
    return -A*omega*np.abs(np.sin(data*omega + phi))

def abs_cos_variable_frequency_derivative(data, A, omega, phi):
    return -A*omega*np.abs(np.sin(data*omega + phi))

def voigt_profile_fit_derivative(data, C, mu, sigma, gamma):

    z = ((data - mu) + 1j * gamma) / (sigma * np.sqrt(2))
    w = scipy.special.wofz(z)  # Faddeeva function
    dw_dx = (1j / (sigma * np.sqrt(2 * np.pi))) - (z * w) / (sigma**2 * np.sqrt(np.pi))
    return C*np.real(dw_dx)

def N_voigt_profile_fit_derivative(data, *params):
    sum = 0
    for i in range(int(len(params)/4)):
        sum += voigt_profile_fit_derivative(data, *params[4*i:4*i+4])
    return sum

def v_theory_derivative(v_exp, C, mu_e_mu_g_ratio, isomer_shift):   #only used in TO5 (Mössbauer)
    return 1.0

def reflektometrie_fit_derivative(data, C, mu, omega_1, phi_1, omega_2, phi_2):                 # used in Roetngen: reflektrometrie
    return C*np.exp(data*mu)*(np.cos(data*omega_1+ phi_1)*np.cos(data*omega_2 + phi_2)*mu \
                              - omega_1*np.sin(data*omega_1+ phi_1)*np.cos(data*omega_2 + phi_2) \
                              -np.cos(data*omega_1+ phi_1)*omega_2*np.sin(data*omega_2 + phi_2))

def N_voigt_profile_fit_plus_exp_background_derivative(data, *params):
    return exponential_fit_derivative(data, params[0], params[1])  + N_voigt_profile_fit_derivative(data, *params[2:])

def N_voigt_profile_fit_plus_polynomial_0_deg_exp_background_derivative(data, *params):
    polynomial_0_deg_background_derivative =  0
    return polynomial_0_deg_background_derivative + exponential_fit_derivative(data, params[1], params[2])  + N_voigt_profile_fit_derivative(data, *params[3:])

def N_voigt_profile_fit_plus_polynomial_3rd_deg_exp_background_derivative(data, *params):
    polynomial_3rd_deg_background_derivative =  params[1] + 2*params[2]*data + 3*params[3]*data**2
    return polynomial_3rd_deg_background_derivative + exponential_fit_derivative(data, params[4], params[5])  + N_voigt_profile_fit_derivative(data, *params[6:])

def real_beam_box_derivative(data, mu, I_0,  p, D):
    ydata = np.zeros_like(data)
    for i in range(len(data)):
        if (data[i] - mu) > -D and (data[i] - mu) < -p:
            ydata[i] = 1.0
        elif (data[i] - mu) >= -p and (data[i] - mu) <= p:
            ydata[i] = -((data[i] - mu))/p
        elif (data[i] - mu) > p and (data[i] - mu) < D:
            ydata[i] = -1.0
    return I_0*ydata

def Malus_law_derivative(data, N_min, N_max, phi_0):  
    return -2*np.pi/180*(N_max - N_min)*np.sin((data - phi_0)*np.pi/180)*np.cos((data - phi_0)*np.pi/180)



def get_fit_function_derivative(func, xdata, params):
    """Return the derivative of the given function."""
    if func.__name__ == linear_fit.__name__:
        return linear_fit_derivative(xdata, *params)
    
    elif func.__name__ == gain_fit.__name__:
        return gain_fit_derivative(xdata, *params)
    
    elif func.__name__ == const_fit.__name__:
        return const_fit_derivative(xdata, *params)
    
    elif func.__name__ == exponential_fit.__name__:
        return exponential_fit_derivative(xdata, *params)
    
    elif func.__name__ == log_linear_fit.__name__:
        return log_linear_fit_derivative(xdata, *params)
    
    elif func.__name__ == inverse_exponential_fit.__name__:
        return inverse_exponential_fit_derivative(xdata, *params)
    
    elif func.__name__ == gaussian_fit.__name__:
        return gaussian_fit_derivative(xdata, *params)
    
    elif func.__name__ == double_gaussian_fit.__name__:
        return double_gaussian_fit_derivative(xdata, *params)
    
    elif func.__name__ == N_gaussian_fit.__name__:
        return N_gaussian_fit_derivative(xdata, *params)
    
    elif func.__name__ == gaussian_ccdf_fit.__name__:
        return gaussian_ccdf_fit_derivative(xdata, *params)
    
    elif func.__name__ == linear_m1_bias_fit.__name__:
        return linear_m1_bias_fit_derivative(xdata, *params)
    
    elif func.__name__ == lorentz_distribution_fit.__name__:
        return lorentz_distribution_fit_derivative(xdata, *params)
    
    elif func.__name__ == N_lorentz_distribution_fit.__name__:
        return N_lorentz_distribution_fit_derivative(xdata, *params)
    
    elif func.__name__ == cos_fit.__name__:
        return cos_fit_derivative(xdata, *params)
    
    elif func.__name__ == abs_cos_fixed_frequency_fit.__name__:
        return abs_cos_fixed_frequency_fit_derivative(xdata, *params)
    
    elif func.__name__ == abs_cos_variable_frequency_fit.__name__:
        return abs_cos_variable_frequency_derivative(xdata, *params)
    
    elif func.__name__ == lorentz_distribution_plus_C_fit.__name__:
        return lorentz_distribution_plus_C_fit_derivative(xdata, *params)
    
    elif func.__name__ == N_lorentz_distribution_plus_C_fit.__name__:
        return N_lorentz_distribution_plus_C_fit_derivative(xdata, *params)
    
    elif func.__name__ == N_lorentz_distribution_plus_c_times_1_over_sin_2.__name__:
        return N_lorentz_distribution_plus_c_times_1_over_sin_2_derivative(xdata, *params)
    
    elif func.__name__ == voigt_profile_fit.__name__:
        return voigt_profile_fit_derivative(xdata, *params)
    
    elif func.__name__ == N_voigt_profile_fit.__name__:
        return N_voigt_profile_fit_derivative(xdata, *params)
    
    elif func.__name__ == v_theory.__name__:
        return v_theory_derivative(xdata, *params)
    
    elif func.__name__ == reflektometrie_fit.__name__:
        return reflektometrie_fit_derivative(xdata, *params)
    
    elif func.__name__ == N_voigt_profile_fit_plus_exp_background.__name__:
        return N_voigt_profile_fit_plus_exp_background_derivative(xdata, *params)
    
    elif func.__name__ == N_voigt_profile_fit_plus_polynomial_0_deg_exp_background.__name__:
        return N_voigt_profile_fit_plus_polynomial_0_deg_exp_background_derivative(xdata, *params)
    
    elif func.__name__ == N_voigt_profile_fit_plus_polynomial_3rd_deg_exp_background.__name__:
        return N_voigt_profile_fit_plus_polynomial_3rd_deg_exp_background_derivative(xdata, *params)
    
    elif func.__name__ == real_beam_box.__name__:
        return real_beam_box_derivative(xdata, *params)
    
    elif func.__name__ == Malus_law.__name__:
        return Malus_law_derivative(xdata, *params)
    else:
        raise ValueError("Function not supported for derivative calculation.")
      


def chi_squared_dof(residuals, xerr, yerr, dof):
    return np.sum(residuals**2 /(yerr**2 + xerr**2))/dof

def s(x1, err1, x2, err2):
    return abs(x1-x2)/np.sqrt(err1**2 + err2**2)

def combine_weighted(data, err):
    weights = 1/err**2

    mean = np.sum(data*weights)/np.sum(weights)
    mean_std = 1/np.sqrt(np.sum(weights))
    return [mean, mean_std]

###############################################################################################################################   MISCELLANEOUS FUNCTIONS
def zero_phase_filter(signal, cutoff_freq, sample_rate, filter_type='low', order=4):
    """
    Applies zero-phase digital filtering to a signal using forward-backward filtering
    
    Parameters:
    -----------
    signal : np.ndarray
        Input signal to be filtered
    cutoff_freq : float or [float, float]
        Cutoff frequency(ies) in Hz (for bandpass/bandstop as [low, high])
    sample_rate : float
        Sampling rate of the signal in Hz
    filter_type : str
        Filter type: 'low', 'high', 'bandpass', or 'bandstop'
    order : int
        Butterworth filter order
    
    Returns:
    --------
    filtered_signal : np.ndarray
        Zero-phase filtered signal
    """
    
    # Calculate Nyquist frequency
    nyquist = 0.5 * sample_rate
    
    # Normalize cutoff frequency(ies)
    normal_cutoff = np.array(cutoff_freq) / nyquist
    
    # Design Butterworth filter
    b, a = butter(order, normal_cutoff, btype=filter_type, analog=False)
    
    # Apply zero-phase filtering using forward-backward filtfilt
    return filtfilt(b, a, signal)


def curve_fit_to_odr_wrapper(func):
    """Wrap a curve_fit-style function to work with ODR."""
    def odr_style(B, x):
        return func(x, *B)
    return odr_style




##################################################################################################################    PRINTING FUNCTIONS: CURVE FIT RESULTS LATEX

def print_curve_fit_settings(fit_function_name, algorithm, p0):

    fit_function_name = fit_function_name.replace("_", " ")

    print("% Label and Table 1: Settings")
    print("\\begin{minipage}{0.2\\textwidth}")
    print("\\textbf{Settings}")
    print("\\end{minipage}%")
    print("\\begin{minipage}{0.75\\textwidth}")
    print("\\begin{tabular}{l|l}")
    print("\\hline")
    print("\\hline")
    print(f"Fit function & {fit_function_name} \\\\")
    print("\\hline")
    print(f"Algorithm & {algorithm} \\\\")
    print("\\hline")

    str_row = 'Starting Guess & ('
    for i in range(len(p0)):
        str_row += f"{p0[i]:.4g}" 
        if i < len(p0)-1:
            str_row += ', '
    str_row += ')' + '\\\\'
    print(str_row)
    print("\\hline")
    print("\\hline")
    print("\\end{tabular}")
    print("\\end{minipage}")
    print("")

    return None


def print_curve_fit_popt_general(chi2_dof, popt, popt_std, popt_labels):


    print("% Label and Table 2: Results")
    print("\\begin{minipage}{0.2\\textwidth}")
    print("\\textbf{Results}")
    print("\\end{minipage}%")


    if len(popt) > 5:

    
        noof_rows_1 = int(len(popt)/2)
        noof_rows_2 = len(popt) - noof_rows_1

        print("\\begin{minipage}{0.3\\textwidth}")
        print("\\begin{tabular}{c|c}")
        print("\\hline")
        print("\\hline")
        print(f"$\\chi^2/N_{{\\text{{dof}}}}$ & {chi2_dof:.4g} \\\\")
        print("\\hline")

        for i in range(noof_rows_1):
            str_row = popt_labels[i]
            str_row += ' & ' + f"{popt[i]:.4g}" + ' $\\pm$ ' + f"{popt_std[i]:.4g}"
            str_row += '\\\\'
            print(str_row)

        print("\\hline")
        print("\\hline")
        print("\\end{tabular}")
        print("\\end{minipage}")

        print("\\hspace{2cm}")

        print("\\begin{minipage}{0.3\\textwidth}")
        print("\\begin{tabular}{c|c}")
        print("\\hline")
        print("\\hline")

        for i in range(noof_rows_2):
            str_row = popt_labels[i+noof_rows_1]
            str_row += ' & ' + f"{popt[i+ noof_rows_1]:.4g}" + ' $\\pm$ ' + f"{popt_std[i+noof_rows_1]:.4g}"
            str_row += '\\\\'
            print(str_row)

        print("\\hline")
        print("\\hline")
        print("\\end{tabular}")
        print("\\end{minipage}")
        print("")


    else: 
        print("\\begin{minipage}{0.75\\textwidth}")
        print("\\begin{tabular}{c|c}")
        print("\\hline")
        print("\\hline")
        print(f"$\\chi^2/N_{{\\text{{dof}}}}$ & {chi2_dof:.4g} \\\\")
        print("\\hline")

        for i in range(len(popt)):
            str_row = popt_labels[i]
            str_row += ' & ' + f"{popt[i]:.4g}" + ' $\\pm$ ' + f"{popt_std[i]:.4g}"
            str_row += '\\\\'
            print(str_row)

        print("\\hline")
        print("\\hline")
        print("\\end{tabular}")
        print("\\end{minipage}")
        print("")

    return None



def print_correlation_table(pcorr, popt_labels):


    def print_corr_table_column_labels(popt_labels):
        str_column = ' '
        for i in range(1,len(popt_labels)):
            str_column += ' & ' + popt_labels[i]
        str_column += ' \\\\'
        print(str_column)
        return None

    def print_noof_columns_table(popt_labels):

        str_column_label = '\\begin{tabular}{l|'
        for i in range(len(popt_labels)-1):
            if i == len(popt_labels)-1:
                str_column_label+= 'c'
            else: 
                str_column_label += 'c'

        str_column_label += '}'
        print(str_column_label)
        return None

    print("% Label and Table 3: Correlation")
    print("\\begin{minipage}{0.2\\textwidth}")
    print("\\textbf{Correlation}")
    print("\\end{minipage}%")
    print("\\begin{minipage}{0.75\\textwidth}")

    print_noof_columns_table(popt_labels)

    print("\\hline")
    print("\\hline")

    print_corr_table_column_labels(popt_labels)
    print("\\hline")


    for i in range(len(popt_labels)-1):
        str_row = popt_labels[i]  # Start the row with the label
        for j in range(1,len(popt_labels)):
            if j > i:
                str_row += ' & ' + f"{pcorr[i][j]:.4g}"  # Append the correlation value
            else:
                str_row += ' & '  # Append an empty cell for lower triangular and diagonal
        str_row += '\\\\'  # End the row
        print(str_row)


    print("\\hline")
    print("\\hline")
    print("\\end{tabular}")
    print("\\end{minipage}")
    print("")

    return None

def print_curve_fit_settings_popt_pcorr_general(fit_function_name, algorithm, p0, chi2_dof, popt, popt_std, pcorr):

    if fit_function_name == linear_fit.__name__:   
        popt_labels = ["$a$", "$b$"]

    elif fit_function_name == gain_fit.__name__:   
        popt_labels = ["$a$"]

    elif fit_function_name == const_fit.__name__:
        popt_labels = ["$b$"]

    elif fit_function_name == gaussian_fit.__name__ or fit_function_name == gaussian_ccdf_fit.__name__:  #single gaussian fit
        popt_labels = ["$C_1$", "$\\mu_1$", "$\\sigma_1$"]

    elif fit_function_name == exponential_fit.__name__ or fit_function_name == log_linear_fit.__name__ or  \
            fit_function_name == inverse_exponential_fit.__name__: 
        popt_labels = ["$C$", "$\\mu$"]

    elif fit_function_name == linear_m1_bias_fit.__name__:
        popt_labels = ["b"]

    elif fit_function_name == lorentz_distribution_fit.__name__:
        popt_labels = ["$A$", "$x_0$", "$\\Gamma$"]

    elif fit_function_name == N_lorentz_distribution_fit.__name__:
        popt_labels = []
        for i in range(int(len(popt)/3)):
            popt_labels.append('$C_' + str(i+1) + '$')
            popt_labels.append('$\\mu_' + str(i+1) + '$')
            popt_labels.append('$\\Gamma_' + str(i+1) + '$')

    elif fit_function_name == N_gaussian_fit.__name__:
        popt_labels = []
        for i in range(int(len(popt)/3)):
            popt_labels.append('$C_' + str(i+1) + '$')
            popt_labels.append('$\\mu_' + str(i+1) + '$')
            popt_labels.append('$\\sigma_' + str(i+1) + '$')

    elif fit_function_name == lorentz_distribution_plus_C_fit.__name__:
        popt_labels = ["$C$", "$A$", "$\\mu$", "$\\Gamma$"]

    elif fit_function_name == N_lorentz_distribution_plus_C_fit.__name__:
        popt_labels = ["$C$"]
        for i in range(int((len(popt)-1)/3)):
            popt_labels.append('$A_' + str(i+1) + '$')
            popt_labels.append('$\\mu' + str(i+1) + '$')
            popt_labels.append('$\\Gamma_' + str(i+1) + '$')

    elif fit_function_name == cos_fit.__name__:
        popt_labels = ["$A$", "$\\omega$", "$\\phi$"]

    elif fit_function_name == N_lorentz_distribution_plus_c_times_1_over_sin_2.__name__:
        popt_labels = ["$A$", "$B$"]
        for i in range(int((len(popt)-2)/3)):
            popt_labels.append('$C_' + str(i+1) + '$')
            popt_labels.append('$\\mu_' + str(i+1) + '$')
            popt_labels.append('$\\Gamma_' + str(i+1) + '$')

    elif fit_function_name == abs_cos_fixed_frequency_fit.__name__:
        popt_labels = ["$A$", "$\\phi$"]

    elif fit_function_name == abs_cos_variable_frequency_fit.__name__:
        popt_labels = ["$A$", "$\\omega$", "$\\phi$"]

    elif fit_function_name == voigt_profile_fit.__name__:
        popt_labels = ["$C$", "$\\mu$", "$\\sigma$", "$\\gamma$"]

    elif fit_function_name == N_voigt_profile_fit.__name__:
        popt_labels = []
        for i in range(int(len(popt)/4)):
            popt_labels.append('$C_' + str(i+1) + '$')
            popt_labels.append('$\\mu_' + str(i+1) + '$')
            popt_labels.append('$\\sigma_' + str(i+1) + '$')
            popt_labels.append('$\\gamma_' + str(i+1) + '$')

    elif fit_function_name == v_theory.__name__:
        popt_labels = ["$Hc\\mu_a/E_0$", "$\\mu_e/\\mu_a$", "$\\delta E_{\\mathrm{iron}}$"]

    elif fit_function_name == reflektometrie_fit.__name__:
        popt_labels = ["$C$", "$\\mu$", "$\\omega_1$", "$\\phi_1$", "$\\omega_2$", "$\\phi_2$"]

    elif fit_function_name == N_voigt_profile_fit_plus_exp_background.__name__:
        popt_labels = ["$B$", "$\\mu$"]
        for i in range(int((len(popt)-2)/4)):
            popt_labels.append('$C_' + str(i+1) + '$')
            popt_labels.append('$\\mu_' + str(i+1) + '$')
            popt_labels.append('$\\sigma_' + str(i+1) + '$')
            popt_labels.append('$\\gamma_' + str(i+1) + '$')

    elif fit_function_name == N_voigt_profile_fit_plus_polynomial_0_deg_exp_background.__name__:
        popt_labels = ["$a_0$", "$B$", "$\\mu$"]
        for i in range(int((len(popt)-3)/4)):
            popt_labels.append('$C_' + str(i+1) + '$')
            popt_labels.append('$\\mu_' + str(i+1) + '$')
            popt_labels.append('$\\sigma_' + str(i+1) + '$')
            popt_labels.append('$\\gamma_' + str(i+1) + '$')

    elif fit_function_name == N_voigt_profile_fit_plus_polynomial_3rd_deg_exp_background.__name__:
        popt_labels = ["$a_0$", "$a_1$", "$a_2$", "$a_3$", "$B$", "$\\mu$"]
        for i in range(int((len(popt)-6)/4)):
            popt_labels.append('$C_' + str(i+1) + '$')
            popt_labels.append('$\\mu_' + str(i+1) + '$')
            popt_labels.append('$\\sigma_' + str(i+1) + '$')
            popt_labels.append('$\\gamma_' + str(i+1) + '$')

    elif fit_function_name == real_beam_box.__name__:
        popt_labels = ["$\\mu$", "$I_0$", "$p$", "$D$"]

    elif fit_function_name == Malus_law.__name__:
        popt_labels = ["$N_{\\mathrm{min}}$", "$N_{\\mathrm{max}}$", "$\\phi_0$"]


    print_curve_fit_settings(fit_function_name, algorithm, p0)    #print curve fit settings

    print("\\vspace{0.1cm} % Add vertical space between tables")
    print("")
    
    print_curve_fit_popt_general(chi2_dof, popt, popt_std, popt_labels)   #print curve fit results

    print("\\vspace{0.1cm} % Add vertical space between tables")
    print("")

    print_correlation_table(pcorr, popt_labels)   #print correlation table

    print("\\caption{Curvefit settings, optimized parameters and correlation between parameters}")
    print("\\label{fig:fit_results}")
    print("\\end{figure}") 

    return None



##############################################################################################################################################   PLOTTING FUNCTIONS




def general_curve_fit(xdata, ydata, xerr, yerr, model, p0, optimize_starting_guess,  # data and optimization settings
                     fit1_label, plot_uncertainties, xlabels, plot1_ylabel, plot2_ylabel, plot1_title,  #plot settings
                     plot2_title, plot1_legend_loc, plot2_legend_loc, file_name, peak_index = None, peak1_label = None, peak2_label = None,
                     exclude_0_err = False, selected_indices = None, custom_plot = False, xlims = None, ylims = None):
    
    ###################################################################################
    # perform chi2 curvefit (xerr = 0, yerr) or orthigonal distance regression (xerr, yerr)
    #get better starting values by fitting with (xerr = 0, yerr = 0)
    
    

    if np.any(selected_indices):

        all_xmin = np.min(xdata)
        all_xmax = np.max(xdata)
        all_xmodel = np.linspace(all_xmin, all_xmax, 1000)

        excluded_xdata = np.delete(xdata, selected_indices)
        excluded_xerr = np.delete(xerr, selected_indices)
        excluded_ydata = np.delete(ydata, selected_indices)
        excluded_yerr = np.delete(yerr, selected_indices)

        sorted_indices = np.argsort(excluded_xdata)

        excluded_xdata = excluded_xdata[sorted_indices]
        excluded_xerr = excluded_xerr[sorted_indices]
        excluded_ydata = excluded_ydata[sorted_indices]
        excluded_yerr = excluded_yerr[sorted_indices]


        xdata = xdata[selected_indices]
        xerr = xerr[selected_indices]
        ydata = ydata[selected_indices]
        yerr = yerr[selected_indices]
    



    if optimize_starting_guess:
        p0 = scipy.optimize.curve_fit(model, xdata= xdata, ydata = ydata,  p0 = p0)[0]
    if(np.all(xerr) == 0):
        algorithm = 'curve\\_fit'

        if exclude_0_err == True:   #exclude 0-count channels when True
            exclude_indices = []
            for i in range(len(xdata)):
                if yerr[i] == 0:
                    exclude_indices.append(i)

            xdata = np.delete(xdata, exclude_indices)
            ydata = np.delete(ydata, exclude_indices)
            xerr = np.delete(xerr, exclude_indices)
            yerr = np.delete(yerr, exclude_indices)

        popt, pcov = scipy.optimize.curve_fit(model, xdata= xdata, ydata = ydata,  p0 = p0, \
                                           sigma = yerr, absolute_sigma = True)



    else: 
        algorithm = 'ODR'

        if exclude_0_err == True:   #exclude 0-count channels when True
            exclude_indices = []
            for i in range(len(xdata)):
                if yerr[i] == 0 or xerr[i] == 0:
                    exclude_indices.append(i)

            xdata = np.delete(xdata, exclude_indices)
            ydata = np.delete(ydata, exclude_indices)
            xerr = np.delete(xerr, exclude_indices)
            yerr = np.delete(yerr, exclude_indices)

        data = RealData(xdata, ydata, sx=xerr, sy=yerr)

        # Create a Model object
        ODR_model = Model(curve_fit_to_odr_wrapper(model))

        # Create ODR object
        odr = ODR(data, ODR_model, beta0=p0)
        result = odr.run()
        popt, pcov = result.beta, result.cov_beta 

    #evaluate results
    popt_std = np.sqrt(np.diag(pcov))
    pcorr = pcov / np.outer(popt_std, popt_std)
    residuum = ydata - model(xdata, *popt)
    chi2_dof = chi_squared_dof(residuum, xerr*get_fit_function_derivative(model, xdata, popt), yerr, len(residuum)-len(popt))

    #print results
    print('-----   curve fit: ', algorithm, '-------')
    print('Fitted Function:  ', model.__name__)

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
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(6,7), gridspec_kw={'height_ratios': [2, 1]})

    if model.__name__ == v_theory.__name__:
        xmodel = np.append(np.arange(1,7), np.arange(1,7))
    else:
        xmodel = np.linspace(np.min(xdata), np.max(xdata), 1000)
    ymodel = model(xmodel, *popt)


    #########################################################    TOP PLOT
    #PLOT DATA
    if plot_uncertainties:  # WHEN PLOTTING UNCERTAINTIES
        if(np.all(xerr == 0)):
            ax1.errorbar(xdata, ydata, yerr=yerr, fmt='.', elinewidth=1.5,  markersize = 10, capsize=0, color = '#1f77b4', label = 'Data')
            if np.any(selected_indices):
                ax1.errorbar(excluded_xdata, excluded_ydata, yerr=excluded_yerr, fmt='.', elinewidth=1.5,  markersize = 10, capsize=0, color = '#1f77b4')
        else:
            ax1.errorbar(xdata, ydata, xerr=xerr, yerr=yerr, fmt='.', elinewidth=1.5,  markersize = 10, capsize=0, color = '#1f77b4', label = 'Data')
            if np.any(selected_indices):
                ax1.errorbar(excluded_xdata, excluded_ydata, xerr=excluded_xerr, yerr=excluded_yerr, fmt='.', elinewidth=1.5,  markersize = 10, capsize=0, color = '#1f77b4')


    else:
        ax1.scatter(xdata, ydata, s=8, color = '#1f77b4', label = 'Data')

        if np.any(selected_indices):
            ax1.scatter(excluded_xdata, excluded_ydata, s=8, color = '#1f77b4')


    # PLOT FIT


    if(model.__name__ == double_gaussian_fit.__name__):  # if double gaussian, plot individual distributions and their sum

        #determine which gaussian is background gaussian:
        if peak1_label == None and peak2_label == None:
            if peak_index == None:  #background index not specified
                peak_index = np.argmin(np.array([popt[2], popt[5]])) #determine peak index via smaller variance
                if peak_index == 0:
                    ax1.plot(xmodel, gaussian_fit(xmodel, *popt[3:6]), color = '#006400', linewidth = 2, label = 'Background Gaussian fit')
                    ax1.plot(xmodel, gaussian_fit(xmodel, *popt[0:3]), color = '#9467bd', linewidth = 2, label = 'Peak Gaussian fit')
                    
                else:
                    ax1.plot(xmodel, gaussian_fit(xmodel, *popt[0:3]), color = '#006400', linewidth = 2, label = 'Background Gaussian fit')
                    ax1.plot(xmodel, gaussian_fit(xmodel, *popt[3:6]), color = '#9467bd', linewidth = 2, label = 'Peak Gaussian fit')


            else:  #peak index is specified
                if peak_index == 0:
                    ax1.plot(xmodel, gaussian_fit(xmodel, *popt[3:6]), color = '#006400', linewidth = 2, label = 'Background Gaussian fit')
                    ax1.plot(xmodel, gaussian_fit(xmodel, *popt[0:3]), color = '#9467bd', linewidth = 2, label = 'Peak Gaussian fit')
                    
                else:
                    ax1.plot(xmodel, gaussian_fit(xmodel, *popt[0:3]), color = '#006400', linewidth = 2, label = 'Background Gaussian fit')
                    ax1.plot(xmodel, gaussian_fit(xmodel, *popt[3:6]), color = '#9467bd', linewidth = 2, label = 'Peak Gaussian fit')
        else: 
            ax1.plot(xmodel, gaussian_fit(xmodel, *popt[0:3]), color = '#006400', linewidth = 2, label = peak1_label)
            ax1.plot(xmodel, gaussian_fit(xmodel, *popt[3:6]), color = '#9467bd', linewidth = 2, label = peak2_label)

        ax1.plot(xmodel, ymodel, color = 'red', linewidth = 2, label = 'Sum Gaussian fit')

    elif model.__name__ == v_theory.__name__:
        ax1.scatter(xmodel, ymodel, color = 'red', s = 100, marker = 'x', linewidth = 2, label = fit1_label)

    else: #for all other models, just plot one line
        ax1.plot(xmodel, ymodel, color = 'red', linewidth = 2, label = fit1_label)

        #if np.any(selected_indices):
            #ax1.plot(all_xmodel, model(all_xmodel, *popt), linestyle = '--', color = 'red', linewidth = 2)



    ax1.set_ylabel(plot1_ylabel, fontsize = 12)
    ax1.set_xlabel(xlabels, fontsize = 12)
    ax1.tick_params(axis='both', labelsize=12)
    ax1.legend(fontsize = 12, loc = plot1_legend_loc)
    ax1.set_title(plot1_title, fontsize = 13)
    ax1.grid(True)

    # if loglinear fit (exponential fit, but logarithmic y axis)
    if(model.__name__ == log_linear_fit.__name__):
        ax1.set_yscale('log')


    if xlims is not None:
        ax1.set_xlim(xlims)
    if ylims is not None:
        ax1.set_ylim(ylims)



    ############################# Bottom plot: Residuals
    if(np.all(xerr == 0)):
        residuum_std = yerr
    else:  #calulate combined reisdual errors (sqrt((delY/delX)**2 *deltaX**2 + deltaY**2)
        residuum_std = np.sqrt((get_fit_function_derivative(model, xdata, popt)*xerr)**2 + yerr**2)



    if plot_uncertainties:  # WHEN PLOTTING UNCERTAINTIES
        ax2.errorbar(xdata, residuum, yerr=residuum_std, fmt='.', elinewidth=1.5,  markersize = 10, capsize=0, label = 'Residuals', color = '#9467bd')
    else: 
            ax2.scatter(xdata, residuum, s = 8, label = 'Residuals', color = '#9467bd')

    ax2.plot(xdata, np.zeros_like(xdata), '--', linewidth = 3, color = 'black')
    ax2.set_ylabel(plot2_ylabel, fontsize = 12)
    ax2.set_xlabel(xlabels, fontsize = 12)
    ax2.tick_params(axis='both', labelsize=12)
    ax2.legend(fontsize = 12, loc = plot2_legend_loc)
    ax2.grid(True)
    ax2.set_title(plot2_title + ' ' + r'($\chi^2/N_{\mathrm{dof}}$ = ' + f"{chi2_dof:.4g}" + ')', fontsize = 13)


    plt.tight_layout()

    file_name += '.pdf'


    #print results to latex
    if(model.__name__ == double_gaussian_fit.__name__): 
        print_fit_data_to_latex(model.__name__, algorithm, p0, chi2_dof, popt, popt_std,  pcorr, peak_index)
    else:
        print_fit_data_to_latex(model.__name__, algorithm, p0, chi2_dof, popt, popt_std,  pcorr)


    if custom_plot == False:
        plt.savefig(file_name, format='pdf', bbox_inches='tight')
        plt.show()
        return popt, popt_std, pcorr
    else:
        return popt, popt_std, pcorr, ax1



def print_fit_data_to_latex(fit_function_name, algorithm, p0,  #settings
                            chi2_dof, popt, popt_std,   #results
                            pcorr, peak_index = None):    #correlation matrix
    

    if fit_function_name ==  double_gaussian_fit.__name__:      #double gaussian fit

        popt_labels = ["$C_1$", "$\\mu_1$", "$\\sigma_1$", "$C_2$", "$\\mu_2$", "$\\sigma_2$"]

        print_curve_fit_settings(fit_function_name, algorithm, p0)    #print curve fit settings

        print("\\vspace{0.1cm} % Add vertical space between tables")
        print("")
        print("% Label and Table 2: Results")
        print("\\begin{minipage}{0.2\\textwidth}")
        print("\\textbf{Results}")
        print("\\end{minipage}%")
        print("\\begin{minipage}{0.75\\textwidth}")
        print("\\begin{tabular}{c|c|c}")
        print("\\hline")
        print("\\hline")
        print(f"& $\\chi^2/N_{{\\text{{dof}}}}$ & {chi2_dof:.4g} \\\\")
        print("\\hline")

        ### print optimized values
        if peak_index == 0:
            for i in range(len(popt)):
                str_row = ''
                if i == 0:
                    str_row += "\multirow{3}{*}{Peak}"

                if i == 3:
                    print('\\hline')
                    str_row += "\multirow{3}{*}{Background}"
                str_row += '& ' + popt_labels[i]
                str_row += ' & ' + f"{popt[i]:.4g}" + ' $\\pm$ ' + f"{popt_std[i]:.4g}"
                str_row += '\\\\'
                print(str_row)

        else: 
            for i in range(len(popt)):
                str_row = ''
                if i == 0:
                    str_row += "\multirow{3}{*}{Background}"

                if i == 3:
                    print('\\hline')
                    str_row += "\multirow{3}{*}{Peak}"
                str_row += '& ' + popt_labels[i]
                str_row += ' & ' + f"{popt[i]:.4g}" + ' $\\pm$ ' + f"{popt_std[i]:.4g}"
                str_row += '\\\\'
                print(str_row)

        print("\\hline")
        print("\\hline")
        print("\\end{tabular}")
        print("\\end{minipage}")
        print("")
        print("\\vspace{0.1cm} % Add vertical space between tables")
        print("")
        
        print_correlation_table(pcorr, popt_labels)   #print correlation table

        print("\\caption{Curvefit settings, optimized parameters and correlation between parameters}")
        print("\\label{fig:fit_results}")
        print("\\end{figure}") 


    else: 
        print_curve_fit_settings_popt_pcorr_general(fit_function_name, algorithm, p0, chi2_dof, popt, popt_std, pcorr)   #print results table (general form)
                                                                                                                                      #for special functions unique treatment







   

     