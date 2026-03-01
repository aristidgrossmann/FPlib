import numpy as np

from fplib.models.ModelTemplate import ModelTemplate



def chi_squared_dof(residuals, xerr, yerr, dof):
    return np.sum(residuals**2 /(yerr**2 + xerr**2))/dof

def curve_fit_to_odr_wrapper(func):
    """Wrap a curve_fit-style function to work with ODR."""
    def odr_style(B, x):
        return func(x, *B)
    return odr_style


#################################    PRINTING ANALYSIS TABLE TO LATEX   ##############################

def print_curve_fit_settings(modelName, algorithm, p0, file_name):

    print("\t" + "% Label and Table 1: Settings")
    print("\t" + "\\begin{minipage}{0.2\\textwidth}")
    print("\t\t" + "\\textbf{Settings}")
    print("\t" + "\\end{minipage}%")
    print("\t" + "\\begin{minipage}{0.75\\textwidth}")
    print("\t\t" + "\\begin{tabular}{l|l}")
    print("\t\t\t" + "\\hline")
    print("\t\t\t" + "\\hline")
    print("\t\t\t" + f"Fit function & {modelName} \\\\")
    print("\t\t\t" + "\\hline")
    print("\t\t\t" + f"Algorithm & {algorithm} \\\\")
    print("\t\t\t" + "\\hline")

    str_row = 'Starting Guess & ('
    for i in range(len(p0)):
        str_row += f"{p0[i]:.4g}" 
        if i < len(p0)-1:
            str_row += ', '
    str_row += ')' + '\\\\'
    print("\t\t\t" + str_row)
    print("\t\t\t" + "\\hline")
    print("\t\t\t" + "\\hline")
    print("\t\t" + "\\end{tabular}")
    print("\t" + "\\end{minipage}")
    print("")

    return None


def print_correlation_table(pcorr, popt_labels):


    def print_corr_table_column_labels(popt_labels):
        str_column = ''
        for i in range(1,len(popt_labels)):
            str_column += ' & ' + popt_labels[i]
        str_column += ' \\\\'
        print("\t\t\t" + str_column)
        return None

    def print_noof_columns_table(popt_labels):

        str_column_label = '\\begin{tabular}{l|'
        for i in range(len(popt_labels)-1):
            if i == len(popt_labels)-1:
                str_column_label+= 'c'
            else: 
                str_column_label += 'c'

        str_column_label += '}'
        print("\t\t" + str_column_label)
        return None

    print("\t" + "% Label and Table 3: Correlation")
    print("\t" + "\\begin{minipage}{0.2\\textwidth}")
    print("\t\t" + "\\textbf{Correlation}")
    print("\t" + "\\end{minipage}%")
    print("\t" + "\\begin{minipage}{0.75\\textwidth}")

    print_noof_columns_table(popt_labels)

    print("\t\t\t" + "\\hline")
    print("\t\t\t" + "\\hline")

    print_corr_table_column_labels(popt_labels)
    print("\t\t\t" + "\\hline")


    for i in range(len(popt_labels)-1):
        str_row = popt_labels[i]  # Start the row with the label
        for j in range(1,len(popt_labels)):
            if j > i:
                str_row += ' & ' + f"{pcorr[i][j]:.4g}"  # Append the correlation value
            else:
                str_row += ' & '  # Append an empty cell for lower triangular and diagonal
        str_row += '\\\\'  # End the row
        print("\t\t\t" + str_row)


    print("\t\t\t" + "\\hline")
    print("\t\t\t" + "\\hline")
    print("\t\t" + "\\end{tabular}")
    print("\t" + "\\end{minipage}")
    print("")

    return None


def print_fit_data_to_latex(model: type[ModelTemplate], algorithm, p0, file_name, chi2_dof, popt, popt_std, pcorr, peak_index = None):

    print("\\begin{figure}[H]")
    print("\t" + "\\centering")
    print("\t" + "\includegraphics{" + file_name + ".pdf}")
    print('')

    print("\t" + "\\begin{equation*}")
    print("\t\t" + model.getModelEquationLatex())
    print("\t" + "\\end{equation*}")
    print("")

    print_curve_fit_settings(model.getName(), algorithm, p0, file_name)    #print curve fit settings

    print("\t" + "\\vspace{0.1cm} % Add vertical space between tables")
    print("")
    
    model.print_curve_fit_popt_general(chi2_dof, popt, popt_std, peak_index)

    print("\t" + "\\vspace{0.1cm} % Add vertical space between tables")
    print("")

    print_correlation_table(pcorr, model.getModelParameterLabels())   #print correlation table

    print("\t" + "\\caption{Curvefit settings, optimized parameters and correlation between parameters}")
    print("\t" + "\\label{fig:fit_results}")
    print("\\end{figure}") 

    return None
