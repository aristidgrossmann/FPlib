from fplib.models.ModelTemplate import ModelTemplate
import numpy as np

class DoubleGaussian(ModelTemplate):

    MODEL_NAME = 'Double Gaussian'
    MODEL_EQUATION_LATEX = "f(x) = C_1*\\exp(-\\frac{(x - \\mu_1)^2}{2 \\sigma_1^2}) + C_2*\\exp(-\\frac{(x - \\mu_2)^2}{2 \\sigma_2^2})"
    MODEL_PARAMETER_LABELS = ["$C_1$", "$\\mu_1$", "$\\sigma_1$", "$C_2$", "$\\mu_2$", "$\\sigma_2$"]
    
    @staticmethod
    def modelFunction(data, const1, mu1, sigma1, const2, mu2, sigma2):
        gauss1 = const1*np.exp(-(data-mu1)**2 / (2*sigma1**2))
        gauss2 = const2*np.exp(-(data-mu2)**2 / (2*sigma2**2))
        return gauss1 + gauss2
    
    @classmethod
    def modelFunctionDerivative(cls, data, const1, mu1, sigma1, const2, mu2, sigma2):
        gauss1_derivative = const1 * np.exp(-(data - mu1)**2 / (2 * sigma1**2)) * (-(data - mu1) / sigma1**2)
        gauss2_derivative = const2 * np.exp(-(data - mu2)**2 / (2 * sigma2**2)) * (-(data - mu2) / sigma2**2)
        return gauss1_derivative + gauss2_derivative
    
    @classmethod
    def print_curve_fit_popt_general(cls, chi2_dof, popt, popt_std, peak_index) -> None:
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
                    str_row += "\\multirow{3}{*}{Peak}"

                if i == 3:
                    print('\\hline')
                    str_row += "\\multirow{3}{*}{Background}"
                str_row += '& ' + cls.MODEL_PARAMETER_LABELS[i]
                str_row += ' & ' + f"{popt[i]:.4g}" + ' $\\pm$ ' + f"{popt_std[i]:.4g}"
                str_row += '\\\\'
                print(str_row)

        else: 
            for i in range(len(popt)):
                str_row = ''
                if i == 0:
                    str_row += "\\multirow{3}{*}{Background}"

                if i == 3:
                    print('\\hline')
                    str_row += "\\multirow{3}{*}{Peak}"
                str_row += '& ' + cls.MODEL_PARAMETER_LABELS[i]
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

        return None
        
