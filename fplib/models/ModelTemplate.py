from abc import ABC, abstractmethod
import numpy as np

class ModelTemplate(ABC):

    # must be defined when inheriting
    MODEL_NAME :str = 'modelTemplate'
    MODEL_EQUATION_LATEX: str = 'modelTemplate'
    MODEL_PARAMETER_LABELS: list[str] = []

    #override for 'log'
    X_SCALE = 'linear'
    Y_SCALE = 'linear'
    
    @classmethod
    def __init_subclass__(cls):
        super().__init_subclass__()
        for attr in ["MODEL_NAME", "MODEL_EQUATION_LATEX", "MODEL_PARAMETER_LABELS"]:
            if not hasattr(cls, attr):
                raise TypeError(f"{cls.__name__} must define {attr}")
            
    @staticmethod
    @abstractmethod  #must be defined when inheriting
    def modelFunction(data, *params):
        pass

    @classmethod
    @abstractmethod  #must be defined when inheriting
    def modelFunctionDerivative(cls, data, *params):
        pass


    ##########  DO NOT TOUCH  #########
    @classmethod
    def getName(cls) -> str:
        return cls.MODEL_NAME
    
    @classmethod
    def getModelEquationLatex(cls) -> str:
        return cls.MODEL_EQUATION_LATEX
    
    @classmethod
    def getModelParameterLabels(cls) -> list[str]:
        return cls.MODEL_PARAMETER_LABELS
    
    @classmethod
    def getXscale(cls):
        return cls.X_SCALE
    
    @classmethod
    def getYscale(cls):
        return cls.Y_SCALE

    ##########  CHANGE ONLY WHEN FORMATTING IS OFF  #########
    @classmethod
    def print_curve_fit_popt_general(cls, chi2_dof, popt, popt_std, peak_index = None) -> None:

        print("\t" + "% Label and Table 2: Results")
        print("\t" + "\\begin{minipage}{0.2\\textwidth}")
        print("\t\t" + "\\textbf{Results}")
        print("\t" + "\\end{minipage}%")

        if len(popt) > 5:
            noof_rows_1 = int(len(popt)/2)
            noof_rows_2 = len(popt) - noof_rows_1

            print("\t" + "\\begin{minipage}{0.3\\textwidth}")
            print("\t\t" + "\\begin{tabular}{c|c}")
            print("\t\t\t" + "\\hline")
            print("\t\t\t" + "\\hline")
            print("\t\t\t" + f"$\\chi^2/N_{{\\text{{dof}}}}$ & {chi2_dof:.4g} \\\\")
            print("\t\t\t" + "\\hline")

            for i in range(noof_rows_1):
                str_row = cls.MODEL_PARAMETER_LABELS[i]
                str_row += ' & ' + f"{popt[i]:.4g}" + ' $\\pm$ ' + f"{popt_std[i]:.4g}"
                str_row += '\\\\'
                print("\t\t\t" + str_row)

            print("\t\t\t" + "\\hline")
            print("\t\t\t" + "\\hline")
            print("\t\t" + "\\end{tabular}")
            print("\t" + "\\end{minipage}")

            print("\t" + "\\hspace{2cm}")

            print("\t" + "\\begin{minipage}{0.3\\textwidth}")
            print("\t\t" + "\\begin{tabular}{c|c}")
            print("\t\t\t" + "\\hline")
            print("\t\t\t" + "\\hline")

            for i in range(noof_rows_2):
                str_row = cls.MODEL_PARAMETER_LABELS[i+noof_rows_1]
                str_row += ' & ' + f"{popt[i+ noof_rows_1]:.4g}" + ' $\\pm$ ' + f"{popt_std[i+noof_rows_1]:.4g}"
                str_row += '\\\\'
                print("\t\t\t" + str_row)

            print("\t\t\t" + "\\hline")
            print("\t\t\t" + "\\hline")
            print("\t\t" + "\\end{tabular}")
            print("\t" + "\\end{minipage}")
            print("")

        else: 
            print("\t" + "\\begin{minipage}{0.75\\textwidth}")
            print("\t\t" + "\\begin{tabular}{c|c}")
            print("\t\t\t" + "\\hline")
            print("\t\t\t" + "\\hline")
            print("\t\t\t" + f"$\\chi^2/N_{{\\text{{dof}}}}$ & {chi2_dof:.4g} \\\\")
            print("\t\t\t" + "\\hline")

            for i in range(len(popt)):
                str_row = cls.MODEL_PARAMETER_LABELS[i]
                str_row += ' & ' + f"{popt[i]:.4g}" + ' $\\pm$ ' + f"{popt_std[i]:.4g}"
                str_row += '\\\\'
                print("\t\t\t" + str_row)

            print("\t\t\t" + "\\hline")
            print("\t\t\t" + "\\hline")
            print("\t\t" + "\\end{tabular}")
            print("\t" + "\\end{minipage}")
            print("")

        return None