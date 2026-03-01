from fplib.ModelTemplate import ModelTemplate
import numpy as np

class Exponential(ModelTemplate):

    MODEL_NAME = 'exponential'
    MODEL_EQUATION_LATEX = "f(x) = C \\exp(-\\mu x)"
    MODEL_PARAMETER_LABELS = ["$C$", "$\\mu$"]

    @staticmethod
    def modelFunction(data, C, mu):
        return C*np.exp(-data*mu)
    
    @staticmethod
    def modelFunctionDerivative(data, C, mu):
        return -mu*C*np.exp(-data*mu)