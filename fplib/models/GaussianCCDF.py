from fplib.ModelTemplate import ModelTemplate
from fplib.models.Gaussian import Gaussian
import numpy as np
import scipy.special

class GaussianCCDF(ModelTemplate):

    MODEL_NAME = 'Gaussian CCDF'
    MODEL_EQUATION_LATEX = "f(x) = C \\left(1 - \\Phi(x, \\mu, \\sigma) \\right)"
    MODEL_PARAMETER_LABELS = ["$C$", "$\\mu$", "$\\sigma$"]
    
    @staticmethod
    def modelFunction(data, C, mu, sigma):
        return C*0.5*(1 - scipy.special.erf((data - mu) / (sigma * np.sqrt(2))))
    
    @staticmethod
    def modelFunctionDerivative(data, C, mu, sigma):
        return -Gaussian.modelFunction(data, C, mu, sigma)
        
