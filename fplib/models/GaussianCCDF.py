from fplib.models.ModelTemplate import ModelTemplate
from fplib.models.Gaussian import Gaussian
import numpy as np
import scipy.special

class GaussianCCDF(ModelTemplate):

    MODEL_NAME = 'Gaussian CCDF'
    MODEL_EQUATION_LATEX = "f(x) = \\frac{C}{2} \\left(1 - \\Phi(\\mu, \\sigma) \\right)"
    MODEL_PARAMETER_LABELS = ["$C$", "$\\mu$", "\\sigma"]
    
    @staticmethod
    def modelFunction(data, const, mu, sigma):
        return const*0.5*(1 - scipy.special.erf((data - mu) / (sigma * np.sqrt(2))))
    
    @classmethod
    def modelFunctionDerivative(cls, data, const, mu, sigma):
        return -Gaussian.modelFunction(data, const, mu, sigma)
        
