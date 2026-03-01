from fplib.models.ModelTemplate import ModelTemplate
import numpy as np

class Gaussian(ModelTemplate):

    MODEL_NAME = 'Gaussian'
    MODEL_EQUATION_LATEX = "f(x) = C*\\exp(-\\frac{(x - \\mu)^2}{2 \\sigma^2})"
    MODEL_PARAMETER_LABELS = ["$C$", "$\\mu$", "\\sigma"]
    
    @staticmethod
    def modelFunction(data, C, mu, sigma):
        gauss = C*np.exp(-(data-mu)**2 / (2*sigma**2))
        return gauss
    
    @classmethod
    def modelFunctionDerivative(cls, data, C, mu, sigma):
        return C * np.exp(-(data - mu)**2 / (2 * sigma**2)) * (-(data - mu) / sigma**2)
        
