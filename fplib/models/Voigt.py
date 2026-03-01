from fplib.models.ModelTemplate import ModelTemplate
import numpy as np
import scipy.special

class Voigt(ModelTemplate):

    MODEL_NAME = 'Voigt'
    MODEL_EQUATION_LATEX = "Convolution of Gaussian and Lorentzian"
    MODEL_PARAMETER_LABELS = ["$C$", "$\\mu$", "$\\sigma$", "$\\gamma$"]
    
    @staticmethod
    def modelFunction(data, C, mu, sigma, gamma):
        return C*scipy.special.voigt_profile(data - mu, sigma, gamma)
    
    @classmethod
    def modelFunctionDerivative(cls, data, C, mu, sigma, gamma):
        z = ((data - mu) + 1j * gamma) / (sigma * np.sqrt(2))
        w = scipy.special.wofz(z)  # Faddeeva function
        dw_dx = (1j / (sigma * np.sqrt(2 * np.pi))) - (z * w) / (sigma**2 * np.sqrt(np.pi))
        return C*np.real(dw_dx)

        
