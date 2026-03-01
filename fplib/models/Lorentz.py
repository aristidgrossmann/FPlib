from fplib.models.ModelTemplate import ModelTemplate
import numpy as np

class Lorentz(ModelTemplate):

    MODEL_NAME = 'Lorentz Distribution'
    MODEL_EQUATION_LATEX = "f(x) = \\frac{C \\gamma}{2\\pi (x - x_0)^2 + \\gamma^2/4}"
    MODEL_PARAMETER_LABELS = ["$C$", "$x_0", "\\gamma"]
    
    @staticmethod
    def modelFunction(data, C, x0, gamma):
        return C/(2*np.pi)*gamma/((data - x0)**2 + gamma**2/4)
    
    @classmethod
    def modelFunctionDerivative(cls, data, C, x0, gamma):
        return - C*gamma*(data-x0)/(np.pi*((data - x0)**2 + gamma**2/4)**2)
        
