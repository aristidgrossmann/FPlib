from fplib.ModelTemplate import ModelTemplate
import numpy as np

class DampedOscillation(ModelTemplate):

    MODEL_NAME = 'Damped Oscillation'
    MODEL_EQUATION_LATEX = "f(x) = C \\exp(- \\frac{t}{\\tau}) \\cos(\\omega t + \\phi)"
    MODEL_PARAMETER_LABELS = ["$C$", "$\\tau", "$\\omega","$\\phi"]

    @staticmethod
    def modelFunction(data, C, tau, omega, phi):
        return C*np.exp(-data/tau)*np.cos(omega*data + phi)
    
    @staticmethod
    def modelFunctionDerivative(data, C, tau, omega, phi):
        return - 1/tau *DampedOscillation.modelFunction(data, C, tau, omega, phi) - omega* C*np.exp(-data/tau)*np.sin(omega*data + phi)