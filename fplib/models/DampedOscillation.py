from fplib.ModelTemplate import ModelTemplate
import numpy as np

class DampedOscillation(ModelTemplate):

    MODEL_NAME = 'Damped Oscillation'
    MODEL_EQUATION_LATEX = "f(x) = A\\exp(- \\frac{t}{\\tau}) \\cos(\\omega t + \\phi)"
    MODEL_PARAMETER_LABELS = ["$A$", "$\\tau$", "$\\omega$","$\\phi$"]

    @staticmethod
    def modelFunction(data, A, tau, omega, phi):
        return A*np.exp(-data/tau)*np.cos(omega*data + phi)
    
    @staticmethod
    def modelFunctionDerivative(data, A, tau, omega, phi):
        return - 1/tau *DampedOscillation.modelFunction(data, A, tau, omega, phi) - omega* A*np.exp(-data/tau)*np.sin(omega*data + phi)