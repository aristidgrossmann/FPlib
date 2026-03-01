from fplib.models.ModelTemplate import ModelTemplate
import numpy as np

class AbsCosine(ModelTemplate):

    MODEL_NAME = 'Cosine'
    MODEL_EQUATION_LATEX = "f(x) = A \\cos(\\omega x + \\phi)"
    MODEL_PARAMETER_LABELS = ["$A$", "$\\omega$", "$\\phi$"]

    @staticmethod
    def modelFunction(data, A, omega, phi):
        return A*np.cos(data*omega + phi)

    @staticmethod
    def modelFunctionDerivative(data, A, omega, phi):
        return -A*omega*np.sin(data*omega + phi)
