from fplib.models.ModelTemplate import ModelTemplate
import numpy as np

class Constant(ModelTemplate):

    MODEL_NAME = 'Constant'
    MODEL_EQUATION_LATEX = "$f(x) = C$"
    MODEL_PARAMETER_LABELS = ["$C$"]

    @staticmethod
    def modelFunction(data, C):
        return C*np.ones_like(data)

    @classmethod
    def modelFunctionDerivative(cls, data, C):
        return np.zeros_like(data)
