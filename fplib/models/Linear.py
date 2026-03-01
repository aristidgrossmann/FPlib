from fplib.ModelTemplate import ModelTemplate
import numpy as np

class Linear(ModelTemplate):

    MODEL_NAME = 'linear'
    MODEL_EQUATION_LATEX = "f(x) = ax"
    MODEL_PARAMETER_LABELS = ["$a$"]

    @staticmethod
    def modelFunction(data, a):
        return a*data

    @staticmethod
    def modelFunctionDerivative(data, a):
        return a*np.ones_like(data)
