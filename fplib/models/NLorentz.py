from fplib.models.ModelTemplate import ModelTemplate
from fplib.models.Lorentz import Lorentz

import numpy as np

class NLorentz(ModelTemplate):

    MODEL_NAME = 'N Lorentz distributions'
    MODEL_EQUATION_LATEX = "f(x) = \\sum_{i = 1}^{N} \\frac{C_i \\gamma_i}{2\\pi (x - x_{0,i})^2 + \\gamma_i^2/4}"
    MODEL_PARAMETER_LABELS = []  #same order as Lorentz (is filled at runtime)
    
    @staticmethod
    def modelFunction(data, *params):
        if (len(params))%3 != 0:
            raise ValueError("Wrong number of parameters")
        sum = 0
        for i in range(int(len(params)/3)):
            sum += Lorentz.modelFunction(data, *params[3*i:3*i+3])
        return sum
    
    @classmethod
    def modelFunctionDerivative(cls, data, *params):

        #write parameter labels depending on length(params)
        cls.MODEL_PARAMETER_LABELS = []
        for i in range(int(len(params)/3)):
            cls.MODEL_PARAMETER_LABELS.append('$C_' + str(i+1) + '$')
            cls.MODEL_PARAMETER_LABELS.append('$x_{0, ' + str(i+1) + '}$')
            cls.MODEL_PARAMETER_LABELS.append('$\\gamma_' + str(i+1) + '$')

        sum = 0
        for i in range(int(len(params)/3)):
            sum += Lorentz.modelFunctionDerivative(data, *params[3*i:3*i+3])
        return sum
    

            
