from fplib.models.ModelTemplate import ModelTemplate
from fplib.models.Gaussian import Gaussian

import numpy as np

class NGaussians(ModelTemplate):

    MODEL_NAME = 'N Gaussians'
    MODEL_EQUATION_LATEX = "f(x) = \\sum_{i = 1}^{N} C_i*\\exp(-\\frac{(x - \\mu_i)^2}{2 \\sigma_i^2})"
    MODEL_PARAMETER_LABELS = []  #same order as Gaussian (is filled at runtime)
    
    @staticmethod
    def modelFunction(data, *params):
        if (len(params))%3 != 0:
            raise ValueError("Wrong number of parameters")
        sum = 0
        for i in range(int(len(params)/3)):
            sum += Gaussian.modelFunction(data, *params[3*i:3*i+3])
        return sum
    
    @classmethod
    def modelFunctionDerivative(cls, data, *params):

        #write parameter labels depending on length(params)
        cls.MODEL_PARAMETER_LABELS = []
        for i in range(int(len(params)/3)):
            cls.MODEL_PARAMETER_LABELS.append('$C_' + str(i+1) + '$')
            cls.MODEL_PARAMETER_LABELS.append('$\\mu_' + str(i+1) + '$')
            cls.MODEL_PARAMETER_LABELS.append('$\\sigma_' + str(i+1) + '$')

        #compute derivative
        sum = 0
        for i in range(int(len(params)/3)):
            sum += Gaussian.modelFunctionDerivative(data, *params[3*i:3*i+3])
        return sum
    

            
