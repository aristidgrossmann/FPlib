from fplib.models.ModelTemplate import ModelTemplate
from fplib.models.Voigt import Voigt

import numpy as np

class NVoigt(ModelTemplate):

    MODEL_NAME = 'N Voigt distributions'
    MODEL_EQUATION_LATEX = "f(x) = \\sum_{i = 1}^{N} Voigt(C_i, \\mu_i, \\sigma_i, \\gamma_i)"
    MODEL_PARAMETER_LABELS = []
    
    @staticmethod
    def modelFunction(data, *params):
        if (len(params))%4 != 0:
            raise ValueError("Wrong number of parameters")
        sum = 0
        for i in range(int(len(params)/4)):
            sum += Voigt.modelFunction(data, *params[4*i:4*i+4])
        return sum
    
    @classmethod
    def modelFunctionDerivative(cls, data, *params):

        #write parameter labels depending on length(params)
        cls.MODEL_PARAMETER_LABELS = []
        for i in range(int(len(params)/4)):
            cls.MODEL_PARAMETER_LABELS.append('$C_' + str(i+1) + '$')
            cls.MODEL_PARAMETER_LABELS.append('$\\mu_' + str(i+1) + '$')
            cls.MODEL_PARAMETER_LABELS.append('$\\sigma_' + str(i+1) + '$')
            cls.MODEL_PARAMETER_LABELS.append('$\\gamma_' + str(i+1) + '$')

        sum = 0
        for i in range(int(len(params)/4)):
            sum += Voigt.modelFunctionDerivative(data, *params[4*i:4*i+4])
        return sum
    

            
