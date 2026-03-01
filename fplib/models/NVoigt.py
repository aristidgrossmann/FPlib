from fplib.models.Voigt import Voigt

class NVoigt(Voigt):

    MODEL_NAME = 'N Voigt distributions'
    MODEL_EQUATION_LATEX = "f(x) = \\sum_{i = 1}^{N} C_i \\mathrm{Voigt}(x, \\mu_i, \\sigma_i, \\gamma_i)"
    MODEL_PARAMETER_LABELS = []  #same order as Voigt (is filled at runtime)
    
    @staticmethod
    def modelFunction(data, *params):
        if (len(params))%4 != 0:
            raise ValueError("Wrong number of parameters")
        sum = 0
        for i in range(int(len(params)/4)):
            sum += Voigt.modelFunction(data, *params[4*i:4*i+4])
        return sum
    
    @staticmethod
    def modelFunctionDerivative(data, *params):
        sum = 0
        for i in range(int(len(params)/4)):
            sum += Voigt.modelFunctionDerivative(data, *params[4*i:4*i+4])
        return sum
    
    @classmethod
    def getModelParameterLabels(cls, params):
        #write parameter labels depending on length(params)
        cls.MODEL_PARAMETER_LABELS = []
        for i in range(int(len(params)/4)):
            cls.MODEL_PARAMETER_LABELS.append('$C_' + str(i+1) + '$')
            cls.MODEL_PARAMETER_LABELS.append('$\\mu_' + str(i+1) + '$')
            cls.MODEL_PARAMETER_LABELS.append('$\\sigma_' + str(i+1) + '$')
            cls.MODEL_PARAMETER_LABELS.append('$\\gamma_' + str(i+1) + '$')

        return cls.MODEL_PARAMETER_LABELS
    
    
    

            
