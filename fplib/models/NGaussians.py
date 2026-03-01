from fplib.models.Gaussian import Gaussian

class NGaussians(Gaussian):

    MODEL_NAME = 'N Gaussians'
    MODEL_EQUATION_LATEX = "f(x) = \\sum_{i = 1}^{N} C_i e^{-\\frac{(x - \\mu_i)^2}{2 \\sigma_i^2}}"
    MODEL_PARAMETER_LABELS = []  #same order as Gaussian (is constructed at runtime)
    
    @staticmethod
    def modelFunction(data, *params):
        if (len(params))%3 != 0:
            raise ValueError("Wrong number of parameters")
        sum = 0
        for i in range(int(len(params)/3)):
            sum += Gaussian.modelFunction(data, *params[3*i:3*i+3])
        return sum
    
    @staticmethod
    def modelFunctionDerivative(data, *params):
        sum = 0
        for i in range(int(len(params)/3)):
            sum += Gaussian.modelFunctionDerivative(data, *params[3*i:3*i+3])
        return sum
    
    @classmethod
    def getModelParameterLabels(cls, params):
        #write parameter labels depending on length(params)
        cls.MODEL_PARAMETER_LABELS = []
        for i in range(int(len(params)/3)):
            cls.MODEL_PARAMETER_LABELS.append('$C_' + str(i+1) + '$')
            cls.MODEL_PARAMETER_LABELS.append('$\\mu_' + str(i+1) + '$')
            cls.MODEL_PARAMETER_LABELS.append('$\\sigma_' + str(i+1) + '$')

        return cls.MODEL_PARAMETER_LABELS
    

            
