from fplib.models.ModelTemplate import ModelTemplate

class NDegreePolynomial(ModelTemplate):

    MODEL_NAME = 'N Degree Polynomial'
    MODEL_EQUATION_LATEX = "f(x) = \\sum_{k = 0}^{N} a_{k} x^k"
    MODEL_PARAMETER_LABELS = []  # a_0, a_1, ... (is filled at runtime)
    
    @staticmethod
    def modelFunction(data, *params):
        if (len(params)) == 0:
            raise ValueError("Wrong number of parameters")
        sum = 0
        for k in range(int(len(params))):
            sum += params[k]*data**k
        return sum
    
    @staticmethod
    def modelFunctionDerivative(data, *params):
        sum = 0
        for k in range(int(len(params))):
            sum += k*params[k]*data**(k-1)
        return sum
    
    @classmethod
    def getModelParameterLabels(cls, params):
        #write parameter labels depending on length(params)
        cls.MODEL_PARAMETER_LABELS = []
        for i in range(int(len(params))):
            cls.MODEL_PARAMETER_LABELS.append('$a_' + str(i) + '$')
        return cls.MODEL_PARAMETER_LABELS
    

            
