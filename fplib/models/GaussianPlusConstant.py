from fplib.models.Gaussian import Gaussian

class GaussianPlusConstant(Gaussian):

    MODEL_NAME = Gaussian.getName() + " + Constant"
    MODEL_EQUATION_LATEX = Gaussian.getName() + " + A"
    MODEL_PARAMETER_LABELS = Gaussian.getModelParameterLabels().append("$A$")
    
    @staticmethod
    def modelFunction(data, C, mu, sigma, A):
        return Gaussian.modelFunction(data, C, mu, sigma) + A
    
    @staticmethod
    def modelFunctionDerivative(data, C, mu, sigma, A):
        return super().modelFunctionDerivative(data, C, mu, sigma)
        
