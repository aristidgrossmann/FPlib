from fplib.models.Voigt import Voigt

class VoigtPlusConstant(Voigt):

    MODEL_NAME = Voigt.getName() + " + Constant"
    MODEL_EQUATION_LATEX = Voigt.getName() + " + A"
    MODEL_PARAMETER_LABELS = Voigt.getModelParameterLabels().append("$A$")
    
    @staticmethod
    def modelFunction(data, C, mu, sigma, gamma, A):
        return Voigt.modelFunction(data, C, mu, sigma, gamma) + A
    
    @staticmethod
    def modelFunctionDerivative(data, C, mu, sigma, gamma, A):
        return super().modelFunctionDerivative(data, C, mu, sigma, gamma)
        
