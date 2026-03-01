from fplib.models.Lorentz import Lorentz

class LorentzPlusConstant(Lorentz):

    MODEL_NAME = Lorentz.getName() + " + Constant"
    MODEL_EQUATION_LATEX = Lorentz.getName() + " + A"
    MODEL_PARAMETER_LABELS = Lorentz.getModelParameterLabels().append("$A$")
    
    @staticmethod
    def modelFunction(data, C, x0, gamma, A):
        return Lorentz.modelFunction(data, C, x0, gamma) + A
    
    @staticmethod
    def modelFunctionDerivative(data, C, x0, gamma, A):
        return Lorentz.modelFunctionDerivative(data, C, x0, gamma)
        
