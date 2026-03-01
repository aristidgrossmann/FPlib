from fplib.models.ModelTemplate import ModelTemplate

class Linear(ModelTemplate):

    MODEL_NAME = 'linear'
    MODEL_EQUATION_LATEX = "$f(x) = ax $"
    MODEL_PARAMETER_LABELS = ["$a$"]

    @staticmethod
    def modelFunction(data, a):
        return a*data

    @classmethod
    def modelFunctionDerivative(cls, data, a):
        return a
