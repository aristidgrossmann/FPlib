from fplib.models.ModelTemplate import ModelTemplate

class AffineLinear(ModelTemplate):

    MODEL_NAME = 'Affine-Linear'
    MODEL_EQUATION_LATEX = "$f(x) = ax + b$"
    MODEL_PARAMETER_LABELS = ["$a$", "$b$"]

    @staticmethod
    def modelFunction(data, a, b):
        return a*data + b

    @classmethod
    def modelFunctionDerivative(cls, data, a, b):
        return a
