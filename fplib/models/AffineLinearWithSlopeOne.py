from fplib.models.AffineLinear import AffineLinear

class AffineLinearWithSlopeOne(AffineLinear):

    MODEL_NAME = AffineLinear.getName() + " with Slope = 1"
    MODEL_EQUATION_LATEX = "f(x) = x + b"
    MODEL_PARAMETER_LABELS = ["$b$"]

    @staticmethod
    def modelFunction(data, b):
        return AffineLinear.modelFunction(data, 1, b)

    @staticmethod
    def modelFunctionDerivative(data, b):
        return AffineLinear.modelFunctionDerivative(data, 1, b)
