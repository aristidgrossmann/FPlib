from fplib.models.ModelTemplate import ModelTemplate
from fplib.models.Exponential import Exponential
import numpy as np

class LogLinear(Exponential):  #same as exponential, but with logarithmic y-scaling 

    MODEL_NAME = 'Log-Linear'

    #override scaling to 'log': 
    X_SCALE = 'linear'
    Y_SCALE = 'log'