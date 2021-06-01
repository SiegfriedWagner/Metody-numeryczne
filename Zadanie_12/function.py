import numpy as np

period = 4


def fun(x):
    return np.abs(np.mod(x, period) - 2) - 1
