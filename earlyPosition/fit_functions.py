"""
Define some functions you might want to use for fitting
"""
import numpy as np


def func0(t, c, m):
    """
    Linear function
    y= c + m*t
    """
    return( c + m * t)


def func1(t, c, a, tau, m):
    """
    Linear function with exponential
    y = c + m*t - a*e^(-t/tau)
    """
    return( (c + a) - (a * (np.exp(-t*tau))) + (m*t))


def func2v1(t, c, a, tauA, b, tauB):
    """
    Double exponential function
    y = c - a*e^(-t/tauA) + b*e^(t/tauB)
    """
    c=0
    return( -(a * (np.exp(-t/tauA))) + (b * (np.exp(-t/tauB))))


def func2(t, a, tauA, b, tauB):
    """
    Double exponential function
    y = c - a*e^(-t/tauA) + b*e^(t/tauB)
    """
    c=0
    return( -(a * (np.exp(-t/tauA))) + (b * (np.exp(-t/tauB))))


