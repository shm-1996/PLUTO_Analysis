import numpy as np
import readdbl_timerange as read
import ReferenceUnits as ref
import matplotlib.pyplot as plt
import PhysicalConstantsCGS as constant 
import matplotlib.patches as patches
import os
import sys
import argparse
from scipy.stats import gaussian_kde
from scipy.optimize import curve_fit
from scipy.optimize import OptimizeWarning
import warnings
import pickle
import tqdm

# Pickle Data Handling
def saveObj(obj, name):
    """
    Save a pickle object.

    INPUTS:
    ----------
    obj      - the name of the data object to be saved
    name     - the name of the pickle object

    """

    os.system("touch " + name + ".pkl")
    with open(name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)
        #pickle.dump(obj, f, protocol=2)


def loadObj(name):
    """
    Load a pickle object.

    INPUTS:
    ----------
    name     - the name of the pickle object

    """

    with open(name + '.pkl', 'rb') as f:
        return pickle.load(f)

