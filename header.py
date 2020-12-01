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

