import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from astropy.cosmology import FlatLambdaCDM
from astropy.table import Table
from scipy import stats
import seaborn as sns
from scipy.stats import skewnorm
import math as m
from numpy import polyfit as pfit
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
plt.rcParams["figure.figsize"] = [10, 8]
plt.rcParams.update({'font.size': 18})

#%%
## Table 1 from Perlmutter et al 1999
file = Table.read('C:/Users/hchri/OneDrive/Documents/GitHub/SN_project/Table1_perlmutter99.csv')
df = file.to_pandas()

#%%
# get luminosity distance
dL = cosmo.luminosity_distance(df.z)
dL_err = cosmo.luminosity_distance(df.z_err)
