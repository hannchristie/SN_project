#import necessary packages
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from astropy.cosmology import FlatLambdaCDM
from astropy.table import Table
from scipy import stats
import seaborn as sns
from scipy.stats import skewnorm
from scipy.optimize import least_squares
from scipy.optimize import curve_fit
import scipy.integrate as integrate 
from scipy.integrate import quad
import math 
from numpy import polyfit as pfit
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
plt.rcParams["figure.figsize"] = [10, 8]
plt.rcParams.update({'font.size': 18})

#%%
# read in data into two separate dataframes
## Table 1 from Perlmutter et al 1999
file = Table.read('C:/Users/hchri/OneDrive/Documents/GitHub/SN_project/Table1_perlmutter99.csv')
df = file.to_pandas()

## Table 2 from Permutter et al 1999
## Table 1 from Perlmutter et al 1999
file1 = Table.read('C:/Users/hchri/OneDrive/Documents/GitHub/SN_project/Table2_perlmutter99.csv')
df1 = file1.to_pandas()

#%%
# define constants
c = 3*10**5 #[km/s]
hubble = 70 #[km s^-1 Mpc^-1]

#%%
# find the parameters for the relation between apparent magnitude and redshift
# define the curve fit function 
def curve_func(x, a, b): # takes in the variable and any necessary parameters
    return a*np.log10(x) + b # returns the expected function

a, b = curve_fit(curve_func, df1.z, df1.mb_corr)

#%%
# checking values here
print(a,b)

#print(slope*df1.z+intercept)

#print(df1.z)

#%%
print(np.log10(c*df1.z))

#%%
# find the value of 
# isolate for the value of "fancy" M, log rules...
fancy_M = a[1] - 5*np.log10(c)
print(fancy_M)


#%%
# try with the linear fit
slope, intercept, r, p, sterr = stats.linregress(np.log10(c*df1.z), df1.mb_corr)
print(slope, intercept)
#%%
# plot the apparent magnitude as a function of redshfit
x_vals = np.linspace(min(df1.z), max(df1.z))
#plt.scatter(df1.z, df1.mb_corr)
plt.scatter(np.log10(c*df1.z), df1.mb_corr)
#plt.errorbar(np.log10(c*df1.z), df1.mb_corr,  df1.mb_err, df1.z_err, fmt = 'None')
plt.plot(x_vals, slope*np.log10(c*x_vals) + intercept, label = 'lsq fit', c = 'black', linestyle = '--')

plt.ylabel('m_B Corrected')
plt.xlabel('redshift')
#plt.xlim(0, 0.2)
#plt.xscale('log')
plt.ylim(14, 20)
plt.legend()

#%%
# plot the apparent magnitude as a function of redshfit
x_vals = np.linspace(min(df1.z), max(df1.z))
plt.scatter(np.log10(c*df1.z), df1.mb_corr)
#plt.errorbar(df1.z, df1.mb_corr,  df1.mb_err, df1.z_err, fmt = 'None')
#plt.plot(x_vals, a[0]*x_vals + fancy_M, label = '5.03*$log_{10}(z)$ - 4.693', c = 'black', linestyle = '--')
plt.plot(np.log10(c*(x_vals)), slope*np.log10(c*(x_vals)) + intercept, label = 'LSQ fit', c = 'black', linestyle = '--')
plt.ylabel('m_B Corrected')
plt.xlabel('$log_{10}(cz)$')
#plt.xlim(0, 0.2)
#plt.xscale('log')
plt.ylim(14, 20)
plt.legend()

#%%

#%%
# combine the two data sets in order to do fits with them together
z_both = [df1.z, df.z]
zerr_both = [df1.z_err, df.z_err]
redshift = pd.concat(z_both)
redshift_err = pd.concat(zerr_both)

print(redshift) # checking to make sure the concat function worked properly
#%%
mb_both = [df1.mb_corr, df.mb_eff]
app_mag = pd.concat(mb_both)
app_err = [df1.mb_err, df.mb_err]
appmag_err = pd.concat(app_err)

#%%
# plot the apparent magnitude as a funciton of redhsift
# plot the two different data sets with different colours
x_vals = np.linspace(min(df1.z), max(df.z))
plt.scatter(df1.z, df1.mb_corr, c = 'firebrick', label = 'Low Redshift SN')
plt.errorbar(df1.z, df1.mb_corr,  df1.mb_err, df1.z_err, fmt = 'None', c = 'black', alpha = 0.5)
#plt.scatter(df.z, df.mb_eff, c = 'seagreen', label = 'High Redshift SN')
plt.errorbar(df.z, df.mb_eff,  df.mb_err, df.z_err, fmt = 'None', c = 'black', alpha = 0.5)#
#lt.plot(x_vals, slope*x_vals+intercept, label = '5.03*$log_{10}(z)$ + 14.73', c = 'black', linestyle = '--')
#plt.fill_between(x_vals, slope*x_vals+(intercept-sd), slope*x_vals+(intercept+sd), alpha = 0.5)
plt.ylabel('m_B Corrected')
plt.xlabel('redshift')
#plt.xscale('log')
#plt.xticks([0.02, 0.05, 0.1, 0.2, 0.5, 1.0])
plt.legend()

#%%
# intialize the "fancy" D
DL = hubble*redshift
print(DL)

#%%
# define necessary functions
# too many functions, I wanted to simplify the process
'''
# define the function to integrate 
def int_func(z_prime, app_m, Om, Ode): # takes in a variable z_prime, and two constants
    function = (np.sqrt(Om*(1+z_prime)**3+Ode))**(-1) 
    return function # oututs the function to to integrate (given in 5)

# define the integrator fucntion
def integral(func, z, Om, Ode): #takes in the function to integrate, the variable, and any constants
    dL = quad((np.sqrt(Om*(1+z)**3+Ode))**(-1) , 0, z)
    return dL #outputs the integrated function

# function to calcuate apparent magnitude from the modified abs magnitude and the calculated value of D_L
def app_mag(arg): # takes in the calculated value of d_L
    return (fancy_M + 5*np.log10(arg)) # outputs the apparent magnitude

# define the curve fit function 
def curve_func(x, a, b): # takes in the variable and any necessary parameters
    return a*np.log(x) + b # returns the expected function

guess = [.5, .5]
def lsq_func(omega, z, m):
    return ((c*(1+z))*quad((np.sqrt(omega[0]*(1+z)**3+omega[1]))**(-1) , 0, z)) - m
   
for i in range(0, len(redshift)):
    z = redshift[i]
    abs_m = mb_both[i]
    matter, darkE =  least_squares(lsq_func, guess, args = (z, abs_m))

'''
#%%


def fit_func(z, Om, Ode):
    #Ode = 1.0 - Om  # Enforce the constraint Om + Ode = 1
    new_mb = []
    func = lambda zprime: (np.sqrt(Om*(1+zprime)**3 + Ode)**(-1))
    for z in redshift:
        integral = quad(func, 0, z)
        int_arg = (c*(1+z)*integral[0])
        new_mb.append(-3.31 + 5*np.log10(int_arg))
    return new_mb
    
#%%
 # initial guess between 0 and 1
# start the lsq fitting
params, covar = curve_fit(fit_func, redshift, app_mag, p0 = [.2, .8])
print(params)


#%%
test_fit = fit_func(redshift, params[0], params[1])


#%%
optA = [0, 1]
optB = [.5, .5]
optC = [1, 0]
optD = [1.5, -0.5]

paramsA, covarA = curve_fit(fit_func, redshift, app_mag, p0 = optA)
paramsB, covarB= curve_fit(fit_func, redshift, app_mag, p0 = optB)
paramsC, covarC = curve_fit(fit_func, redshift, app_mag, p0 =optC)
paramsD, covard = curve_fit(fit_func, redshift, app_mag, p0 =optD)

fitA = fit_func(redshift, optA[0], optA[1])
fitB = fit_func(redshift, optB[0], optB[1])
fitC = fit_func(redshift, optC[0], optC[1])
fitD = fit_func(redshift, optD[0], optD[1])

print(paramsA)
#%%
z_sort = sorted(redshift)
x_vals = np.linspace(min(df1.z), max(df.z), 1000)

def curve(z, a, b):
    return a*np.log10(z)+b
curve_pA, cvarA = curve_fit(curve, redshift, fitA)
curve_pB, cvarB = curve_fit(curve, redshift, fitB)
curve_pC, cvarC = curve_fit(curve, redshift, fitC)
curve_pD, carD = curve_fit(curve, redshift, fitD)
#%%
print(curve_pA, cvarA)
#%%

# plot the apparent magnitude as a funciton of redhsift
# plot the two different data sets with different colours

plt.scatter(df1.z, df1.mb_corr, c = 'firebrick', label = 'Low Redshift SN')
plt.errorbar(df1.z, df1.mb_corr,  df1.mb_err, df1.z_err, fmt = 'None', c = 'black', alpha = 0.5)
plt.scatter(df.z, df.mb_eff, c = 'seagreen', label = 'High Redshift SN')
plt.errorbar(df.z, df.mb_eff,  df.mb_err, df.z_err, fmt = 'None', c = 'black', alpha = 0.5)#
plt.plot(z_sort, sorted(fitA), label = '[0, 1]')
#plt.plot(z_sort, sorted(test_fit), label = 'No Restrictions', c = 'purple', linewidth = 2)
plt.plot(z_sort, sorted(fitB), label = '[.5, .5]')
plt.plot(z_sort, sorted(fitC), label = '[1, 0]')
plt.plot(z_sort, sorted(fitD), label = '[1.5, -.5]')
#lt.plot(x_vals, slope*x_vals+intercept, label = '5.03*$log_{10}(z)$ + 14.73', c = 'black', linestyle = '--')
#plt.fill_between(x_vals, slope*x_vals+(intercept-sd), slope*x_vals+(intercept+sd), alpha = 0.5)
plt.ylabel('m_B Corrected')
plt.xlabel('redshift')
#plt.xscale('log')
#plt.xticks([0.02, 0.05, 0.1, 0.2, 0.5, 1.0])
plt.legend()


#%%

def dist_mod(appm, absM):
    return appm - absM

abs_mag = fancy_M +5*np.log10(hubble) - 25
print(abs_mag)

d_modA = dist_mod(sorted(fitA), abs_mag)
d_modB = dist_mod(sorted(fitB), abs_mag)
d_modC = dist_mod(sorted(fitC), abs_mag)
d_modD = dist_mod(sorted(fitD), abs_mag)



