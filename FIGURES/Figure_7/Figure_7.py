# Generate Figure 7a, 7b of Van Zalinge et al. (2017), On determining the point of no return in climate change, Earth System Dynamics.
# ============================================================
# PACKAGE IMPORT
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

from datetime import datetime,timedelta
from scipy.optimize import leastsq
from scipy.stats import norm

# Graphics functions =========================================
def figsize(scale):
    '''Give fraction of width, get figsize in inches with height scaled according to golden_mean'''
    fig_width_pt = 426.79135 # Get this from LaTeX using \the\textwidth
    inches_per_pt = 1.0/72.27                       # Convert pt to inch
    golden_mean = (np.sqrt(5.0)-1.0)/2.0            # Aesthetic ratio (you could change this)
    fig_width = fig_width_pt*inches_per_pt*scale    # width in inches
    fig_height = fig_width*golden_mean              # height in inches
    fig_size = [fig_width,fig_height]
    return fig_size

def figsize2(width,height):
    '''Give fraction of width and heigt, get figsize in inches'''
    textwidth_pt = 426.79135 
    textheight_pt = 426.79135 
    inches_per_pt = 1.0/72.27 
    return textwidth_pt * inches_per_pt * width,textheight_pt * inches_per_pt * height

# I make my own newfig and savefig functions
def newfig(width,nr=1,nc=1,**kwargs):
    '''use function figsize'''
    fig,ax = plt.subplots(nrows=nr,ncols=nc,figsize=figsize(width),**kwargs)
    return fig, ax

def newfig2(width,height,nr=1,nc=1,**kwargs):
    '''use function figsize2'''
    fig,ax = plt.subplots(nrows=nr,ncols=nc,figsize=figsize2(width,height),**kwargs)
    return fig, ax

def savefig(fig,filename):
    fig.savefig('{}.pgf'.format(filename),bbox_inches='tight',dpi=400,rasterized=True)
    fig.savefig('{}.pdf'.format(filename),bbox_inches='tight',dpi=400,rasterized=True)

# ============================================================
# Further Functions

# Mitigation scenarios
def mitscen(emis,tstart,efold):
	'''return emissions under a reduction scenario that starts at t=tstart and decays expoentially with efolding timescale efold
	assumes timestep of Delta t = 1 yr and tstart is counted from the start with t_start = 0
	'''
	emisreduce = np.zeros(emis.shape)
	index = np.arange(emis.shape[0])
	emisreduce[:] = emis
	emisreduce[index>=tstart] = emis[tstart] * np.exp(-(index[tstart:]-tstart)/efold)

	return emisreduce

# Radiative forcing
def getRadForcCO2(C):
	'''radiative forcing due to CO2 relative to 1750'''
	alpha = 5.35
	C0 = 278 # ppm # pre-industrial CO2 (1750)
	return alpha * np.log(C/C0)

# Compute temperature response
def getDeltaTCO2_v1(EMIS_CO2_GtC,CO2_0,G_CO2,G_Temp,radForc_factor):
	'''
	compute Temperature response from emissions
	radiative forcing is scaled up by a factor of 		'radForc_factor' 
	'''
	MA = 28.97 # kg kmol^{-1} # mean molecular weight of air
	MC = 12.011 # kg kmol^{-1} # mean molecular weight of C
	TM = 5.1352e18 # kg # total mass of atmosphere
	factorC = (1/MC) * (MA/TM) * 1e18 # conversion from Gt to ppm
	EMIS_CO2_ppm = factorC * EMIS_CO2_GtC

	# Compute CO2 concentration as function of time
	DeltaC_CO2_ppm = np.convolve(G_CO2,EMIS_CO2_ppm,mode='full')[:EMIS_CO2_ppm.shape[0]]
	C_CO2_ppm = CO2_0 + DeltaC_CO2_ppm
	# Compute radiative forcing relative to 1750
	radForc_CO2 = getRadForcCO2(C_CO2_ppm) * radForc_factor

	# Compute temperature perturbation
	DeltaT_CO2 = np.convolve(G_Temp,radForc_CO2,mode='full')[:radForc_CO2.shape[0]]

	res = {
		'EMIS_GtC':EMIS_CO2_GtC,
		'EMIS_ppm':EMIS_CO2_ppm,
		'C_CO2_ppm':C_CO2_ppm,
		'radForc_CO2':radForc_CO2,
		'DeltaT_CO2':DeltaT_CO2
		}
	return res

# Joos et al., 2013 IRF for emissions to concentrations
def IRFexpFit(a0,a1,a2,a3,tau1,tau2,tau3,time):
    '''exponential fit from Joos et al., 2011'''
    return a0 + a1 * np.exp(-time/tau1) + a2 * np.exp(-time/tau2) + a3 * np.exp(-time/tau3)

# Fit exponential
def fitGfunc4(time,a1,t1):
	return a1 * np.exp(-time/t1)

def devGGfit(params,time,Garr,Gfunc):
	''' deviation between fit and data'''
	return Gfunc(time,*params) - Garr

# ============================================================
# ============================================================
# DATA

# Load RCP Emission data
emisrcp26 = pd.read_table('../../DATA/RCP_Data/RCP3PD_EMISSIONS.DAT',skiprows=37,delim_whitespace=True,index_col=0,parse_dates=True,infer_datetime_format=True)
emisrcp26.index = pd.period_range('%i-01-01'%emisrcp26.index[0],'%i-01-01' % emisrcp26.index[-1],freq='a')

emisrcp45 = pd.read_table('../../DATA/RCP_Data/RCP45_EMISSIONS.DAT',skiprows=37,delim_whitespace=True,index_col=0,parse_dates=True,infer_datetime_format=True)
emisrcp45.index = pd.period_range('%i-01-01'%emisrcp45.index[0],'%i-01-01' % emisrcp45.index[-1],freq='a')

emisrcp60 = pd.read_table('../../DATA/RCP_Data/RCP6_EMISSIONS.DAT',skiprows=37,delim_whitespace=True,index_col=0,parse_dates=True,infer_datetime_format=True)
emisrcp60.index = pd.period_range('%i-01-01'%emisrcp60.index[0],'%i-01-01' % emisrcp60.index[-1],freq='a')

emisrcp85 = pd.read_table('../../DATA/RCP_Data/RCP85_EMISSIONS.DAT',skiprows=37,delim_whitespace=True,index_col=0,parse_dates=True,infer_datetime_format=True)
emisrcp85.index = pd.period_range('%i-01-01'%emisrcp85.index[0],'%i-01-01' % emisrcp85.index[-1],freq='a')

# ============================================================
# GREEN FUNCTIONS

# CARBON GREEN FUNCTION FROM JOOS ET AL., 2013
a0,a1,a2,a3,tau1,tau2,tau3=2.17278e-01, 2.24037e-01, 2.82381e-01, 2.76303e-01, 3.94409e+02, 3.65393e+01, 4.30365e+00
time1000 = np.arange(1000)
# For mean response
G_CO2 = IRFexpFit(a0,a1,a2,a3,tau1,tau2,tau3,time1000)

# Variance from Joos Data
pd100 = np.loadtxt('../../DATA/Carbon_Response/IRF_PD100_SMOOTHED_CO2.dat')
G_CO2_var = pd100[:1000,18]**2

# PLASIM Data
# add pre-industrial reference before perturbation at 0-th timestep for each member
T0 = 287.804719046 # K # reference temperature, pre-industrial
Tvar0 = 0.00843298290214 # K^2 # pre-industrial temperature variance
GMST_2co2 = np.zeros((201,201))
GMST_2co2[1:,:] = np.load('../../DATA/PlaSim_Ensembles/2co2_ym.npy')
GMST_2co2[0,:] = T0
GMST_1pct = np.zeros((201,201))
GMST_1pct[1:,:] = np.load('../../DATA/PlaSim_Ensembles/co2_1pr_ym.npy')
GMST_1pct[0,:] = T0

# ensemble mean
EXP_2co2 = GMST_2co2.mean(axis=1)
EXP_1pct = GMST_1pct.mean(axis=1)

# ensemble standard deviation
VAR_2co2 = GMST_2co2.var(axis=1)
VAR_1pct = GMST_1pct.var(axis=1)

VAR_2co2[0] = Tvar0
VAR_1pct[0] = Tvar0

# CO2 Forcing in ppm for abrupt and smooth scenarios
C0 = 278 # ppm
CO2_2co2 = C0 * np.ones(EXP_2co2.shape[0])
CO2_2co2[1:] = 2 * C0

CO2_1pct = C0 * 1.01**(np.arange(EXP_1pct.shape[0]))
CO2_1pct[CO2_1pct>2*C0] = 2*C0

# Radiative forcing in W m^-2 from CO2 Forcing
F2co2 = getRadForcCO2(CO2_2co2)
F1pct = getRadForcCO2(CO2_1pct)

# Compute Data based green functions
Garrd = 1/F2co2[1] * np.diff(EXP_2co2)
Garrd = np.array(list(Garrd)+[Garrd[-1]])

Garr_vard = 1/F2co2[1] * np.diff(VAR_2co2)
Garr_vard = np.array(list(Garr_vard)+[Garr_vard[-1]])

# For long simulations (longer than the simulations on which the Green function is based) we need to extend the Green function. We do so by exponential fits to the data and thereby generate response function for simulations of up to 500 years duration. The data-based functions are used for as long as possible

# Green function for mean response
time = np.arange(Garrd.shape[0])
par7 = leastsq(devGGfit,[1,2],args=(time[1:],Garrd[1:],fitGfunc4))[0]

time500 = np.arange(500)
G_Temp_i = fitGfunc4(time500,*par7)
G_Temp_i[0] = Garrd[0] 

# Assign data-based Green function for t<=200 years
G_Temp = np.copy(G_Temp_i)
G_Temp[:Garrd.shape[0]] = Garrd

# Green function for variance of response
# Variance: Simple white noise with mean and variance from data-based G
time500 = np.arange(500)
Garr_var = np.random.normal(loc=Garr_vard.mean(),scale=Garr_vard.std(),size=len(time500))
Garr_var[:Garr_vard.shape[0]] = Garr_vard

# ============================================================
# SIMULATION RESULTS

# Compute temperature distributions for mitigation scenarios, based on RCP emissions, with exponential emission reduction starting in a given year

CO2_0 = 278.05158 # ppm # pre-industrial CO2 concentration

# RCP2.6
muT2100 = pd.Series(index=np.arange(2005,2099).astype('str'))
varT200_Tonly = pd.Series(index=np.arange(2005,2099).astype('str'))

for tstart in range(2005,2100):
	emisco226 = mitscen(emisrcp26['FossilCO2'][:'2100'].values,tstart-1765,25)
	res = getDeltaTCO2_v1(emisco226,CO2_0,G_CO2,G_Temp,0.6)
	muT2100[str(tstart)] = res['DeltaT_CO2'][-1]

	varT200_Tonly[str(tstart)] = Tvar0 + np.convolve(Garr_var,res['radForc_CO2'],mode='full')[:res['radForc_CO2'].shape[0]][-1]

T0926 = norm.ppf(0.9,muT2100,varT200_Tonly**0.5) 

# RCP4.5
muT2100 = pd.Series(index=np.arange(2005,2099).astype('str'))
varT200_Tonly = pd.Series(index=np.arange(2005,2099).astype('str'))

for tstart in range(2005,2100):
	emisco245 = mitscen(emisrcp45['FossilCO2'][:'2100'].values,tstart-1765,25)
	res = getDeltaTCO2_v1(emisco245,CO2_0,G_CO2,G_Temp,0.6)
	muT2100[str(tstart)] = res['DeltaT_CO2'][-1]

	varT200_Tonly[str(tstart)] = Tvar0 + np.convolve(Garr_var,res['radForc_CO2'],mode='full')[:res['radForc_CO2'].shape[0]][-1]

T0945 = norm.ppf(0.9,muT2100,varT200_Tonly**0.5)

# RCP6.0
muT2100 = pd.Series(index=np.arange(2005,2099).astype('str'))
varT200_Tonly = pd.Series(index=np.arange(2005,2099).astype('str'))

for tstart in range(2005,2100):
	emisco260 = mitscen(emisrcp60['FossilCO2'][:'2100'].values,tstart-1765,25)
	res = getDeltaTCO2_v1(emisco260,CO2_0,G_CO2,G_Temp,0.6)
	muT2100[str(tstart)] = res['DeltaT_CO2'][-1]

	varT200_Tonly[str(tstart)] = Tvar0 + np.convolve(Garr_var,res['radForc_CO2'],mode='full')[:res['radForc_CO2'].shape[0]][-1]

T0960 = norm.ppf(0.9,muT2100,varT200_Tonly**0.5)

# RCP8.5
muT2100 = pd.Series(index=np.arange(2005,2099).astype('str'))
varT200_Tonly = pd.Series(index=np.arange(2005,2099).astype('str'))

for tstart in range(2005,2100):
	emisco285 = mitscen(emisrcp85['FossilCO2'][:'2100'].values,tstart-1765,25)
	res = getDeltaTCO2_v1(emisco285,CO2_0,G_CO2,G_Temp,0.6)
	muT2100[str(tstart)] = res['DeltaT_CO2'][-1]

	varT200_Tonly[str(tstart)] = Tvar0 + np.convolve(Garr_var,res['radForc_CO2'],mode='full')[:res['radForc_CO2'].shape[0]][-1]

T0985 = norm.ppf(0.9,muT2100,varT200_Tonly**0.5) 


# Point of no return
T0926 = pd.Series(T0926,index=muT2100.index)
T0945 = pd.Series(T0945,index=muT2100.index)
T0960 = pd.Series(T0960,index=muT2100.index)
T0985 = pd.Series(T0985,index=muT2100.index)


# ============================================================
# PLOTTING
# ============================================================

# FIGURE 7a

fig,ax = newfig(1.0)
T0926.plot(ax=ax,label='RCP2.6')
T0945.plot(ax=ax,label='RCP4.5')
T0960.plot(ax=ax,label='RCP6.0')
T0985.plot(ax=ax,label='RCP8.5')
pd.Series(2*np.ones(T0985.shape),index=muT2100.index).plot(c='k',linestyle='--',label='')
ax.set_xlabel('year starting to reduce emissions')
ax.set_ylabel(r'$90\%$ Warming in 2100 relative to pre-industrial')
ax.set_title(r'PLASIM warming for exponential emissions decrease (25 yr)')
ax.legend()
#savefig(fig,'PLASIM_DeltaT_emisreduce6_2')
plt.show()


# ============================================

#FIGURE 7b

#PofnoR26 = T0926.index[np.where(T0926>=2)[0][0]]
PofnoR45 = int(T0945.index[np.where(T0945>=2)[0][0]])
PofnoR60 = int(T0960.index[np.where(T0960>=2)[0][0]])
PofnoR85 = int(T0985.index[np.where(T0985>=2)[0][0]])

# Plot concentrations for mitigation starting at point of no return
fig,ax = newfig(1.0)
years = np.arange(1765,2101)

# RCP2.6
resn26 = getDeltaTCO2_v1(emisrcp26['FossilCO2'][:'2100'].values,CO2_0,G_CO2,G_Temp,0.6)
ax.plot(years,resn26['C_CO2_ppm'],c='C0')

# RCP4.5
resn45 = getDeltaTCO2_v1(emisrcp45['FossilCO2'][:'2100'].values,CO2_0,G_CO2,G_Temp,0.6)
ax.plot(years,resn45['C_CO2_ppm'],c='C1')

emisco245 = mitscen(emisrcp45['FossilCO2'][:'2100'].values,PofnoR45-1765,25)
resy45 = getDeltaTCO2_v1(emisco245,CO2_0,G_CO2,G_Temp,0.6)
ax.plot(years,resy45['C_CO2_ppm'],'--',c='C1')

# RCP6.0
resn60 = getDeltaTCO2_v1(emisrcp60['FossilCO2'][:'2100'].values,CO2_0,G_CO2,G_Temp,0.6)
ax.plot(years,resn60['C_CO2_ppm'],c='C2')

emisco260 = mitscen(emisrcp60['FossilCO2'][:'2100'].values,PofnoR60-1765,25)
resy60 = getDeltaTCO2_v1(emisco260,CO2_0,G_CO2,G_Temp,0.6)
ax.plot(years,resy60['C_CO2_ppm'],'--',c='C2')

# RCP8.5
resn85 = getDeltaTCO2_v1(emisrcp85['FossilCO2'][:'2100'].values,CO2_0,G_CO2,G_Temp,0.6)
ax.plot(years,resn85['C_CO2_ppm'],c='C3')

emisco285 = mitscen(emisrcp85['FossilCO2'][:'2100'].values,PofnoR85-1765,25)
resy85 = getDeltaTCO2_v1(emisco285,CO2_0,G_CO2,G_Temp,0.6)
ax.plot(years,resy85['C_CO2_ppm'],'--',c='C3')

ax.plot(PofnoR45,resy45['C_CO2_ppm'][np.where(years==PofnoR45)],'o',label=r'RCP4.5, $\pi_t = %i$' % PofnoR45,c='C1')
ax.plot(PofnoR60,resy60['C_CO2_ppm'][np.where(years==PofnoR60)],'o',label=r'RCP6.0, $\pi_t = %i$' % PofnoR60,c='C2')
ax.plot(PofnoR85,resy85['C_CO2_ppm'][np.where(years==PofnoR85)],'o',label=r'RCP8.5, $\pi_t = %i$' % PofnoR85,c='C3')

ax.set_xlim(2000,2100)
ax.set_ylim(350,600)
ax.ticklabel_format(axis='x',scilimits=(-5,5))  
ax.legend()
ax.set_xlabel('Year')
ax.set_ylabel(r'CO$_2$ (ppm)')

#savefig(fig,'PLASIM_PofnoR_CO2')
plt.show()

# ============================================================

# End of script
