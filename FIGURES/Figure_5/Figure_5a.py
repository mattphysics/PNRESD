# Generate Figure 5a of Van Zalinge et al. (2017), On determining the point of no return in climate change, Earth System Dynamics.
import numpy as N
import matplotlib.pyplot as plt
from numpy import trapz
import matplotlib as mpl
mpl.rcParams['lines.linewidth'] = 2

#############
#DATA PLASIM#
#############
# Instantaneous doubling
GMST_2CO2=N.load('../../DATA/PlaSim_Ensembles/2co2_ym.npy')
# Smooth increase of 1% per year
GMST_1PR=N.load('../../DATA/PlaSim_Ensembles/co2_1pr_ym.npy')
len_ens=200.

#####################################
#FORCING PROFILE 1 PERCENT EACH YEAR#
#####################################
CO2=N.zeros(len_ens)
for i in N.arange(len_ens):
	if i<=70:
		CO2[i]=1./70.*i
	else:
		CO2[i]=1.
		
#####################################
#EXPECTATION VALUE MINUS BEGIN VALUE#
#####################################
EXP_1pr=N.zeros(len_ens)
EXP_2CO2=N.zeros(len_ens)

for i in N.arange(len_ens):
	EXP_1pr[i]=N.mean(GMST_1PR[i])-287.804719046
	EXP_2CO2[i]=N.mean(GMST_2CO2[i])-287.804719046
EXP_2CO2[0]=0
EXP_1pr[0]=0

################
#GREEN FUNCTION#
################
Green=N.zeros(len_ens)
Green=N.gradient(EXP_2CO2)

# N.save('/Users/user/Thesis_Brenda/Codes/Figure_14/DATA/Green_exp.npy',Green)
#############
#CONVOLUTION#
#############
exp_LRT=N.zeros(len_ens)
for t in N.arange(len_ens):
	f=N.zeros(len_ens)
	for sigma in N.arange(len_ens):
		if t-sigma>=0:
			f[sigma]=CO2[t-sigma]
		else:
			f[sigma]=0
	exp_LRT[t]=trapz(Green*f,dx=1)

######
#PLOT#
######
plt.plot(exp_LRT+287.804719046,label="LRT",color='#0066FF')
plt.plot(EXP_1pr+287.804719046, label="PLASIM",color='#F37C51')

####################
#CUSTOMIZING LAYOUT#
####################
plt.xlabel('Years',fontsize='20')
plt.ylabel('Exp (K)',fontsize='20')
plt.title("Expectation value GMST",fontsize='20')

plt.xticks(fontsize='16')
plt.yticks(fontsize='16')

plt.gcf().subplots_adjust(bottom=0.15)
plt.legend(loc=4, fontsize='16')
plt.show()
