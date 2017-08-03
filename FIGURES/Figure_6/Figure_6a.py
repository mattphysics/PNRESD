# Generate Figure 6a of Van Zalinge et al. (2017), On determining the point of no return in climate change, Earth System Dynamics.
import numpy as N
import matplotlib.pyplot as plt
import matplotlib.colors
import matplotlib as mpl
mpl.rcParams['lines.linewidth'] = 2.

#########################
#THE POINT OF NO RETURNS#
#########################
PofnoR45=2055
PofnoR60=2057
PofnoR85=2047
PofnoR26=2097

###############
#RCP SCENARIOS#
###############
file = open("../../DATA/RCP_Data/RCP45.dat","r")
lines = file.readlines()
RCP45=N.array(lines,'f')

file = open("../../DATA/RCP_Data/RCP85.dat","r")
lines = file.readlines()
RCP85=N.array(lines,'f')

file = open("../../DATA/RCP_Data/RCP60.dat","r")
lines = file.readlines()
RCP60=N.array(lines,'f')

file = open("../../DATA/RCP_Data/RCP26.dat","r")
lines = file.readlines()
RCP26=N.array(lines,'f')

###################
#EXPONENTIAL DECAY#
###################
def f(t,C_0):
	y=(C_0-400.)*N.exp(-t/25.)+400.
	return y

tt=N.arange(2000,2200)

############
#PLOT RCP85#
############
plt.plot(tt,RCP85,color='#9D0444')

plt.plot(tt[PofnoR85-2000:],f(N.arange(len(tt[PofnoR85-2000:])),RCP85[PofnoR85-2000]),color='#9D0444',linestyle=':')
plt.scatter(tt[PofnoR85-2000:][0],RCP85[PofnoR85-2000], marker='o',s=90, color='#9D0444',linestyle=':',label='RCP8.5   $\pi_t=2047$')

############
#PLOT RCP60#
############
plt.plot(tt,RCP60,color='#F37C51')

plt.plot(tt[PofnoR60-2000:],f(N.arange(len(tt[PofnoR60-2000:])),RCP60[PofnoR60-2000]),color='#F37C51',linestyle=':')
plt.scatter(tt[PofnoR60-2000:][0],RCP60[PofnoR60-2000], marker='o',s=90, color='#F37C51',linestyle=':',label='RCP6.0   $\pi_t=2057$')

############
#PLOT RCP45#
############
plt.plot(tt,RCP45,color='#0066FF')

plt.plot(tt[PofnoR45-2000:],f(N.arange(len(tt[PofnoR45-2000:])),RCP45[PofnoR45-2000]),color='#0066FF',linestyle=':')
plt.scatter(tt[PofnoR45-2000:][0],RCP45[PofnoR45-2000], marker='o',s=90, color='#0066FF',linestyle=':',label='RCP4.5   $\pi_t=2055$')


############
#PLOT RCP26#
############
plt.plot(tt,RCP26,color='#148F77')

plt.plot(tt[PofnoR26-2000:],f(N.arange(len(tt[PofnoR26-2000:])),RCP26[PofnoR26-2000]),color='#148F77',linestyle=':')
plt.scatter(tt[PofnoR26-2000:][0],RCP26[PofnoR26-2000], marker='o',s=90, color='#148F77',linestyle=':',label='RCP2.6   $\pi_t=2097$')

####################
#CUSTOMIZING LAYOUT#
####################
plt.xticks(fontsize='16')
plt.yticks(fontsize='16')
plt.xlabel("Year", fontsize=20)
plt.ylabel('CO$_2$eq (ppmv)',fontsize='20')

plt.gcf().subplots_adjust(bottom=0.13)

#plt.legend(scatterpoints=1, loc=1, fontsize='16')
plt.legend(scatterpoints=1, fontsize='16', loc=3, ncol=2)#, mode='expand', scatterpoints=1, fontsize='16')
plt.axis([2000,2150,250,750])
plt.show()
