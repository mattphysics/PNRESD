# Generate Figure 6b of Van Zalinge et al. (2017), On determining the point of no return in climate change, Earth System Dynamics.
import numpy as N
import matplotlib.pyplot as plt
import matplotlib.colors
import matplotlib as mpl
mpl.rcParams['lines.linewidth'] = 2.

#########################
#THE POINT OF NO RETURNS#
#########################
PofnoR08=2057
PofnoR09=2055
PofnoR099=2051

###############
#RCP SCENARIOS#
###############
file = open("../../DATA/RCP_Data/RCP45.dat","r")
lines = file.readlines()
RCP45=N.array(lines,'f')

###################
#EXPONENTIAL DECAY#
###################
def f(t,C_0):
	y=(C_0-400.)*N.exp(-t/25.)+400.
	return y

tt=N.arange(2000,2200)

####################
#PROBABILITY IS 0.8#
####################
plt.plot(tt,RCP45[:200],color='#0066FF')

plt.plot(tt[PofnoR08-2000:],f(N.arange(len(tt[PofnoR08-2000:])),RCP45[PofnoR08-2000]),color='#148F77',linestyle=':')
plt.scatter(tt[PofnoR08-2000:][0],RCP45[PofnoR08-2000], marker='o',s=90, color='#148F77',label='$\\beta_T=0.80$    $\pi_t=2057$')

####################
#PROBABILITY IS 0.9#
####################
plt.plot(tt[PofnoR09-2000:],f(N.arange(len(tt[PofnoR09-2000:])),RCP45[PofnoR09-2000]),color='#F37C51',linestyle=':')
plt.scatter(tt[PofnoR09-2000:][0],RCP45[PofnoR09-2000], marker='o',s=90, color='#F37C51',label='$\\beta_T=0.90$    $\pi_t=2055$')


#####################
#PROBABILITY IS 0.99#
#####################
plt.plot(tt[PofnoR099-2000:],f(N.arange(len(tt[PofnoR099-2000:])),RCP45[PofnoR099-2000]),color='#9D0444',linestyle=':')
plt.scatter(tt[PofnoR099-2000:][0],RCP45[PofnoR099-2000], marker='o',s=90, color='#9D0444',label='$\\beta_T=0.99$    $\pi_t=2051$')



####################
#CUSTOMIZING LAYOUT#
####################
plt.xticks(fontsize='16')
plt.yticks(fontsize='16')

plt.xlabel("Year", fontsize=20)
plt.ylabel('CO$_2$eq (ppmv)',fontsize='20')

plt.gcf().subplots_adjust(bottom=0.13)
plt.legend(scatterpoints=1,loc=4, fontsize='16')
plt.axis([2000,2200,250,600])

plt.show()
