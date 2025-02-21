
#standard imports and setups used during EAE127
import math
import numpy as np 
import os
import matplotlib.pyplot as plt 
import matplotlib.lines as mlines
#Disable Python Warning Output
#(NOTE: Only for production, comment out for debugging)
import warnings
warnings.filterwarnings('ignore')
### PLOTTING DEFAULTS BOILERPLATE (OPTIONAL) #########################
#SET DEFAULT FIGURE APPERANCE
import seaborn as sns #Fancy plotting package 
#No Background fill, legend font scale, frame on legend
sns.set_theme(style='whitegrid', font_scale=1.5, rc={'legend.frameon': True})
#Mark ticks with border on all four sides (overrides 'whitegrid')
sns.set_style('ticks')
#ticks point in
sns.set_style({"xtick.direction": "in","ytick.direction": "in"})
#fix invisible marker bug
sns.set_context(rc={'lines.markeredgewidth': 0.1})
#restore default matplotlib colormap
mplcolors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd',
'#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']
sns.set_palette(mplcolors)

#Get color cycle for manual colors
colors = sns.color_palette()
#SET MATPLOTLIB DEFAULTS
#(call after seaborn, which changes some defaults)
params = {
#FONT SIZES
'axes.labelsize' : 30, #Axis Labels
'axes.titlesize' : 30, #Title
'font.size' : 28, #Textbox
'xtick.labelsize': 22, #Axis tick labels
'ytick.labelsize': 22, #Axis tick labels
'legend.fontsize': 15, #Legend font size
'font.family' : 'serif',
'font.fantasy' : 'xkcd',
'font.sans-serif': 'Helvetica',
'font.monospace' : 'Courier',
#AXIS PROPERTIES
'axes.titlepad' : 2*6.0, #title spacing from axis
'axes.grid' : True, #grid on plot
'figure.figsize' : (8,8), #square plots
'savefig.bbox' : 'tight', #reduce whitespace in saved figures
#LEGEND PROPERTIES
'legend.framealpha' : 0.5,
'legend.fancybox' : True,
'legend.frameon' : True,
'legend.numpoints' : 1,
'legend.scatterpoints' : 1,
'legend.borderpad' : 0.1,
'legend.borderaxespad' : 0.1,
'legend.handletextpad' : 0.2,
'legend.handlelength' : 1.0,
'legend.labelspacing' : 0,
}
import matplotlib #type:ignore
matplotlib.rcParams.update(params) #update matplotlib defaults, call after￿
### END OF BOILERPLATE ##################################################
colors = sns.color_palette() #color cycle


k = 1.4
Cp = 0.24 #BTU/(lbm * Rankine)
R = 53.35 #ft * lbf / (lbm * Rankine)

mdot = 192 #lbm/s
deltaH = 17900 #Btu/lbm
Tt4 = 2350 #Rankine

Pa = 6.20027 #ambient pressure
Tstp = 440.248 #ambient temperature

def getThrust (piC = 17, statTempRatio = 1,Tt4 = 2350, Ma = 0.88):
    
    Ta = statTempRatio * Tstp

    #initial calculations
    aa = (k*R*Ta*32.17)**(1/2)
    ua = Ma*aa

    #compressor
    Pt2 = Pa * (1 + ((k-1)/2)*Ma**2)**(k/(k-1))
    Pt3 = piC*Pt2
    Tt2 = Ta * (1 + ((k-1)/2)*Ma**2)

    tauC = piC**((k-1)/k)
    Tt3 = tauC * Tt2

    #burner
    mdotf = mdot * Cp * (Tt4 - Tt3) / deltaH

    #turbine
    Pt4 = Pt3 #ideal assumption
    Tt5 = Tt4 - (Tt3-Tt2)

    tauT = Tt5/Tt4
    piT = tauT**(k/(k-1))

    Pt5 = Pt4 * piT
    #nozzle
    Pt8 = Pt5
    P8 = Pa
    M8 = (((Pt8/P8)**((k-1)/k) - 1) * (2/(k-1)))**(1/2)

    Tt8 = Tt5
    T8 = Tt8 * (1 + ((k-1)/2)*M8**2)**(-1)

    a8 = (k*R*T8*32.17)**(1/2)
    u8 = M8*a8

    F = mdot*(u8 - ua) / 32.17
    TSFC = (mdotf/F) * 3600 
    ndTSFC = (TSFC / 3600) * 32.17 * aa
    specificF = F/(mdot*aa) * 32.17

    return F, ndTSFC, specificF


n = 200
x = []
STarr = []
ndTSFCarr = []
for i in np.linspace(5,20,n):
    x.append(i)
    F, ndTSFC, specificF = getThrust(piC=i)
    STarr.append(specificF)
    ndTSFCarr.append(ndTSFC)

plt.figure()
plt.plot(x,STarr)
plt.xlabel("Compressor Ratio")
plt.ylabel("Non-dimensional Thrust")

plt.figure()
plt.plot(x,ndTSFCarr)
plt.xlabel("Compressor Ratio")
plt.ylabel("Non-dimensional TSFC (x10^-3)")
plt.show()

x = []
STarr = []
ndTSFCarr = []
for i in np.linspace(0.5,1.2,n):
    x.append(i)
    F, ndTSFC, specificF = getThrust(statTempRatio = i)
    STarr.append(specificF)
    ndTSFCarr.append(ndTSFC)

plt.figure()
plt.plot(x,STarr)
plt.xlabel("Static Temperature Ratio")
plt.ylabel("Non-dimensional Thrust")

plt.figure()
plt.plot(x,ndTSFCarr)
plt.xlabel("Static Temperature Ratio")
plt.ylabel("Non-dimensional TSFC (x10^-3)")
plt.show()

x = []
STarr = []
ndTSFCarr = []
for i in np.linspace(4,8,n):
    x.append(i)
    Tt4in = Tstp * i
    F, ndTSFC, specificF = getThrust(Tt4=Tt4in)
    STarr.append(specificF)
    ndTSFCarr.append(ndTSFC)

plt.figure()
plt.plot(x,STarr)
plt.xlabel("Burner Temperature Ratio")
plt.ylabel("Non-dimensional Thrust")

plt.figure()
plt.plot(x,ndTSFCarr)
plt.xlabel("Burner Temperature Ratio")
plt.ylabel("Non-dimensional TSFC (x10^-3)")
plt.show()

x = []
STarr = []
ndTSFCarr = []
for i in np.linspace(0,5,n):
    x.append(i)
    F, ndTSFC, specificF = getThrust(Ma=i)
    STarr.append(specificF)
    ndTSFCarr.append(ndTSFC)

plt.figure()
plt.plot(x,STarr)
plt.xlabel("Mach Number")
plt.ylabel("Non-dimensional Thrust")

plt.figure()
plt.plot(x,ndTSFCarr)
plt.xlabel("Mach Number")
plt.ylabel("Non-dimensional TSFC (x10^-3)")
plt.show()