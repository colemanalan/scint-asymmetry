import numpy as np
import matplotlib.pyplot as plt
import re
import math
import os
from matplotlib.offsetbox import AnchoredText
from scipy.optimize import curve_fit
from matplotlib import rc, rcParams
rc('text', usetex=True)
rc('font', size = 18.0)
rc('font', **{'family': 'serif'})
rc('mathtext', **{'fontset': 'dejavuserif'})
 
def hillas6(T, P1, P2, P3, P4, P5, P6):
    return P1 * np.power( (T-P2)/(P3-P2) , (P3-P2)/(P4+P5*T+P6*T**2)) * np.exp((P3-T)/(P4+P5*T+P6*T**2))

listt  = ['data/prot_lgE_16.0_sin2_0.4/DAT000041.long']

fig, ax = plt.subplots(1,1)
for l in range(0, len(listt)) :
    long_file = listt[l]
    dir_ = os.path.dirname(long_file)
    file_ =os.path.basename(long_file)
    nr_ = file_[3:-5]
    ### getting the zenith angle and energy
    lst_file = dir_+'/SIM'+file_[3:-5] + '.log'
    lst_file_ = open(lst_file, 'r')
    for line in lst_file_:
        if re.search('THETA OF INCIDENCE IS FIXED TO', line):
            zenith = float(line.split(' ')[12])
        if re.search('PRIMARY ENERGY IS FIXED AT', line):
            en = math.log10(float(line.split(' ')[12+5])  )
            en = 10**(en-6)
    long_file_ = open(listt[l], 'r')
    ### --------------------------------------------
    for line in long_file_:
        if re.search(' LONGITUDINAL DISTRIBUTION IN', line):
            try: 
                nmax = float(line.split(' ' )[6])
            except:
                nmax = float(line.split(' ' )[7])
        if re.search('PARAMETERS         =', line):
            pp = list(filter(None, line.split(" ")[-10:]))
            par1 = [float(pp[0]), float(pp[1]), float(pp[2]), float(pp[3]), float(pp[4]), float(pp[5])]
    
    data = np.genfromtxt(long_file, skip_header=2, max_rows=nmax, names=['DEPTH', 'GAMMAS', 'POSITRONS', 'ELECTRONS', 'MUp', 'MUm', 'HADRONS', 'CHARGED', 'NUCLEI', 'CHERENKOV'])
    rangeP = np.arange(0., 1000., 0.01)
    plt.scatter(data['DEPTH'],  data['CHARGED'], label="charged", c="#659ece")
    plt.plot(rangeP, hillas6(rangeP, *par1), color="#3e6c93", label="Gaisser-Hillas", lw=0.8)
    
    text = r'\noindent H,'+r' E$_{\mathrm{MC}}$ = '+str(round(en,2))+' PeV\n'+r'$\theta$ = '+str(round(zenith,2))+r'$^{\circ}$'+'\n X$_{\mathrm{max}}$ = '+str(round(par1[2], 1))+' g/cm$^2$'
    anchored_text = AnchoredText( text, loc=2, prop=dict(fontsize=10), frameon=False) 
    ax.add_artist(anchored_text)
    del data

plt.ticklabel_format(axis='y',style='sci',scilimits=(1,4))
plt.xlabel(r'Slant depth/ g\,cm$^{-2}$')
plt.ylabel(r'Number of charged particles')  
plt.legend( loc='upper right', fontsize=12)
plt.tight_layout()
fig.savefig('plots/LongH_'+str(nr_)+'_'+str(round(en,2))+'.png', dpi=400)



