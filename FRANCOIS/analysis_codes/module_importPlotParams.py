## !!! XXX : you can comment the line     matplotlib.use('Agg')   in order to have interactive plots.
## if you run on a cluster, you must use  matplotlib.use('Agg')   to desactivate X server use.

import matplotlib
#matplotlib.use('Agg')   ## use a non-interactive backend (allows saving saveFigFormat files) (it does not use the standard screen output (X server))
import matplotlib.pyplot as plt 


def fsaveFigFormat():
    saveFigFormat='.svg'
    saveFigFormat='.pdf'
    saveFigFormat='.png'
    return saveFigFormat

def savefig_perso(name):
    plt.tight_layout()
    plt.savefig(name)
    return 0

def fparams():
    params = {'backend': 'pdf',
          'axes.labelsize': 25,
          'axes.labelweight': 'roman',
          'axes.titlesize': 16,
          'legend.fontsize': 16,
          'xtick.labelsize': 18,
          'ytick.labelsize': 18,
          'font.size': 12,
          'text.usetex': False,
#          'font.family': 'serif',
          'font.serif': ['computer modern roman'],
            }
    return params
