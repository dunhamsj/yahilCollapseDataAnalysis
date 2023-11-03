#!/usr/bin/env python3

import numpy as np
from sys import argv
import matplotlib.pyplot as plt
plt.style.use( '../publication.sty' )

import sys
sys.path.append( '../' )

from UtilitiesModule import GetData

#### ========== User Input ==========

ID = 'YahilCollapse_XCFC'

# Specify directory containing plotfiles
plotfileDirectory \
  = '/home/kkadoogan/Work/Codes/thornado/SandBox/AMReX/Applications/' \
      + '{:}/'.format( ID )

# Specify plot file base name
plotfileBaseName = ID + '.plt'

# Specify field to plot
field = 'PF_D'
yLabel = r'$\rho\ \left[\mathrm{g\,cm}^{-3}\right]$'

# Specify to plot in log-scale
useLogYScale = True

maxLevel = -1

verbose = True

useCustomYLimits = False
eps = 1.0e-5
ymin = 1.0 - 0.1*eps
ymax = 1.0 + 0.5*eps

saveFig = False

#### ====== End of User Input =======

Data, DataUnit, X1_C, X2_C, X3_C, dX1, dX2, dX3, xL, xH, nX, time \
  = GetData( plotfileDirectory, plotfileBaseName, field, \
             'spherical', True, argv = argv, \
             MaxLevel = maxLevel, \
             ReturnTime = True, ReturnMesh = True, Verbose = verbose )

fig, ax = plt.subplots( 1, 1 )
ax.set_title( r'$\texttt{{{:}}}$'.format( ID ), fontsize = 15 )

ax.plot( X1_C[:,0,0], Data[:,0,0], \
         label = r'$t={:.2f}\,\mathrm{{ms}}$'.format( time ) )

if( useLogYScale ): ax.set_yscale( 'log' )

ax.set_xlabel( r'$r\ \left[\mathrm{km}\right]$' )
ax.set_ylabel( yLabel )

ax.set_xlim( 1.0, 2.0e5 )
ax.set_xscale( 'log' )

if useCustomYLimits: ax.set_ylim( ymin, ymax )

ax.legend()

if saveFig:

    figName = '/home/kkadoogan/fig.YahilCollapse_TvA.png'
    plt.savefig( figName, dpi = 300 )
    print( '\n  Saved {:}'.format( figName ) )

else:

    plt.show()

plt.close()

import os
os.system( 'rm -rf __pycache__ ' )
