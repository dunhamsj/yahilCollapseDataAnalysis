#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
plt.style.use( '../publication.sty' )

import sys
sys.path.append( '../' )

from UtilitiesModuleHDF import ReadFieldsHDF

field = 'PF_D'
yLabel = r'$\rho\,\left[\mathrm{g\,cm}^{-3}\right]$'
yScale = 1.0
ylim = []
useLogYScale = True

#field = 'PolytropicConstant'
#yLabel = r'$K/K_{\mathrm{central}}$'
#yScale = 6.0e27 / ( 7.0e9 )**( 1.30 )
#epsMin = 1.0e-5
#epsMax = 1.0e-5
#ylim = [ 1.0 - epsMin, 1.0 + epsMax ]
#useLogYScale = False

def ComputeCellAverage( names, iSS, field ):

    rC = names['X1_C'][1]
    nX = rC.shape[0]
    del rC

    r  = names['X1'][1]
    nN = r.shape[0] // nX
    del r

    if nN == 2:
        WeightsX = np.array( [ 0.5, 0.5 ], np.float64 )

    SqrtGm = names['GF_Sg'][1][iSS,0,0,:]

    uK = np.empty( nX, np.float64 )

    if field == 'PolytropicConstant':

        p   = names['AF_P' ][1][iSS,0,0,:]
        rho = names['PF_D' ][1][iSS,0,0,:]
        Gm  = names['AF_Gm'][1][iSS,0,0,:]

        pK   = np.empty( nX, np.float64 )
        rhoK = np.empty( nX, np.float64 )
        GmK  = np.empty( nX, np.float64 )

        for iX1 in range( nX ):

            iLo = nN * iX1
            iHi = iLo + nN

            pK[iX1] \
              = np.sum( WeightsX * p[iLo:iHi] * SqrtGm[iLo:iHi] ) \
                  / np.sum( WeightsX * SqrtGm[iLo:iHi] )
            rhoK[iX1] \
              = np.sum( WeightsX * rho[iLo:iHi] * SqrtGm[iLo:iHi] ) \
                  / np.sum( WeightsX * SqrtGm[iLo:iHi] )
            GmK[iX1] \
              = np.sum( WeightsX * Gm[iLo:iHi] * SqrtGm[iLo:iHi] ) \
                  / np.sum( WeightsX * SqrtGm[iLo:iHi] )

        uK = pK / rhoK**GmK

    else:

        uN = names[field][1][iSS,0,0,:]

        for iX1 in range( nX ):

            iLo = nN * iX1
            iHi = iLo + nN

            uK[iX1] \
              = np.sum( WeightsX * uN[iLo:iHi] * SqrtGm[iLo:iHi] ) \
                  / np.sum( WeightsX * SqrtGm[iLo:iHi] )

    return uK
# END ComputeCellAverage( names, iSS, field )

THORNADO_DIR = HOME + 'Work/Codes/thornado/'

plotfileDirectory \
  = THORNADO_DIR + 'SandBox/YahilCollapse_XCFC/Output/YahilCollapse'

snapshots = [ 0, 1100 ]

names \
  = ReadFieldsHDF \
      ( plotfileDirectory, \
        snapshots, \
        CoordinateSystem = 'SPHERICAL', \
        UsePhysicalUnits = True )

rC   = names['X1_C'][1]
Time = names['Time'][1]

fig, ax = plt.subplots( 1, 1 )

rhoK0 = ComputeCellAverage( names, 0, field )
rhoK1 = ComputeCellAverage( names, -1, field )

ax.set_title( 'Yahil Collapse, Native thornado/Poseidon, Piecewise-Uniform Mesh' )
ax.semilogx( rC, rhoK0 / yScale, '.' )
ax.semilogx( rC, rhoK1 / yScale, '.' )
#np.savetxt( 'Yahil_native.dat', np.vstack( ( rC, rhoK / yScale ) ) )

#ax.set_xlim( 1.0e2, 2.0e5 )

#ax.legend()

if useLogYScale:
    ax.set_yscale( 'log' )

ax.set_xlabel( r'$r/\mathrm{km}$' )
ax.set_ylabel( yLabel )

if len( ylim ) > 0 :
    ax.set_ylim( ylim )

xRef = [ 5.0000e+4, 2.50000e+4, 1.250000e+4, \
         6.2500e+3, 3.12500e+3, 1.562500e+3, \
         7.8125e+2, 3.90625e+2, 1.953125e+2 ]
for xx in xRef:
    ax.axvline( xx, color = 'r', alpha = 0.5 )

plt.show()

#figName = 'fig.YahilCollapse_native.png'
#plt.savefig( figName, dpi = 300 )
#print( '\n  Saved {:}'.format( figName ) )

import os
os.system( 'rm -rf __pycache__' )
