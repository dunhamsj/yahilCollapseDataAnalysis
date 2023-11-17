#!/usr/bin/env python3

import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
plt.style.use( '../publication.sty' )

import sys
sys.path.append( '../' )

from UtilitiesModuleHDF import ReadFieldsHDF
from setGlobalVariables import *

############################ User Input ############################

THORNADO_DIR = HOME + 'Work/Codes/thornado/'

RootPath = THORNADO_DIR + 'SandBox/YahilCollapse_XCFC/'
suffix = 'Output/'

Problem = 'YahilCollapse'

RunTime = 10.0 # seconds

Fields           = [ 'PF_V1' ]
Dimension        = 'X1'
SnapshotRange    = [0,625]
plotEvery        = 1
WriteFile        = True

if Fields[0] == 'PF_V1':

  UseCustomLimits_Y = False
  UseLogScale_Y     = False
  yScale = [ 2.99792458e5 ]
  labels = [ r'$v/c$' ]

elif Fields[0] == 'PF_D':

  UseCustomLimits_Y = False#; yMin = 1.0; yMax = 1.0e15
  UseLogScale_Y     = True
  yScale = [ 1.0e0 ]
  labels = [ r'$\rho\ \left[\mathrm{g\,cm}^{-3}\right]$' ]

elif Fields[0] == 'PolytropicConstant':

  UseCustomLimits_Y = True; yMin = 1.0 - 1.0e-4; yMax = 1.0 + 1.0e-4
  UseLogScale_Y     = False
  yScale = [ 6.0e27 / 7.0e9**1.30 ]
  labels = [ r'$K/K_{\mathrm{exact}}$' ]

elif Fields[0] == 'AF_P':

  UseCustomLimits_Y = False
  UseLogScale_Y     = True
  yScale = [ 1.0 ]
  labels = [ r'$p\ \left[\mathrm{erg\,cm}^{-3}\right]$' ]

else:

  UseCustomLimits_Y = False
  UseLogScale_Y     = False
  yScale = [ 1.0e0 ]
  labels = Fields

figTitle = 'Yahil Collapse, 8192 elements, plF, slF'

############################

print( '\n  makeMovie_native.py' )
print( '  -------------------' )
print( '    plotfileDirectory: ', RootPath + suffix )
print( '               Fields: ', Fields )
print( '            WriteFile: ', WriteFile )
print()

nFields = len( Fields )

nSS = SnapshotRange[1] - SnapshotRange[0] + 1

Snapshots \
  = np.linspace( SnapshotRange[0], SnapshotRange[1], nSS, \
                 dtype = np.int64 )
Snapshots = np.copy( Snapshots[::plotEvery] )
nSS = Snapshots.shape[0]

PathToData = RootPath + suffix + Problem

tmp = Fields[0]
for iFd in range( 1, nFields ):
    tmp += '_' + Fields[iFd]
SaveFileAs = 'mov.{:}_XCFC_{:}_native.mp4'.format( Problem, tmp )

Names = ReadFieldsHDF \
          ( PathToData, Snapshots, \
            CoordinateSystem = 'SPHERICAL', \
            UsePhysicalUnits = True, \
            UseGeometryFields = True )
TimeUnit = Names['Time'][0]
Time     = Names['Time'][1]

XC   = np.array( Names['X1_C'][1] )
dX   = np.diff( XC )
xlim = [ xL, xH ]

YN = np.empty( nFields, object )

for iFd in range( nFields ):
  YN[iFd] = Names[Fields[iFd]][1][:,0,0,:]

if not UseCustomLimits_Y:

  yMin = +np.inf
  yMax = -np.inf

  for iFd in range( nFields ):
    yMin = min( yMin, YN[iFd].min() ) / yScale[iFd]
    yMax = max( yMax, YN[iFd].max() ) / yScale[iFd]

ylim = [ yMin, yMax ]

################ Plotting information

# Animation program adapted from
# https://jakevdp.github.io/blog/2012/08/18/matplotlib-animation-tutorial/
fig, ax = plt.subplots( 1, 1 )
ax.set_title( '{:}'.format( figTitle ) )

time_text = plt.text( 0.1, 0.8, '', transform = ax.transAxes )

ax.set_xlim( xlim )
ax.set_xscale( 'log' )
ax.set_xlabel( xLabel )

ax.xaxis.set_tick_params \
  ( which = 'both', top = True, left = False, bottom = True, right = False  )
ax.yaxis.set_tick_params \
  ( which = 'both', top = False, left = True, bottom = False, right = True  )

ax.set_ylim( ylim )
if UseLogScale_Y:
  if( yMin < 0.0 ):
    ax.set_yscale( 'symlog' )
  else:
    ax.set_yscale( 'log' )

# colorblind-friendly palette: https://gist.github.com/thriveth/8560036
color = ['#377eb8', '#ff7f00', '#4daf4a', \
         '#f781bf', '#a65628', '#984ea3', \
         '#999999', '#e41a1c', '#dede00']

linesN = np.empty( nFields, object )
linesK = np.empty( nFields, object )
for iFd in range( nFields ):
  linesN[iFd], = ax.plot( [], [], '.', color = color[iFd], \
                          label = labels[iFd] + ' (N)' )
  linesK[iFd], = ax.plot( [], [], '.', color = color[iFd+1], \
                          label = labels[iFd] + ' (K)' )

ax.grid()
ax.legend()

XN = np.array( Names['X1'][1] )

nNodes = XN.shape[0] // XC.shape[0]
nX1    = XN.shape[0] // nNodes

if   nNodes == 1:
  wq = np.array( [ 1.0 ], np.float64 )
elif   nNodes == 2:
  wq = np.array( [ 0.5, 0.5 ], np.float64 )
elif nNodes == 3:
  wq = np.array( [ 5.0, 8.0, 5.0 ], np.float64 ) / 18.0
else:
  exit( 'Not available for nNodes = {:}'.format( nNodes ) )

def computeCellAverage( Names, nNodes, iFd, nSS ):

  SqrtGm = Names['GF_SqrtGm'][1][:,0,0,:]

  uK = np.empty( (nSS,nX1), np.float64 )

  for iSS in range( nSS ):

    print( '\r  Computing cell average: {:}/{:}'.format( iSS+1, nSS ), \
           end = '\r' )

    for iX1 in range( nX1 ):

      iLo = nNodes * iX1
      iHi = iLo + nNodes

      vK = np.sum( wq * SqrtGm[iSS,iLo:iHi] )

      if Fields[iFd] == 'PolytropicConstant':

        p   = np.sum( wq * Names['AF_P' ][1][iSS,0,0,iLo:iHi] \
                             * SqrtGm[iSS,iLo:iHi] ) / vK
        rho = np.sum( wq * Names['PF_D' ][1][iSS,0,0,iLo:iHi] \
                             * SqrtGm[iSS,iLo:iHi] ) / vK
        Gmm = np.sum( wq * Names['AF_Gm'][1][iSS,0,0,iLo:iHi] \
                             * SqrtGm[iSS,iLo:iHi] ) / vK

        uK[iSS,iX1] = p / rho**Gmm

      else:

        uK[iSS,iX1] \
          = np.sum( wq * YN[iFd][iSS,iLo:iHi] * SqrtGm[iSS,iLo:iHi] ) / vK

  return uK

print()
YK = np.empty( nFields, object )
for iFd in range( nFields ):
  YK[iFd] = computeCellAverage( Names, nNodes, iFd, nSS )
  print()

if WriteFile:

  header  = 'data[0,1:]  = XC [km]\n'
  header += 'data[1:,0]  = Time [ms]\n'
  header += 'data[1:,1:] = uK'

  YKK = np.empty( (nSS+1,nX1+1), np.float64 )
  YKK[0,0 ]  = np.nan
  YKK[0,1:]  = XC
  YKK[1:,0]  = Time

  for iFd in range( nFields ):

    YKK[1:,1:] = YK[iFd]

    fileName = '{:}_native_{:}.dat'.format( Problem, Fields[iFd] )
    np.savetxt( fileName, YKK, header = header )

    print( '  Wrote file ', fileName )

  #os.system( 'rm -rf __pycache__' )
  #exit()

# Intialize each new frame
def InitializeFrame():

  ret = []
  for iFd in range( nFields ):

      linesN[iFd].set_data([],[])
      ret.append( linesN[iFd] )

      linesK[iFd].set_data([],[])
      ret.append( linesK[iFd] )

  time_text.set_text('')
  ret.append( time_text )
  ret = ( ret )

  return ret

# Animation function
def UpdateFrame(t):

  print( '\r          Updating frame: {:}/{:}'.format( t+1, nSS ), end = '\r' )

  ret = []

  for iFd in range( nFields ):

    linesN[iFd].set_data( XN, YN[iFd][t] / yScale[iFd] )
    ret.append( linesN[iFd] )

    linesK[iFd].set_data( XC, YK[iFd][t] / yScale[iFd] )
    ret.append( linesK[iFd] )

  time_text.set_text( r'$t={:.3e}\ \mathrm{{ms}}$'.format( Time[t] ) )
  ret.append( time_text )
  ret = ( ret )

  return ret

print()
# Call the animator
anim = animation.FuncAnimation \
         ( fig, \
           UpdateFrame, \
           init_func = InitializeFrame, \
           frames    = nSS, \
           interval  = 1000, \
           blit      = True )

anim.save( SaveFileAs, fps = max( 1, int( nSS / RunTime ) ), dpi = 300 )
print()

os.system( 'rm -f *.pyc' )
os.system( 'rm -rf __pycache__' )
