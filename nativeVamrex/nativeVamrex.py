#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
plt.style.use( '../publication.sty' )

import sys
sys.path.append( '../' )

from MakeDataFile import MakeDataFile, ReadHeader
from setGlobalVariables import *

"""

Creates a directory with structure as laid out
in MakeDataFile.py and makes a movie from it

Usage:
  $ python3 makeMovie1D.py

"""

#### ========== User Input ==========

rootName = 'YahilCollapse'

# ID to be used for naming purposes
ID = '{:}_XCFC'.format( rootName )

# Directory containing AMReX plotfiles
#plotfileDirectory = 'thornado/SandBox/AMReX/'
plotfileDirectory \
  = HOME + 'Work/Codes/thornado/\
SandBox/AMReX/Applications/YahilCollapse_XCFC/'

# plotfile base name (e.g., Advection1D.plt######## -> Advection1D.plt )
plotfileBaseName = ID + '.plt'

# Field to plot
Field  = 'PF_V1'
FieldT = 'PF_V1'
dataT = np.loadtxt( '{:}_native_{:}.dat'.format( rootName, FieldT ) )
X1_C = np.copy( dataT[0,1:] )
X2_C = [np.pi/2.0]
X3_C = [np.pi]
LeafElementLocations = []
for i in range( X1_C.shape[0] ):
    LeafElementLocations.append \
      ( np.array( [ X1_C[i], X2_C[0], X3_C[0] ], np.float64 ) )

figTitle = 'Yahil Collapse, Multi-Level Mesh\n{:}'.format( Field )

if   Field == 'PF_V1':

  UseCustomLimits = False
  UseLogScale_Y   = False
  yScale  = 2.99792458e5
  yLabel0 = r'$v/c$'
  yLabel1 = r'$\left|v_{\mathrm{a}}-v_{\mathrm{t}}\right|/$' \
              + r'\frac{1}{2}\left|v_{\mathrm{a}}+v_{\mathrm{t}}\right|$'

elif Field == 'PF_D':

  UseLogScale_Y   = True
  UseCustomLimits = True; yMin = 1.0; yMax = 1.0e15
  yScale  = 1.0e0
  yLabel0 = r'$\rho\,\left[\mathrm{g\,cm}^{-3}\right]$'
  yLabel1 = r'$\left|\rho_{\mathrm{a}}-\rho_{\mathrm{t}}\right|/$' \
              + r'\frac{1}{2}\left|\rho_{\mathrm{a}}+\rho_{\mathrm{t}}\right|$'

elif Field == 'PolytropicConstant':

  UseLogScale_Y   = False
  UseCustomLimits = False
  yScale  = 6.0e27 / 7.0e9**1.30
  yLabel0 = r'$K/K_{\mathrm{exact}}$'
  yLabel1 = r'$\left|K_{\mathrm{a}}-K_{\mathrm{t}}\right|/$' \
              + r'\frac{1}{2}\left|K_{\mathrm{a}}+K_{\mathrm{t}}\right|$'

elif Field == 'AF_P':

  UseLogScale_Y   = True
  UseCustomLimits = False
  yScale  = 1.0e0
  yLabel0 = r'$p\ \left[\mathrm{erg\,cm}^{-3}\right]$'
  yLabel1 = r'$\left|p_{\mathrm{a}}-p_{\mathrm{t}}\right|/$' \
              + r'\frac{1}{2}\left|p_{\mathrm{a}}+p_{\mathrm{t}}\right|$'

else:

  UseLogScale_Y   = False
  UseCustomLimits = False
  yScale  = 1.0e0
  yLabel0 = ''
  yLabel1 = ''

# Only use every <plotEvery> plotfile
plotEvery = 1

# First and last snapshots and number of snapshots to include in movie
SSi = -1#-1 -> SSi = 0
SSf = -1#-1 -> plotfileArray.shape[0] - 1
nSS = -1#-1 -> plotfileArray.shape[0]

# Max level of refinement to include
MaxLevel = -1

# Include initial conditions in movie?
ShowIC = True

PlotMesh = True

Verbose = True

vmin = +np.inf
vmax = -np.inf

MovieRunTime = 10.0 # seconds

#### ====== End of User Input =======

DataDirectory = '.{:s}_movieData'.format( ID )
MovieName     = 'mov.{:s}_{:s}_nativeVamrex.mp4'.format( ID, Field )

# Append "/" if not present
if not plotfileDirectory[-1] == '/': plotfileDirectory += '/'
if not DataDirectory    [-1] == '/': DataDirectory     += '/'

plotfileArray \
  = MakeDataFile( Field, plotfileDirectory, DataDirectory, \
                  plotfileBaseName, 'spherical', \
                  SSi = SSi, SSf = SSf, nSS = nSS, \
                  forceChoiceD = False, owD = False, \
                  forceChoiceF = False, owF = False, \
                  UsePhysicalUnits = True, \
                  MaxLevel = MaxLevel, Verbose = Verbose, \
                  LEL = LeafElementLocations )
if nSS < 0:
  plotfileArray = np.copy( plotfileArray[::plotEvery] )
else:
  plotfileArray = np.copy( plotfileArray[SSi:SSf+1:plotEvery] )

nSS = min( plotfileArray.shape[0], dataT[1:,0].shape[0] )

nX = np.shape( X1_C )[0]

timeT = np.copy( dataT[1:nSS+1:plotEvery,0 ] )
dataT = np.copy( dataT[1:nSS+1:plotEvery,1:] )
timeA = np.empty( timeT.shape, np.float64 )
dataA = np.empty( dataT.shape, np.float64 )
dataD = np.empty( dataT.shape, np.float64 )

print()
print( '  Computing relative difference array' )
for iSS in range( nSS ):

    FileDirectory = DataDirectory + plotfileArray[iSS] + '/'

    TimeFile = FileDirectory + '{:}.dat'.format( 'Time' )
    DataFile = FileDirectory + '{:}.dat'.format( Field )

    timeA[iSS] = np.loadtxt( TimeFile )
    dataA[iSS] = np.loadtxt( DataFile ).flatten()

    for iX1 in range( nX ):
        print( '\r  iSS: {:}/{:}'.format( iSS, nSS ), end = '\r' )
        dataD[iSS,iX1] \
          = max( np.abs( dataA[iSS,iX1] - dataT[iSS,iX1] ) \
                   / ( 0.5 * np.abs( dataA[iSS,iX1] + dataT[iSS,iX1] ) ), \
                 1.0e-17 )
print()
print(dataD.min(), dataD.max() )

fig, axs = plt.subplots( 2, 1, figsize = (10,6) )
fig.suptitle( '{:}'.format( figTitle ), fontsize = 15 )

time_textA \
  = axs[0].text( 0.1, 0.9, '', transform = axs[0].transAxes, fontsize = 13 )
time_textT \
  = axs[0].text( 0.1, 0.8, '', transform = axs[0].transAxes, fontsize = 13 )

lineA, = axs[0].plot( [],[], 'k-', lw = 2, label = r'$u_{\mathrm{amrex}}$' )
lineT, = axs[0].plot( [],[], 'r-', lw = 1, label = r'$u_{\mathrm{thrnd}}$' )
lineD, = axs[1].plot( [],[], 'k-', lw = 1 )


def InitializeFrame():

    lineA.set_data([],[])
    lineT.set_data([],[])
    lineD.set_data([],[])
    time_textA.set_text('')
    time_textT.set_text('')

    ret = ( lineA, lineT, lineD, time_textA, time_textT )

    return ret

def UpdateFrame( t ):

    print( '\r    Updating frame {:}/{:}'.format( t+1, nSS ), end = '\r' )

    time_textA.set_text( r'$t_{{\mathrm{{amrex}}}}={:.16e}\ \mathrm{{ms}}\ {:d}$' \
                         .format( timeA[t], t ) )
    time_textT.set_text( r'$t_{{\mathrm{{thrnd}}}}={:.16e}\ \mathrm{{ms}}$' \
                         .format( timeT[t] ) )

    lineA.set_data( X1_C, dataA[t] / yScale )
    lineT.set_data( X1_C, dataT[t] / yScale )
    lineD.set_data( X1_C, dataD[t] )

    ret = ( lineA, lineT, lineD, time_textA, time_textT )

    return ret

axs[0].legend( loc = 3, prop = {'size':12} )
xRef = [ 5.0000e+4, 2.50000e+4, 1.250000e+4, 6.2500e+3, 3.12500e+3, 1.562500e+3, 7.8125e+2, 3.90625e+2, 1.953125e+2 ]
for i in range( 2 ):
    for xx in xRef:
        axs[i].axvline( xx, color = 'b', alpha = 0.3 )

anim = animation.FuncAnimation( fig, UpdateFrame, \
                                init_func = InitializeFrame, \
                                frames = nSS, \
                                blit = True )

if not UseCustomLimits:
  yMin = min( dataA.min(), dataT.min() ) / yScale
  yMax = max( dataA.max(), dataT.max() ) / yScale

axs[1].set_xlabel( xLabel )
axs[0].set_ylabel( yLabel0 )
#axs[1].set_ylabel( yLabel1 )

axs[0].set_xlim( xL, xH )
axs[1].set_xlim( xL, xH )
axs[0].set_ylim( yMin, yMax )
axs[1].set_ylim( 1.0e-17, dataD.max() )

axs[0].set_xscale( 'log' )
axs[1].set_xscale( 'log' )
if UseLogScale_Y: axs[0].set_yscale( 'log' )
axs[1].set_yscale( 'log' )

fps = max( 1, nSS / MovieRunTime )

print( '\n  Making movie' )
print( '  ------------' )
anim.save( MovieName, fps = fps, dpi = 300 )
print()

import os
os.system( 'rm -rf __pycache__ ' )
