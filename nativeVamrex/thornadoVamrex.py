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

figTitle = 'Yahil Collapse (8192 elements)'

# plotfile base name (e.g., Advection1D.plt######## -> Advection1D.plt )
plotfileBaseName = ID + '.plt'

# Field to plot
Field   = 'PolytropicConstant'
FieldT  = 'PolytropicConstant'
yScale  = 6.0e27 / 7.0e9**1.30
yScaleT = 1.0e0
dataT = np.loadtxt( '{:}_native_{:}.dat'.format( rootName, FieldT ) )

yLabel = r'$K/K_{\mathrm{exact}}$'

UseLogScale_Y   = False
UseCustomLimits = False

if Field == 'PF_V1':
  UseCustomLimits = True
  yMin = -0.15
  yMax = 0.01
elif Field == 'PF_D':
  UseLogScale_Y   = True
  UseCustomLimits = True
  yMin = 1.0
  yMax = 1.0e15
elif Field == 'PolytropicConstant':
  UseCustomLimits = True
  yMin = 1.0 - 1.0e-1
  yMax = 1.0 + 1.0e-1
else:
  yMin = dataT[1:,1:].min()
  yMax = dataT[1:,1:].max()

# Only use every <plotEvery> plotfile
plotEvery = 1

# First and last snapshots and number of snapshots to include in movie
SSi = -1 # -1 -> SSi = 0
SSf = -1 # -1 -> plotfileArray.shape[0] - 1
nSS = -1 # -1 -> plotfileArray.shape[0]

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
                  MaxLevel = MaxLevel, Verbose = Verbose )
plotfileArray = np.copy( plotfileArray[::plotEvery] )

if nSS < 0: nSS = plotfileArray.shape[0]

X1_C = np.copy( dataT[0,1:] )
nX   = np.shape( X1_C )[0]

timeT = np.copy( dataT[1::plotEvery,0 ] )
dataT = np.copy( dataT[1::plotEvery,1:] )
timeA = np.empty( timeT.shape, np.float64 )
dataA = np.empty( dataT.shape, np.float64 )

for t in range( nSS ):

    FileDirectory = DataDirectory + plotfileArray[t] + '/'

    TimeFile = FileDirectory + '{:}.dat'.format( 'Time' )
    DataFile = FileDirectory + '{:}.dat'.format( Field )

    timeA[t] = np.loadtxt( TimeFile )
    dataA[t] = np.loadtxt( DataFile ).flatten()

fig, ax = plt.subplots( 1, 1 )
ax.set_title( '{:}'.format( figTitle ), fontsize = 15 )

time_textA = ax.text( 0.1, 0.9, '', transform = ax.transAxes, fontsize = 13 )
time_textT = ax.text( 0.1, 0.8, '', transform = ax.transAxes, fontsize = 13 )

lineA, = ax.plot( [],[], 'k-', lw = 2, label = r'$u_{\mathrm{amrex}}$' )
lineT, = ax.plot( [],[], 'r-', lw = 1, label = r'$u_{\mathrm{thrnd}}$' )

def InitializeFrame():

    lineA.set_data([],[])
    lineT.set_data([],[])
    time_textA.set_text('')
    time_textT.set_text('')

    ret = ( lineA, lineT, time_textA, time_textT )

    return ret

def UpdateFrame( t ):

    print( '\r    Updating frame {:}/{:}'.format( t+1, nSS ), end = '\r' )

    time_textA.set_text( r'$t_{{\mathrm{{amrex}}}}={:.16e}\ \mathrm{{ms}}$' \
                         .format( timeA[t] ) )
    time_textT.set_text( r'$t_{{\mathrm{{thrnd}}}}={:.16e}\ \mathrm{{ms}}$' \
                         .format( timeT[t] ) )

    lineA.set_data( X1_C, dataA[t] / yScale )
    lineT.set_data( X1_C, dataT[t] / yScaleT )

    ret = ( lineA, lineT, time_textA, time_textT )

    return ret

ax.legend( loc = 3, prop = {'size':12} )

anim = animation.FuncAnimation( fig, UpdateFrame, \
                                init_func = InitializeFrame, \
                                frames = nSS, \
                                blit = True )

if not UseCustomLimits:
  yMin = min( dataA.min(), dataT.min() )
  yMax = max( dataA.max(), dataT.max() )

ax.set_xlabel( r'$x/\mathrm{km}$', fontsize = 15 )
ax.set_ylabel( yLabel, fontsize = 15 )

ax.set_xlim( xL, xH )
ax.set_ylim( yMin, yMax )

ax.set_xscale( 'log' )
if UseLogScale_Y: ax.set_yscale( 'log' )

fps = max( 1, nSS / MovieRunTime )

print( '\n  Making movie' )
print( '  ------------' )
anim.save( MovieName, fps = fps, dpi = 300 )
print()

import os
os.system( 'rm -rf __pycache__ ' )
