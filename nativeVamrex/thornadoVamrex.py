#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
plt.style.use( '../publication.sty' )

import sys
sys.path.append( '../' )

from MakeDataFile import MakeDataFile, ReadHeader
from setHomeDirectory import *

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
Field   = 'PF_D'
FieldT  = 'PF_D'
yScale  = 1.0e0
yScaleT = 1.0e0
dataT = np.loadtxt( '{:}_native_{:}.dat'.format( rootName, FieldT ) )
UseLogScale_Y = False
if Field == 'PF_V1':
  yMin = -0.15
  yMax = 0.01
elif Field == 'PF_D':
  UseLogScale_Y = True
  yMin = 1.0
  yMax = 1.0e15
else:
  yMin = dataT[1:,1:].min()
  yMax = dataT[1:,1:].max()

# Plot data in log10-scale?
UseLogScale_Y = True
UseLogScale_X = True

# Unit system of the data
UsePhysicalUnits = True

# Coordinate system (currently supports 'cartesian' and 'spherical' )
CoordinateSystem = 'spherical'

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

UseCustomLimits = False
vmin = 0.0
vmax = 2.0

MovieRunTime = 10.0 # seconds

#### ====== End of User Input =======

DataDirectory = '.{:s}_movieData'.format( ID )
MovieName     = 'mov.{:s}_{:s}.mp4'.format( ID, Field )

# Append "/" if not present
if not plotfileDirectory[-1] == '/': plotfileDirectory += '/'
if not DataDirectory    [-1] == '/': DataDirectory     += '/'

TimeUnits = ''
X1Units   = ''
if UsePhysicalUnits:
    TimeUnits = 'ms'
    X1Units   = 'km'

plotfileArray \
  = MakeDataFile( Field, plotfileDirectory, DataDirectory, \
                  plotfileBaseName, CoordinateSystem, \
                  SSi = SSi, SSf = SSf, nSS = nSS, \
                  forceChoiceD = False, owD = False, \
                  forceChoiceF = False, owF = False, \
                  UsePhysicalUnits = UsePhysicalUnits, \
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

vmin = min( dataA.min(), dataT.min() )
vmax = max( dataA.max(), dataT.max() )

dX1 = X1_C[1] - X1_C[0]
xL  = X1_C[0 ] - 0.5 * dX1
xH  = X1_C[-1] + 0.5 * dX1

fig, ax = plt.subplots( 1, 1 )
ax.set_title( r'$\texttt{{{:}}}$'.format( ID ), fontsize = 15 )

time_textA = ax.text( 0.1, 0.9, '', transform = ax.transAxes, fontsize = 13 )
time_textT = ax.text( 0.1, 0.8, '', transform = ax.transAxes, fontsize = 13 )

ax.set_xlabel( r'$x/\mathrm{km}$', fontsize = 15 )

ax.set_xlim( xL + 0.25 * dX1, xH + 1.0e5 )
ax.set_ylim( yMin, yMax )

ax.set_xscale( 'log' )
if UseLogScale_Y: ax.set_yscale( 'log' )

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

    print('    {:}/{:}'.format( t, nSS ) )

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

fps = max( 1, nSS / MovieRunTime )

print( '\n  Making movie' )
print( '  ------------' )
anim.save( MovieName, fps = fps, dpi = 300 )

import os
os.system( 'rm -rf __pycache__ ' )
