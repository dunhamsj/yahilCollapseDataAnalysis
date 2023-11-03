#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
plt.style.use( '../publication.sty' )

import sys
sys.path.append( '../' )

from MakeDataFile import MakeDataFile, ReadHeader
from setGlobalVariables import *

#### ========== User Input ==========

# ID to be used for naming purposes
ID = 'YahilCollapse_XCFC'

# Directory containing AMReX plotFiles
#plotFileDirectory = 'thornado/SandBox/AMReX/'
plotFileDirectory \
  = HOME + 'Work/Codes/thornado/\
SandBox/AMReX/Applications/YahilCollapse_XCFC/'

# plotFile base name (e.g., Advection1D.plt######## -> Advection1D.plt )
plotFileBaseName = ID + '.plt'

# Field to plot
Field   = 'PF_D'
yScale  = 1.0e0
yLabel  = r'$\rho\,\left[\mathrm{g\,cm}^{-3}\right]$'
UseLogScale_Y = True

# Only use every <plotEvery> plotFile
plotEvery = 10

# First and last snapshots and number of snapshots to include in movie
SSi = -1 # -1 -> SSi = 0
SSf = -1 # -1 -> plotFileArray.shape[0] - 1
nSS = -1 # -1 -> plotFileArray.shape[0]

# Max level of refinement to include
MaxLevel = -1

# Include initial conditions in movie?
ShowIC = True

PlotMesh = True

Verbose = True

UseCustomLimits = True
yMin = 1.0
yMax = 1.0e15

MovieRunTime = 10.0 # seconds

#### ====== End of User Input =======

DataDirectory = '.{:s}_movieData'.format( ID )
MovieName     = 'mov.{:s}_{:s}.mp4'.format( ID, Field )

# Append "/" if not present
if not plotFileDirectory[-1] == '/': plotFileDirectory += '/'
if not DataDirectory    [-1] == '/': DataDirectory     += '/'

plotFileArray \
  = MakeDataFile( Field, plotFileDirectory, DataDirectory, \
                  plotFileBaseName, 'spherical', \
                  SSi = SSi, SSf = SSf, nSS = nSS, \
                  forceChoiceD = False, owD = False, \
                  forceChoiceF = False, owF = False, \
                  UsePhysicalUnits = True, \
                  MaxLevel = MaxLevel, Verbose = Verbose )
plotFileArray = np.copy( plotFileArray[::plotEvery] )
if nSS < 0: nSS = plotFileArray.shape[0]

time = np.empty( nSS, np.float64 )
X1_C = np.empty( nSS, object )
data = np.empty( nSS, object )

vmin = +np.inf
vmax = -np.inf

for t in range( nSS ):

  FileDirectory = DataDirectory + plotFileArray[t] + '/'

  TimeFile   = FileDirectory + '{:}.dat'.format( 'Time' )
  X1DataFile = FileDirectory + '{:}.dat'.format( 'X1' )
  DataFile   = FileDirectory + '{:}.dat'.format( Field )

  time[t] = np.loadtxt( TimeFile )
  X1_C[t] = np.loadtxt( X1DataFile ).flatten()
  data[t] = np.loadtxt( DataFile ).flatten()

  DataShape, DataUnits, MinVal, MaxVal = ReadHeader( DataFile )

  vmin = min( vmin, MinVal )
  vmax = max( vmax, MaxVal )

fig, ax = plt.subplots( 1, 1 )
time_text = ax.text( 0.1, 0.9, '', transform = ax.transAxes, fontsize = 13 )

line, = ax.plot( [],[], 'k-', label = r'$u$' )

def InitializeFrame():

  line.set_data([],[])
  time_text.set_text('')

  ret = ( line, time_text )

  return ret

def UpdateFrame( t ):

  print('    {:}/{:}'.format( t+1, nSS ) )

  time_text.set_text( r'$t={:.16e}\ \mathrm{{ms}}$' \
                      .format( time[t] ) )

  line.set_data( X1_C[t], data[t] / yScale )

  ret = ( line, time_text )

  return ret

anim = animation.FuncAnimation( fig, UpdateFrame, \
                                init_func = InitializeFrame, \
                                frames = nSS, \
                                blit = True )
if not UseCustomLimits:
  yMin = vmin
  yMax = vmax

ax.set_title( r'$\texttt{{{:}}}$'.format( ID ), fontsize = 15 )

ax.set_xscale( 'log' )
if UseLogScale_Y: ax.set_yscale( 'log' )

ax.set_xlabel( xLabel )
ax.set_ylabel( yLabel )

ax.set_xlim( xL, xH )
ax.set_ylim( yMin, yMax )

fps = max( 1, nSS / MovieRunTime )

print( '\n  Making movie' )
print( '  ------------' )
anim.save( MovieName, fps = fps, dpi = 300 )

import os
os.system( 'rm -rf __pycache__ ' )
