#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
plt.style.use( '../publication.sty' )

import sys
sys.path.append( '../' )

from UtilitiesModule import GetData, GetFileArray
from MakeDataFile import MakeDataFile, ReadHeader
from setGlobalVariables import *

class MultiPanel:

  def __init__( self, plotFileDirectory, \
                      snapShots = np.array( [ -1 ], np.int64 ) ):

    self.plotFileDirectory = plotFileDirectory
    self.plotFileBaseName = 'YahilCollapse_XCFC.plt'
    self.snapShots = snapShots

    return

  def GetData( self, Field ):

    plotFileArray \
      = GetFileArray( self.plotFileDirectory, self.plotFileBaseName )

    Data = np.empty( self.snapShots.shape[0], object )
    X1   = np.empty( self.snapShots.shape[0], object )
    dX1  = np.empty( self.snapShots.shape[0], object )
    for iSS in range( self.snapShots.shape[0] ):

      d, DataUnit, r, X2, X3, dr, dX2, dX3, xl, xh, nX, Time \
        = GetData( self.plotFileDirectory, self.plotFileBaseName, Field, \
                   'spherical', True, \
                   argv = ['x',str( plotFileArray[self.snapShots[iSS]] )], \
                   MaxLevel = -1, \
                   ReturnTime = True, ReturnMesh = True, Verbose = False )
      Data[iSS] = d [:,0,0]
      X1  [iSS] = r [:,0,0]
      dX1 [iSS] = dr[:,0,0]

    return Data, DataUnit, X1, dX1, xl, xh, Time

if __name__ == '__main__':

  saveFig = False
  saveFigAs = 'fig.AdibaticCollapse_XCFC_MultiPanel.png'

  plotFileDirectory \
    = HOME + 'Work/Codes/thornado/SandBox/AMReX/Applications/YahilCollapse_XCFC/'

#[                546                1074                1421
#                1454                1464 4624633867356078080]
#MaxDensity: 9.991e+14 g/cm^3
#      tMax: 1.467e+02 ms
#    indMax: 1467

  snapShots = np.array( [ 546, 1074, 1421, 1454, 1464, 1467 ], np.int64 )
  nSS = snapShots.shape[0]

  MP \
    = MultiPanel \
        ( plotFileDirectory, snapShots = snapShots )

  PF_D , unit_PF_D , X1, dX1, xl, xh, Time = MP.GetData( 'PF_D'      )
  GF_b1, unit_GF_b1, X1, dX1, xl, xh, Time = MP.GetData( 'GF_Beta_1' )
  PF_V1, unit_PF_V1, X1, dX1, xl, xh, Time = MP.GetData( 'PF_V1'     )
  AF_P , unit_AF_P , X1, dX1, xl, xh, Time = MP.GetData( 'AF_P'      )
  GF_CF, unit_GF_CF, X1, dX1, xl, xh, Time = MP.GetData( 'GF_Psi'    )
  GF_Al, unit_GF_Al, X1, dX1, xl, xh, Time = MP.GetData( 'GF_Alpha'  )

  color = ['#377eb8', '#ff7f00', '#4daf4a', \
           '#f781bf', '#a65628', '#984ea3', \
           '#999999', '#e41a1c', '#dede00']

  fig, axs = plt.subplots( 2, 2, figsize = ( 12,9 ) )

  for i in range( 2 ):
    for j in range( 2 ):
      axs[i,j].set_xlim( xL, xH )
      axs[i,j].set_xscale( 'log')
      axs[i,j].grid()
      axs[i,j].set_xlabel( xLabel )

  axs[0,0].set_ylabel( r'$\rho\,\left[\mathrm{g\,cm}^{-3}\right]$' )
  axs[1,0].set_ylabel( r'$p\,\left[\mathrm{erg\,cm}^{-3}\right]$' )

  for iSS in range( nSS ):

    axs[0,0].semilogy( X1[iSS], PF_D [iSS], '-', c = color[iSS], \
                       label = 'Time = {:.3e} ms'.format( Time ) )
    axs[1,0].semilogy( X1[iSS], AF_P [iSS], '-', c = color[iSS] )

    if ( iSS == 0 ) :
      axs[0,1].plot( X1[iSS], GF_b1[iSS] / 2.99792458e5, '-' , c = color[iSS], label = r'$\beta^{1}/c$' )
      axs[0,1].plot( X1[iSS], PF_V1[iSS] / 2.99792458e5, '--', c = color[iSS], label = r'$v^{1}/c$' )
      axs[1,1].plot( X1[iSS], GF_CF[iSS], '-' , c = color[iSS], label = r'$\psi$' )
      axs[1,1].plot( X1[iSS], GF_Al[iSS], '--', c = color[iSS], label = r'$\alpha$' )
    else:
      axs[0,1].plot( X1[iSS], GF_b1[iSS] / 2.99792458e5, '-' , c = color[iSS] )
      axs[0,1].plot( X1[iSS], PF_V1[iSS] / 2.99792458e5, '--', c = color[iSS] )
      axs[1,1].plot( X1[iSS], GF_CF[iSS], '-' , c = color[iSS] )
      axs[1,1].plot( X1[iSS], GF_Al[iSS], '--', c = color[iSS] )

  axs[1,1].legend()
  axs[0,1].legend()

  fig.suptitle( r'$\texttt{YahilCollapse_XCFC}$', fontsize = 15 )

  if saveFig:
    plt.savefig( saveFigAs, dpi = 300, bbox_inches = 'tight' )
  else:
    plt.show()

  import os
  os.system( 'rm -rf __pycache__ ' )
