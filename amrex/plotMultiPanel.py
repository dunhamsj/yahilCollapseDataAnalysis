#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
plt.style.use( '../publication.sty' )

import sys
sys.path.append( '../' )

from UtilitiesModule import GetData, GetFileArray
from MakeDataFile import MakeDataFile, ReadHeader

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

    for iSS in range( self.snapShots.shape[0] ):

      Data, DataUnit, X1, X2, X3, dX1, dX2, dX3, xL, xU, nX, Time \
        = GetData( self.plotFileDirectory, self.plotFileBaseName, Field, \
                   'spherical', True, \
                   argv = ['x',str( plotFileArray[self.snapShots[iSS]] )], \
                   MaxLevel = -1, \
                   ReturnTime = True, ReturnMesh = True, Verbose = False )

    return Data[:,0,0], DataUnit, X1[:,0,0], dX1[:,0,0], xL, xU, Time

if __name__ == '__main__':

  saveFig = False
  saveFigAs = 'fig.AdibaticCollapse_XCFC_MultiPanel.png'

  HOME = '/home/kkadoogan/'
  plotFileDirectory \
    = HOME + 'Work/Codes/thornado/SandBox/AMReX/Applications/YahilCollapse_XCFC/'

  MP \
    = MultiPanel \
        ( plotFileDirectory, snapShots = np.array( [ 0 ], np.int64 ) )

  PF_D , unit_PF_D , X1, dX1, xL, xU, Time = MP.GetData( 'PF_D'      )
  GF_b1, unit_GF_b1, X1, dX1, xL, xU, Time = MP.GetData( 'GF_Beta_1' )
  PF_V1, unit_PF_V1, X1, dX1, xL, xU, Time = MP.GetData( 'PF_V1'     )
  AF_P , unit_AF_P , X1, dX1, xL, xU, Time = MP.GetData( 'AF_P'      )
  GF_CF, unit_GF_CF, X1, dX1, xL, xU, Time = MP.GetData( 'GF_Psi'    )
  GF_Al, unit_GF_Al, X1, dX1, xL, xU, Time = MP.GetData( 'GF_Alpha'  )

  fig, axs = plt.subplots( 2, 2, figsize = ( 12,9 ) )

  for i in range( 2 ):
    for j in range( 2 ):
      axs[i,j].set_xlim( 1.0, 2.0e5 )
      axs[i,j].set_xscale( 'log')
      axs[i,j].grid()
      axs[i,j].set_xlabel( 'Radial Coordinate [km]' )

  axs[0,0].semilogy( X1, PF_D, 'b-', label = 'Time = {:.3e} ms'.format( Time ) )
  axs[0,0].set_ylabel( 'PF_D ' )
  axs[0,1].plot    ( X1, GF_b1 / 2.99792458e5, 'b-', label = r'$\beta^{1}/c$' )
  axs[0,1].plot    ( X1, PF_V1 / 2.99792458e5, 'b--', label = r'$v^{1}/c$' )
  axs[1,0].semilogy( X1, AF_P, 'b-' )
  axs[1,0].set_ylabel( 'AF_P ' )
  axs[1,1].plot    ( X1, GF_CF, 'b-', label = r'$\psi$' )
  axs[1,1].plot    ( X1, GF_Al, 'b--', label = r'$\alpha$' )

  MP \
    = MultiPanel \
        ( plotFileDirectory, snapShots = np.array( [ -1 ], np.int64 ) )

  PF_D , unit_PF_D , X1, dX1, xL, xU, Time = MP.GetData( 'PF_D'      )
  GF_b1, unit_GF_b1, X1, dX1, xL, xU, Time = MP.GetData( 'GF_Beta_1' )
  PF_V1, unit_PF_V1, X1, dX1, xL, xU, Time = MP.GetData( 'PF_V1'     )
  AF_P , unit_AF_P , X1, dX1, xL, xU, Time = MP.GetData( 'AF_P'      )
  GF_CF, unit_GF_CF, X1, dX1, xL, xU, Time = MP.GetData( 'GF_Psi'    )
  GF_Al, unit_GF_Al, X1, dX1, xL, xU, Time = MP.GetData( 'GF_Alpha'  )

  axs[0,0].semilogy( X1, PF_D, 'r-', label = 'Time = {:.3e} ms'.format( Time ) )
  axs[0,0].legend()
  axs[0,1].plot    ( X1, GF_b1 / 2.99792458e5, 'r-', label = r'$\beta^{1}/c$' )
  axs[0,1].plot    ( X1, PF_V1 / 2.99792458e5, 'r--', label = r'$v^{1}/c$' )
  axs[0,1].legend()
  axs[1,0].semilogy( X1, AF_P, 'r-' )
  axs[1,1].plot    ( X1, GF_CF, 'r-', label = r'$\psi$' )
  axs[1,1].plot    ( X1, GF_Al, 'r--', label = r'$\alpha$' )
  axs[1,1].legend()
  fig.suptitle( r'$\texttt{YahilCollapse_XCFC}$', fontsize = 15 )

  if saveFig:
    plt.savefig( saveFigAs, dpi = 300, bbox_inches = 'tight' )
  else:
    plt.show()

  import os
  os.system( 'rm -rf __pycache__ ' )
