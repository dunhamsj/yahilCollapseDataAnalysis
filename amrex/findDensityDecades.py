#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

import sys
sys.path.append( '../' )

from UtilitiesModule import GetFileArray
from MakeDataFile import MakeDataFile, ReadHeader
from setGlobalVariables import *

class DensityDecadesAMReX:

  def __init__( self, plotFileDirectory, plotFileBaseName, dataFileDirectory ):

    self.plotFileArray = GetFileArray( plotFileDirectory, plotFileBaseName )
    self.dataFileDirectory = dataFileDirectory

    self.SSi = 0
    self.SSf = self.plotFileArray.shape[0] - 1
    self.nSS = self.plotFileArray.shape[0]

    MakeDataFile( 'PF_D', plotFileDirectory, dataFileDirectory, \
                  plotFileBaseName, 'spherical', \
                  SSi = self.SSi, SSf = self.SSf, nSS = self.nSS, \
                  UsePhysicalUnits = True, \
                  MaxLevel = -1, Verbose = True )

    if( not self.dataFileDirectory[-1] == '/' ): self.dataFileDirectory += '/'

    return

  def FindDensityDecades( self ):

    self.DensityDecades = np.logspace( 10, 15, 6 )
    self.FoundDecade    \
      = np.array( [ False, False, False, False, False, False ], bool )

    self.ind = np.empty( (self.DensityDecades.shape[0]), np.int64 )
    self.t   = np.empty( (self.DensityDecades.shape[0]), np.int64 )

    Density    = 0.0
    MaxDensity = 0.0
    tMax       = 0.0
    indMax     = 0

    for t in range( self.nSS ):

      iSS = self.SSi \
              + np.int64( ( self.SSf - self.SSi  ) / ( self.nSS - 1 ) * t )

      if t % 100 == 0:
        print( '  {:}/{:}'.format( t, self.nSS ) )

      dataDir \
        = self.dataFileDirectory + self.plotFileArray[iSS][-8:] + '/'

      dataFile = dataDir + 'PF_D.dat'
      timeFile = dataDir + 'Time.dat'

      DataShape, DataUnits, MinVal, MaxVal = ReadHeader( dataFile )

      Data = np.loadtxt( dataFile ).reshape( np.int64( DataShape ) )
      Density = Data[0]

      Time = np.loadtxt( timeFile )

      if Density > MaxDensity:

        MaxDensity = Density
        tMax       = Time
        indMax     = t

      for iDec in range( self.DensityDecades.shape[0] ):

        if Density > self.DensityDecades[iDec] \
             and not self.FoundDecade[iDec]:

          self.FoundDecade[iDec] = True
          self.ind        [iDec] = t
          self.t          [iDec] = Time

    string = ''
    string += str( self.ind ) + '\n'
    string += 'MaxDensity: {:.3e} g/cm^3\n'.format( MaxDensity )
    string += '      tMax: {:.3e} ms\n'.format( tMax )
    string += '    indMax: {:d}'.format( indMax )
    print(string)

    with open( 'densityDecades.txt', 'w' ) as f:
      f.write( string )

    self.MaxDensity = MaxDensity
    self.tMax       = tMax
    self.indMax     = indMax

if __name__ == '__main__':

  ID = 'YahilCollapse_XCFC'

  rootDir = HOME + 'Work/Codes/thornado/SandBox/AMReX/Applications/'
  plotFileDirectory = rootDir + 'YahilCollapse_XCFC/'

  plotFileBaseName = ID + '.plt'

  dataFileDirectory = '.{:s}_movieData'.format( ID )

  DD = DensityDecadesAMReX \
         ( plotFileDirectory, plotFileBaseName, dataFileDirectory )
  DD.FindDensityDecades()

  import os
  os.system( 'rm -rf __pycache__' )
