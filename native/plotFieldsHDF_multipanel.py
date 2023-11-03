#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import subprocess
plt.style.use( '../publication.sty' )
plt.rcParams['legend.fontsize'] = 16
plt.rcParams['legend.loc'] = 'upper right'
plt.rcParams['axes.labelsize'] = 16
plt.rcParams['xtick.labelsize'] = 14
plt.rcParams['ytick.labelsize'] = 13
plt.rcParams['figure.titlesize'] = 16
plt.rcParams['figure.figsize'] = [12,9]

import sys
sys.path.append( '../' )

from UtilitiesModuleHDF import ReadFieldsHDF
from setGlobalVariables import *

THORNADO_DIR = HOME + 'Work/Codes/thornado/'

class PlotFieldsHDF:

    def __init__( self, DataDirectory, Snapshots ):

        self.DataDirectory = DataDirectory
        self.Snapshots     = np.array( Snapshots )

        if not self.DataDirectory[-1] == '/': self.DataDirectory += '/'

        self.PathToData = self.DataDirectory + 'YahilCollapse'

        self.nFiles = self.Snapshots.shape[0]

        return


    def GetData( self ):

        self.names \
          = ReadFieldsHDF \
              ( self.PathToData, self.Snapshots, \
                CoordinateSystem = 'SPHERICAL', UsePhysicalUnits = True )

        self.r    = self.names['X1'][1]
        self.rC   = self.names['X1_C'][1]
        self.Time = self.names['Time'][1]
        self.nX   = self.rC.shape[0]
        self.nN   = np.int64( self.r.shape[0] / self.rC.shape[0] )

        self.alpha = np.array( [ 1.0 ], np.float64 )

        if self.nN == 2:

          self.WeightsX = np.array( [ 0.5, 0.5 ], np.float64 )

        if self.nFiles > 1:

            self.alpha = np.linspace( 0.2, 1.0, self.nFiles )

        return


    def ComputeCellAverage( self, iSS, Field ):

        SqrtGm = self.names['GF_Sg'][1][iSS,0,0,:]
        uN     = self.names[Field  ][1][iSS,0,0,:]

        uK = np.empty( (self.nX), np.float64 )

        for iX1 in range( self.nX ):

            iLo = self.nN*iX1
            iHi = self.nN*(iX1+1)

            uK[iX1] \
              = np.sum( self.WeightsX * uN[iLo:iHi] * SqrtGm[iLo:iHi] ) \
                  / np.sum( self.WeightsX * SqrtGm[iLo:iHi] )

        return uK


    def AddPlot( self, ax, iSS, Field, c = 'k-', label = '', lw = 2.0, f = '' ):

        ss = 1.0
        if Field == 'PF_V1': ss = 2.99792458e5
        uK = self.ComputeCellAverage( iSS, Field ) / ss

        scale = 1.0
        if not f == '': scale = (2.99792458e10)**2
        if   f == 'lapse': uK = 1.0 + uK / scale
        elif f == 'psi'  : uK = 1.0 - 0.5 * uK / scale

        ax.plot( self.rC, uK, '.', \
                 alpha = self.alpha[iSS], lw = lw, label = str( label ) )

        return

if __name__ == '__main__':

    SaveFig, suffix = False, 'XCFC'

    #Snapshots = [ 969, 1312, 1420, 1453, 1464, 1467 ]
    Snapshots = [ 595, 887, 1016, 1072, 1097 ]

    tb = 1.465e2

    #PlotVariables = 'Fluid'
    PlotVariables = 'CFA'

    DataDirectory = THORNADO_DIR + 'SandBox/YahilCollapse_XCFC/'
    DataDirectory += 'Output/'

    Plot = PlotFieldsHDF( DataDirectory, Snapshots )
    Plot.GetData()

    FigTitle = 'Yahil Collapse (XCFC)'

    if PlotVariables == 'Fluid':

        Fields = np.array( [ 'PF_D', 'PF_V1', 'AF_P' ], str )

        nRows = 3
        nCols = 1
        fig, axs = plt.subplots( nRows, nCols )
        fig.suptitle( FigTitle )

        for iSS in range( Plot.nFiles ):

            for i in range( Fields.shape[0] ):

                if i == 0:

                    Plot.AddPlot( axs[i], iSS, Field = Fields[i], \
                                  label = r'$t-t_{{b}}={:+.6e}$ ms'.format \
                                    ( Plot.Time[iSS] - tb ) )

                else:

                    Plot.AddPlot( axs[i], iSS, Field = Fields[i] )

        for iF in range( nRows ):

            axs[iF].set_xscale( 'log' )
            axs[iF].set_xlim( xL, xH )
            axs[iF].tick_params( axis = 'both', labelsize = 14 )
            if iF < nRows-1: axs[iF].xaxis.set_ticklabels([])
            axs[iF].grid()

        axs[0].set_yscale( 'log' )
        axs[2].set_yscale( 'log' )

        axs[0].set_ylabel( r'$\rho\,\left[\mathrm{g\ cm}^{-3}\right]$' )
        axs[1].set_ylabel( r'$v/c$' )
        axs[2].set_ylabel( r'$p\,\left[\mathrm{erg\ cm}^{-3}\right]$' )

        axs[0].legend( prop = {'size':10} )
        axs[-1].set_xlabel( 'Radial Coordinate [km]' )

        plt.subplots_adjust( hspace = 0.0 )

        if SaveFig:

            plt.savefig( 'fig.YahilCollapse_CFA_{:}_Fluid.png'.format \
                         ( suffix ), dpi = 300 )

        else:

            plt.show()

        plt.close()

    if PlotVariables == 'CFA':

        Fields = np.array( [ 'GF_al', 'GF_CF', 'GF_b1', 'GF_NP' ], str )

        nRows = 2
        nCols = 1
        fig, axs = plt.subplots( nRows, nCols )
        fig.suptitle( FigTitle )

        for iSS in range( Plot.nFiles ):

            if iSS == Plot.nFiles - 1:

#                Plot.AddPlot( axs[0], iSS, Field = Fields[3], c = 'b--', \
#                              label = r'$\psi_{N}$', lw = 2, f = 'psi' )

                Plot.AddPlot( axs[0], iSS, Field = Fields[1], c = 'k--', \
                              label = r'$\psi$' )

#                Plot.AddPlot( axs[0], iSS, Field = Fields[3], c = 'b-', \
#                              label = r'$\alpha_{N}$', lw = 2, f = 'lapse' )

                Plot.AddPlot( axs[0], iSS, Field = Fields[0], c = 'k-', \
                              label = r'$\alpha$' )

            else:

#                Plot.AddPlot( axs[0], iSS, Field = Fields[3], c = 'b--', \
#                              lw = 2, f = 'psi' )
                Plot.AddPlot( axs[0], iSS, Field = Fields[1], c = 'k--' )
#                Plot.AddPlot( axs[0], iSS, Field = Fields[3], c = 'b-', \
#                              lw = 2, f = 'lapse' )
                Plot.AddPlot( axs[0], iSS, Field = Fields[0], c = 'k-' )

            Plot.AddPlot( axs[1], iSS, Field = Fields[2], \
                          label = r'$t-t_{{b}}={:+.2f}$ ms'.format \
                            ( Plot.Time[iSS] - tb ) )

        for iF in range( nRows ):

            axs[iF].set_xscale( 'log' )
            axs[iF].set_xlim( xL, xH )

        axs[1].set_ylabel( Fields[2] + ' ' + Plot.names[Fields[2]][0] )
        axs[1].set_ylabel( r'$\beta^{1}\,\left[\mathrm{km\ s}^{-1}\right]$' )
        #axs[0].xaxis.set_visible( False )
        axs[0].xaxis.set_ticklabels( [] )
        axs[0].grid()
        axs[1].grid()

        axs[0].legend()
        axs[1].legend()
        axs[-1].set_xlabel( 'Radial Coordinate [km]' )

        plt.subplots_adjust( hspace = 0.0 )

        if SaveFig:

            plt.savefig( 'fig.YahilCollapse_CFA_{:}_CFA.png'.format \
                         ( suffix ), dpi = 300 )

        else:

            plt.show()

        plt.close()

    import os
    os.system( 'rm -rf __pycache__' )
