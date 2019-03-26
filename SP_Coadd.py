#!/usr/bin/env python

""" SP_Coadd - Coadd fits files
    v1.0: 2018-06-18, mdevogele@lowell.edu
"""

import argparse, shlex

import SP_Toolbox as SP
import SP_Deveny_Toolbox as SP_Dev_Tool
import numpy as np
import SP_diagnostics as diag

from SP_CheckInstrument import CheckInstrument

from astropy.io import fits

import cosmic

import _SP_conf

def Coadd(filename,OutName):
    
    telescope, obsparam = CheckInstrument([filename[0]])
    for idx, elem in enumerate(filename): 
        hdulist = fits.open(elem)
        data=hdulist[0].data
        if idx == 0:
            print(idx)
            dataTot = data
        else:
            print(idx)
            dataTot += data
    hdulist[0].data = dataTot
    hdulist.writeto(OutName)
        


if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='Coadd fits files')

    parser.add_argument('images', help='images to process or \'all\'',
                        nargs='+')
    parser.add_argument('-o',
                        help='Name of the fits file containing the coadded image',
                        default='MasterFlat.fits')   

    args = parser.parse_args()
    filenames = args.images    
    OutName = args.o
    
    # call run_the_pipeline only on filenames
    Coadd(filenames,OutName)
    pass
