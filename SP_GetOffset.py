#!/usr/bin/env python

""" SP_GetOffset - Coadd fits files
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

def GetOffset(filename):
    
    telescope, obsparam = CheckInstrument([filename[0]])
    for idx, elem in enumerate(filename): 
        hdulist = fits.open(elem)
        Offset = hdulist[0].header['YOFFSET']
        print(Offset)
        


if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='Coadd fits files')

    parser.add_argument('images', help='images to process or \'all\'',
                        nargs='+')

    args = parser.parse_args()
    filenames = args.images    
    
    # call run_the_pipeline only on filenames
    GetOffset(filenames)
    pass
