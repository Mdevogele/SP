#!/usr/bin/env python

""" SP_DetectSpectra - Detect spectra for Gemini
    v1.0: 2018-06-18, mdevogele@lowell.edu
"""

import argparse, shlex

import SP_Toolbox as SP
import SP_Deveny_Toolbox as SP_Dev_Tool
import numpy as np
import SP_diagnostics as diag

from SP_CheckInstrument import CheckInstrument

import simplejson

from astropy.io import fits

import cosmic

import _SP_conf

def Detect_Spectra(filename,OutName,DetecLim):
    
    Offset = []
    telescope, obsparam = CheckInstrument([filename[0]])
        
    for idx, elem in enumerate(filename): 
        hdulist = fits.open(elem)
        Bin = int(hdulist[1].header['CCDSUM'][0])
        data=hdulist[0].data
        hdulist = fits.open(elem)
        Offset.append(hdulist[0].header['YOFFSET'])
        if idx == 0:
            print(idx)
            dataTot = data
        else:
            print(idx)
            dataTot += data

    hdulist[0].data = dataTot
    FitsName = OutName + '.fits'
    hdulist.writeto(FitsName)
    
    Sig = float(DetecLim)
    
    Spec = SP.Auto_Detect_Spectra(OutName + '.fits',Sig)
    print(Spec)
    
    PltScale = obsparam['pixscale']
    
    Offset = np.array(Offset)
    if telescope == 'GMOSS':
        Off = -Offset/(PltScale*Bin)
    if telescope == 'GMOSN':
        Off = Offset/(PltScale*Bin)
        
    print(Spec)
    
    Spec_Loc = SP.Get_Spectra(Spec,Off)
    
    if Bin == 4:
        Spec_Loc = Spec_Loc+400
    if Bin == 2:
        Spec_Loc = Spec_Loc+800
    

    SpecLocName = OutName + '_Spec_loc.txt'
    with open(SpecLocName,'w') as f:
        simplejson.dump(list(Spec_Loc), f)

    OffsetName = OutName + '_Offset.txt'
    with open(OffsetName,'w') as f:
        simplejson.dump(list(Offset), f)
    
    print(Spec_Loc)
    print(Off)
        
    
        


if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='Coadd fits files')

    parser.add_argument('images', help='images to process or \'all\'',
                        nargs='+')
    parser.add_argument('-o',
                        help='Name of the fits file containing the coadded image',
                        default='MasterFlat.fits')   
    parser.add_argument('-d',
                        help='Detection limit to detect spectra',
                        default='4')  

    args = parser.parse_args()
    filenames = args.images    
    OutName = args.o 
    DetecLim = args.d
    
    # call run_the_pipeline only on filenames
    Detect_Spectra(filenames,OutName,DetecLim)

    pass
