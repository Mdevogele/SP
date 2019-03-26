#!/usr/bin/env python

""" SP_Coadd - Coadd fits files
    v1.0: 2018-06-18, mdevogele@lowell.edu
"""

import argparse, shlex

import SP_Toolbox as SP

from SP_CheckInstrument import CheckInstrument

from astropy.io import fits


def Subtract(filename,OutName,outFile):
    
    telescope, obsparam = CheckInstrument([filename[0]])
    Image1 = filename[0]
    Image2 = SecondFile
    
    hdulist1 = fits.open(Image1)
    data1 = hdulist1[0].data
    
    hdulist2 = fits.open(Image2)
    data2 = hdulist2[0].data    
    
    Result = data1-data2
    
    hdulist1[0].data = Result
    hdulist1.writeto(outFile)
        


if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='Substract two fits file')

    parser.add_argument('images', help='images to process or \'all\'',
                        nargs='+')
    parser.add_argument('-s')   
    parser.add_argument('-o') 

    args = parser.parse_args()
    filenames = args.images    
    SecondFile = args.s
    outFile = args.o
    
    # call run_the_pipeline only on filenames
    Subtract(filenames,SecondFile,outFile)
    pass
