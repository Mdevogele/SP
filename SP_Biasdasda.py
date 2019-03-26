#!/usr/bin/env python

""" SP_Bias - wrapper for creating bias 
    v1.0: 2018-04-17, mdevogele@lowell.edu
"""

import argparse, shlex

import SP_Toolbox as SP
import SP_Deveny_Toolbox as SP_Dev_Tool
import numpy as np

from astropy.io import fits
from SP_CheckInstrument import CheckInstrument


def Create_Bias(filenames,MasterName,Verbose):
    
    if Verbose:
        print('Beginning bias processing')
        print('Processing files:')
        print('index \t filename')
        for idx,elem in enumerate(filenames):
            print('{} \t {}'.format(idx+1,elem))

    BiasFlags = {
    'logfile': 'biasLog.txt','RawPath':'', 'WriteFile': './' + MasterName,
    'verbose': False, 'OverWrite': True, 'AddFits': False, 'IsGMOS': False
    }
  
    if Verbose:
        print('Creating the master bias')
    SP.Create_Bias(filenames,**BiasFlags)    
    if Verbose:
        print('Master bias save to {}'.format(MasterName))
        hdulist = fits.open(MasterName)
        data = hdulist[0].data
        print('Statistics of the Master bias')
        print('Mean: {} \t Median: {} \t std: {}'.format(np.mean(data), np.median(data), np.std(data)))
        print('End of bias processing')
        
    




if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='Processing and creation of master bias')

    parser.add_argument('-v',
                        help='Increase verbosity',
                        action="store_true")    
    parser.add_argument('-o',
                        help='Name of the master bias file',
                        default='MasterBias.fits')    
    parser.add_argument('images', help='images to process or \'all\'',
                        nargs='+')

    args = parser.parse_args()

    Verbose = args.v
    MasterName = args.o
    filenames = args.images    
    fewfdscsddoihwerhoweuriwe
    
    
    print(filenames)
    SP.Check(filenames)

    Create_Bias(filenames,MasterName,Verbose)
    pass