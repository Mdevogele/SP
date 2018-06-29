#!/usr/bin/env python

""" SP_Bias - wrapper for creating bias 
    v1.0: 2018-04-17, mdevogele@lowell.edu
"""

import argparse, shlex

import SP_Toolbox as SP
import SP_Deveny_Toolbox as SP_Dev_Tool
import numpy as np

from astropy.io import fits


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
#    parser.add_argument('-prefix', help='data prefix',
#                        default=None)
#    parser.add_argument('-target', help='primary targetname override',
#                        default=None)
#    parser.add_argument('-filter', help='filter name override',
#                        default=None)
#    parser.add_argument('-fixed_aprad', help='fixed aperture radius (px)',
#                        default=0)
#    parser.add_argument('-source_tolerance',
#                        help='tolerance on source properties for registration',
#                        choices=['none', 'low', 'medium', 'high'],
#                        default='high')
    parser.add_argument('-v',
                        help='Increase verbosity',
                        action="store_true")    
    parser.add_argument('-o',
                        help='Name of the master bias file',
                        default='MasterBias.fits')    
    parser.add_argument('images', help='images to process or \'all\'',
                        nargs='+')

    args = parser.parse_args()
#    prefix = args.prefix
#    man_targetname = args.target
#    man_filtername = args.filter
#    fixed_aprad = float(args.fixed_aprad)
#    source_tolerance = args.source_tolerance
    Verbose = args.v
    MasterName = args.o
    filenames = args.images    
    
    # call run_the_pipeline only on filenames
    Create_Bias(filenames,MasterName,Verbose)
    pass